#!/usr/bin/env python3
import numpy as np
from argparse import ArgumentParser
from load_genomes import GenomeSet
import translate as tr
from utils import patterns, array_from_string
from copy import *
from pdb import set_trace as brk


def translate(patch):
	"""Assume patch is codon aligned. Return a residue for each nt e.g.
	RRRLLL"""
	patch = "".join([chr(c) for c in patch])
	assert len(patch) % 3 == 0
	ret = []

	for i in range(0, len(patch), 3):
		c = tr.CODON_TABLE[patch[i:i+3]]
		for j in range(3):
			ret.append(c)

	return ret


def get_codon_offset(orfs, offset):
	for orf in orfs:
		if orf.start <= offset < orf.end:
			return (offset - orf.start) % 3


def search(genome, orfs, pattern):
	"""Return the offset of each potential occurence of pattern and whether it
	is actually there or not (i.e. False means you could have it there with
	only silent mutations)"""
	apat = array_from_string(pattern)
	n = len(apat)

	for i in range(0, len(genome)):
		start, end = i, i+n

		# We actually found the pattern here
		if np.array_equal(apat, genome[start:end]):
			yield i, True
			continue

		codon_offset = get_codon_offset(orfs, i)

		if codon_offset is None:	# Not in an ORF
			continue

		patch = genome[start-codon_offset:end+3-codon_offset].copy()
		residues = translate(patch)[codon_offset:codon_offset+6]

		# OK now try inserting the pattern here
		patch[codon_offset:codon_offset+6] = apat
		mutated_residues = translate(patch)[codon_offset:codon_offset+6]

		# Is it silent?
		if residues == mutated_residues:
			yield i, False


def summary(genome, orfs, pattern):
	potential, actual = 0, 0
	for _, there in search(genome, orfs, pattern):
		if there:
			actual += 1
		else:
			potential += 1

	return actual, potential


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument("-r", "--orfs", type=str, default="WH1-orfs")
	ap.add_argument('-e', "--exhaustive", action="store_true")
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])
	genome = gs.genomes

	orfs = tr.parse_orfs(args.orfs)

	interesting = ("CGTCTC", "GAGACC", "GGTCTC", "GAGACG")

	for pat in interesting:
		print(pat, *summary(genome, orfs, pat))

	if not args.exhaustive:
		return

	print("Controls")
	for pat in patterns():
		if pat in interesting: continue
		print(pat, *summary(genome, orfs, pat))


if __name__ == "__main__":
	main()
