#!/usr/bin/env python3
import numpy as np
from itertools import product
from argparse import ArgumentParser
from translate import get_residue
from pdb import set_trace as brk


def array_from_string(s):
	"""Make a 1D array from a string like 'AGGTAC'"""
	return np.array([ord(c) for c in s], dtype='int8')


def slice_to_string(sl):
	"""Convert a slice back to a string. We just use this for debugging"""
	lines = []
	for row in sl:
		lines.append("".join([chr(c) for c in row]))
	return "\n".join(lines)


def lines(fname):
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			yield line


def load(fname):
	"""Load everything into a big numpy array. One genome per row, one
	nucleotide per column"""
	row, ret = None, None

	for line in lines(fname):
		if line.startswith('>'):
			if row is not None:
				if ret is None:
					ret = row
				else:
					ret = np.vstack((ret, row))

				row = None
		else:
			l = array_from_string(line)

			if row is None:
				row = l
			else:
				row = np.hstack((row, l))

	return ret


def slices(genomes, length):
	"""Generate all vertical slices length nucleotides long."""
	for i in range(genomes.shape[1] - length - 1):
		yield genomes[:, i:i+length]


def locate(genomes, pattern):
	"""Generate all the slices containing pattern in any of the genomes.
	pattern is something like 'GAGACC'"""
	apat = array_from_string(pattern)

	for s in slices(genomes, len(apat)):
		for row in s:
			if np.array_equal(row, apat):
				yield s


def differences(sl):
	"""How many of the subsequent rows does the first one differ from?"""
	ret = 0
	for row in sl[1:]:
		if not np.array_equal(sl[0], row):
			ret += 1
	return ret


def score(genomes, pattern):
	"""The average number of differences between the first row and the
	subsequent ones for each time this pattern appears"""
	ret, count = 0, 0

	for sl in locate(genomes, pattern):
		ret += differences(sl)
		count += 1

	if count:
		return float(ret) / count
	else:
		return 0.0


def patterns(n=6):
	"""generate all possible patterns of n nucleotides"""
	iterables = []
	for i in range(n):
		iterables.append([b for b in "GATC"])

	for pat in product(*iterables):
		yield "".join(pat)


def codon_slices(genomes):
	"""Generate slices corresponding to codons in the first sequence"""
	i, n = 0, genomes.shape[1]

	count = 0
	start = 0

	while True:
		nt = chr(genomes[0][i])
		i += 1

		if nt != '-':
			count += 1

		if count == 3:
			count = 0
			yield genomes[:, start:i]
			start = i

		if i == n: break


def translate_codon(nts):
	codon = [chr(x) for x in nts]
	codon = "".join([x for x in codon if x != '-'])
	return get_residue(codon)


def eliminate_non_silent(genomes):
	"""Remove any mutations that aren't silent"""
	for codon_sl in codon_slices(genomes):
		reference = translate_codon(codon_sl[0])

		for i in range(1, genomes.shape[1]):
			residue = translate_codon(genomes[i, :])
			if residue != reference:
				# Non-silent mutation. So pretend it isn't a mutation at all
				codon_sl[i, :] = codon_sl[0, :]


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	args = ap.parse_args()

	genomes = load(args.fname[0])
	eliminate_non_silent(genomes)
	return

	for p in patterns():
		print(p, score(genomes, p))


if __name__ == "__main__":
	main()
