#!/usr/bin/env python3
import numpy as np
from argparse import ArgumentParser
from load_genomes import GenomeSet
from utils import patterns, array_from_string
from pdb import set_trace as brk


def locate(genomes, pattern):
	"""Generate all the slices containing pattern in any of the genomes.
	pattern is something like 'GAGACC'"""
	apat = array_from_string(pattern)

	for s in slices(genomes, len(apat)):
		for row in s:
			if np.array_equal(row, apat):
				yield s


def slices(genomes, length):
	"""Generate all vertical slices length nucleotides long."""
	for i in range(genomes.shape[1] - length - 1):
		yield genomes[:, i:i+length]


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

def score_mutated(genomes, pattern):
	"""Alternative score function for a different analysis: score 1 if there is
	a 1nt difference for this pattern"""
	count = 0

	for sl in locate(genomes, pattern):
		for row in sl[1:]:
			count = 0
			for a, b in zip(sl[0], row):
				if a == b:
					count += 1

			# The same in 5 out of the 6 will count as a yes
			if count == 5:
				return 1

	return 0


def get_triplet(row, i):
	"""Find the next triplet (skipping gaps) in the first row starting from i.
	Return the triplet and the position of the end of it"""
	t = ''
	while i < len(row):
		c = chr(row[i])

		if c != '-':
			t += c

		if len(t) == 3:
			return t, i

		i += 1
	raise StopIteration


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])

	total, count = 0, 0

	for p in patterns():
		# print(p, score(gs.genomes, p))
		score = score_mutated(gs.genomes, p)
		print(p, score)
		total += score
		count += 1

	print("Average score {}".format(float(total) / count))


if __name__ == "__main__":
	main()
