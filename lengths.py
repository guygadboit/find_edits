#!/usr/bin/env python3
import numpy as np
from argparse import ArgumentParser
from load_genomes import GenomeSet
from utils import patterns, array_from_string
from random import *
from translate import Translator
from pdb import set_trace as brk


INTERESTING = (
		"CGTCTC",
		"GAGACC",
		"GGTCTC",
		"GAGACG",
		)

ALL_PATTERNS = list(patterns(6))


def fragments(genome, patterns):
	apats = [array_from_string(pat) for pat in patterns]
	n = len(apats[0])
	start = 0

	for i in range(0, len(genome)):
		sl = genome[i:i+n]
		for apat in apats:
			if np.array_equal(apat, sl):
				yield start, i
				start = i

	yield start, len(genome)


def longest_fragment(genome, patterns):
	"""Return the number of fragments and the longest one"""
	max_length = 0

	for count, (start, end) in enumerate(fragments(genome, patterns)):
		length = end - start
		max_length = max(max_length, length)

	return count+1, max_length


def random_patterns(n=2):
	return [ALL_PATTERNS[randint(0, len(ALL_PATTERNS)-1)] for _ in range(n)]


def reverse_complement(pattern):
	COMPLEMENTS = {
			'G': 'C',
			'A': 'T',
			'C': 'G',
			'T': 'A'
			}

	ret = [COMPLEMENTS[c] for c in pattern]
	return "".join(reversed(ret))


def montecarlo(genome, n=1000):
	for _ in range(n):
		# To make it a fair comparison two of the patterns will be reverse
		# complements of the other two.
		patterns = random_patterns(2)
		patterns.extend([reverse_complement(p) for p in patterns])
		yield longest_fragment(genome, patterns)


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument('-i', "--iterations", type=int)
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])
	genome = gs.genomes

	print("BsaI/BsmBI sites", *longest_fragment(genome, INTERESTING))

	n = args.iterations
	if n:
		print("Running for {} iterations".format(n))
		with open("montecarlo.txt", "w") as fp:
			for i, lf in enumerate(montecarlo(gs.genomes, n)):
				print(*lf, file=fp)
				if i % 100 == 0:
					print("{}/{} done".format(i, n))


if __name__ == "__main__":
	main()
