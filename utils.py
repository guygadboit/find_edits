import numpy as np
from itertools import product


def array_from_string(s):
	"""Make a 1D array from a string like 'AGGTAC'"""
	return np.array([ord(c) for c in s], dtype='int8')


def slice_to_string(sl):
	"""Convert a slice back to a string. We just use this for debugging"""
	lines = []
	for row in sl:
		lines.append("".join([chr(c) for c in row]))
	return "\n".join(lines)


def patterns(n=6):
	"""generate all possible patterns of n nucleotides"""
	iterables = []
	for i in range(n):
		iterables.append([b for b in "GATC"])

	for pat in product(*iterables):
		yield "".join(pat)
