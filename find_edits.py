import numpy as np
from itertools import product
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


def main():
	genomes = load("input.fasta")
	for p in patterns():
		print(p, score(genomes, p))


if __name__ == "__main__":
	main()
