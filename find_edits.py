import numpy as np
from pdb import set_trace as brk


def lines(fname):
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			yield line


def load(fname):
	"""Load everything into a big numpy array"""

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
			l = np.array([ord(c) for c in line], dtype='int8')

			if row is None:
				row = l
			else:
				row = np.hstack((row, l))

	return ret


def main():
	genomes = load("input.fasta")
	brk()


if __name__ == "__main__":
	main()




