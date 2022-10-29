import numpy as np
from utils import array_from_string
from pdb import set_trace as brk


def _lines(fname):
	with open(fname) as fp:
		for line in fp:
			line = line.strip()
			yield line


class GenomeSet:
	def __init__(self, fname):
		"""Load a bunch of genomes from a fasta file into a numpy matrix, and
		also keep al ist of their names"""
		self.names = []
		self.genomes = None
		self._load(fname)

	def _load(self, fname):
		"""Load everything into a big numpy array. One genome per row, one
		nucleotide per column"""
		row, genomes = None, None

		for line in _lines(fname):
			if line.startswith('>'):
				self.names.append(line[1:])
				if row is not None:
					genomes = (np.vstack((genomes, row)) if genomes is not None
								else row)
					row = None
			else:
				l = array_from_string(line)

				if row is None:
					row = l
				else:
					row = np.hstack((row, l))

		if row is not None:
			genomes = np.vstack((genomes, row)) if genomes is not None else row

		self.genomes = genomes

	def output_clu(self, mutation_map):
		"""Output a sort of clu format file, but with the residues also
		marked"""

	def num_genomes(self):
		return self.genomes.shape[0]

	def get_name(self, row_num):
		return self.names[row_num]
