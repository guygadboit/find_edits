import numpy as np
from argparse import ArgumentParser
from load_genomes import GenomeSet
from collections import namedtuple
from textwrap import wrap
from copy import *
from pdb import set_trace as brk


CODON_TABLE = {
		"TTT": 'F',	# Phenylalanine
		"TTC": 'F',

		"TTA": 'L',	# Leucine
		"TTG": 'L',
		"CTT": 'L',
		"CTC": 'L',
		"CTA": 'L',
		"CTG": 'L',

		"ATT": 'I',	# Isoleucine
		"ATC": 'I',
		"ATA": 'I',

		"ATG": 'M',	# Methionine

		"GTT": 'V',	# Valine
		"GTC": 'V',
		"GTA": 'V',
		"GTG": 'V',

		"TCT": 'S',	# Serine
		"TCC": 'S',
		"TCA": 'S',
		"TCG": 'S',

		"CCT": 'P',	# Proline
		"CCC": 'P',
		"CCA": 'P',
		"CCG": 'P',

		"ACT": 'T',	# Threonine
		"ACC": 'T',
		"ACA": 'T',
		"ACG": 'T',

		"GCT": 'A',	# Alanine
		"GCC": 'A',
		"GCA": 'A',
		"GCG": 'A',

		"TAT": 'Y',	# Tyrosine
		"TAC": 'Y',

		"TAA": '*',	# Stop
		"TAG": '*',

		"CAT": 'H',	# Histidine
		"CAC": 'H',

		"CAA": 'Q',	# Glutadine
		"CAG": 'Q',

		"AAT": 'N',	# Asparagine
		"AAC": 'N',

		"AAA": 'K',	# Lysine
		"AAG": 'K',

		"GAT": 'D',	# Aspartic acid
		"GAC": 'D',

		"GAA": 'E',	# Glutamic acid
		"GAG": 'E',

		"TGT": 'C',	# Cysteine
		"TGC": 'C',

		"TGA": '*',	# Stop
		"TGG": 'W',	# Tryptophan

		"CGT": 'R',	# Arginine
		"CGC": 'R',
		"CGA": 'R',
		"CGG": 'R',

		"AGT": 'S',	# Serine
		"AGC": 'S',

		"AGA": 'R',	# Arginine (again)
		"AGG": 'R',

		"GGT": 'G',	# Glycine
		"GGC": 'G',
		"GGA": 'G',
		"GGG": 'G',
	}


ORF = namedtuple("ORF", "start end")


def parse_orfs(fname):
	ret = []
	with open(fname) as fp:
		for line in fp:
			start, end = [int(x) for x in line.split()]
			# ORFs seem to be conventionally listed as 1-based
			ret.append(ORF(start-1, end))
	return ret


class ORFIterator:
	def __init__(self, orfs):
		self.orfs = orfs
		self._i = 0
		self.orf_i = 0

		# The genome may contain blanks ('-') as part of its alignment with the
		# others. But our list of ORF indices doesn't take that into account.
		# So we return the real offset plus "_skip", which is the numer of '-'s
		# seen by the thing translating the genome
		self._skip = 0

		self._next_orf()

		# Move to the start of the first ORF
		self.incr()

	def _next_orf(self):
		if self.orf_i < len(self.orfs):
			self.current_orf = self.orfs[self.orf_i]
			self._i = 0
		else:
			self.current_orf = None
		self.orf_i += 1

	def incr(self):
		if self.current_orf is None:
			raise StopIteration

		if self._i < self.current_orf.start:
			self._i = self.current_orf.start
			return self._i + self._skip

		self._i += 1
		if self._i == self.current_orf.end:
			self._next_orf()
			return self.incr()

	def add(self, n):
		for _ in range(n):
			self.incr()

	def set(self, n):
		assert n >= self._i
		while self._i != n:
			self.incr()

	def skip(self):
		self._skip += 1

	def val(self):
		"""Where to the read actual genome from"""
		return self._i + self._skip


def get_residue(codon):
	"""Get the residue of a codon, which can either be a string like 'GGT' or
	an np array"""
	if isinstance(codon, np.ndarray):
		codon = "".join([chr(x) for x in codon])
	return CODON_TABLE.get(codon, 'X')


def next_codon(genome, i):
	"""Find the next codon (skipping gaps) in the genome starting from i.
	Return the actual codon and the position of the end of it"""
	t = ''

	while True:
		c = chr(genome[i.val()])

		if c == '-':
			i.skip()
		else:
			t += c

		if len(t) == 3:
			ret = t, i.val()
			i.incr()
			return ret

		i.incr()


def find_codons(genome, orfs):
	"""Generate the offsets of each codon"""
	i, n = ORFIterator(orfs), len(genome)

	while True:
		start = i.val()
		codon, end = next_codon(genome, i)
		yield codon, (start, end+1)


class Translator:
	def __init__(self, genome, orfs):
		self.genome = genome
		self.translator = find_codons(genome, orfs)
		self.current_residue = None
		self.pos = (0, 0)

	def get_residue(self, offset):
		"""Intended to be used as a coroutine-- you can only call this with
		monotonically increasing offsets, i.e. from a for loop over all the
		offsets basically"""
		if self.pos[0] <= offset < self.pos[1]:
			return self.current_residue

		if offset >= self.pos[1]:
			try:
				codon, self.pos = next(self.translator)
				self.current_residue = get_residue(codon)
				return self.get_residue(offset)
			except StopIteration:
				pass

		return '|'	# We'll use this to mean untranslated


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument("-r", "--orfs", type=str, default="WH1-orfs")
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])

	start_pos = None
	orfs = []
	line = ""
	debugging = False

	orfs = parse_orfs(args.orfs)

	for codon, pos in find_codons(gs.genomes, orfs):
		residue = get_residue(codon)
		if residue == '*':
			continue
		line += residue

		if len(line) == 70:
			print(line)
			line = ""

	print(line)


if __name__ == "__main__":
	main()
