import numpy as np
from argparse import ArgumentParser
from load_genomes import GenomeSet
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


def get_residue(codon):
	"""Get the residue of a codon, which can either be a string like 'GGT' or
	an np array"""
	if isinstance(codon, np.ndarray):
		codon = "".join([chr(x) for x in codon])
	return CODON_TABLE.get(codon, 'X')


def next_codon(genome, i):
	"""Find the next codon (skipping gaps) in the starting from i. Return the
	actual codon and the position of the end of it"""
	t = ''
	while i < len(genome):
		c = chr(genome[i])

		if c != '-':
			t += c

		if len(t) == 3:
			return t, i

		i += 1
	raise StopIteration


def find_codons(genome):
	"""Generate the offsets of each codon"""
	i, n = 0, len(genome)

	count = 0
	start = 0

	reading = False

	while True:
		if not reading:
			# Skip one nt at a time until we find the start codon
			codon, j = next_codon(genome, i)
			r = get_residue(codon)
			i += 1

			if r == 'M':
				yield codon, (i-1, j+1)
				reading = True
				i += 2	# skip to the end of the M codon itself
				continue
		else:	# We are reading
			codon, j = next_codon(genome, i)
			yield codon, (i, j+1)
			i = j + 1

			if codon == '*':
				reading = False

class Translator:
	def __init__(self, genome):
		self.genome = genome
		self.translator = find_codons(genome)
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
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])

	start_pos = None
	orfs = []

	for codon, pos in find_codons(gs.genomes):
		residue = get_residue(codon)

		if residue == 'M':
			if not start_pos:
				start_pos = pos[0]
				print("Start at", start_pos)
		elif residue == '*':
			if start_pos:
				orfs.append((start_pos, pos[1]))
				print("\nEnd at", pos[1])
				start_pos = None

		if start_pos:
			print(residue, end="")

	print("\nORFs")
	for orf in orfs:
		print(*orf)


if __name__ == "__main__":
	main()
