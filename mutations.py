#!/usr/bin/env python3
import numpy as np
from argparse import ArgumentParser
import translate as tr
from translate import find_codons, get_residue
from textwrap import wrap
from load_genomes import GenomeSet
from utils import array_from_string, patterns
from pdb import set_trace as brk


class MutationMap:
	def __init__(self, a, b):
		"""a and b are both aligned genomes"""
		self.a = a
		self.b = b

		# Sets of nt offsets for where the silent and non-silent mutations are
		self.silent = set()
		self.non_silent = set()

		# Just counts for these
		self.num_silent_codon_changes = 0
		self.num_non_silent_codon_changes = 0

		self.find_mutations()

	def find_mutations(self):
		for a_codon, (i, j) in tr.find_codons(self.a):
			b_codon = self.b[i:j]

			a_residue = tr.get_residue(a_codon)
			b_residue = tr.get_residue(b_codon)

			if a_residue == b_residue:
				save_in = self.silent
			else:
				save_in = self.non_silent

			codon_changed = False
			for x in range(i, j):
				if self.a[x] != self.b[x]:
					save_in.add(x)
					codon_changed = True

			if codon_changed:
				if a_residue == b_residue:
					self.num_silent_codon_changes += 1
				else:
					self.num_non_silent_codon_changes += 1

	def _summarize_set(self, s):
		ret = ", ".join([str(x) for x in s])
		return "\n".join(wrap(ret))

	def summary(self):
		print("{} silent mutations".format(len(self.silent)))
		print(self._summarize_set(self.silent))

		print("{} non-silent mutations".format(len(self.non_silent)))
		print(self._summarize_set(self.non_silent))

		print("{} silent codon changes, {} non-silent ones".format(
			self.num_silent_codon_changes,
			self.num_non_silent_codon_changes))

	def silent_mutations_in_sequences(self, patterns):
		"""patterns is a list of things like "GAGACC". Return the number of
		silent mutations in those patterns, the number outside them and the
		odds ratio described in the paper"""
		genomes = np.vstack((self.a, self.b))
		num = 0
		matches = 0

		# Use for computing the odds ratio
		a, b, c, d = 0.0, 0.0, 0.0, 0.0

		for pat in patterns:
			apat = array_from_string(pat)
			n = len(apat)
			for i in range(genomes.shape[1] - n - 1):
				for j in range(2):
					if np.array_equal(genomes[j][i:i+n], apat):
						for k in range(i, i+n):
							if k in self.silent:
								a += 1
							else:
								b += 1

		c = len(self.silent) - a
		d = len(self.non_silent) - b

		if a == 0.0 or c == 0.0:
			OR = -1
		else:
			OR = (a / b) / (c / d)

		return a, b, c, d, OR

	def output_clu(self, name_a, name_b, fp):
		"""Output in a sort of clu format but plus the residues"""
		ta, tb = tr.Translator(self.a), tr.Translator(self.b)

		name_a = name_a.split()[0][:12]
		name_b = name_b.split()[0][:12]

		def w(*args, **kwargs):
			print(*args, **kwargs, end='', file=fp)

		w("Sequence alignment plus residues\n\n\n")

		for i in range(0, len(self.a), 60):
			end = min(len(self.a), i+60)
			for name, genome in (
					(name_a, self.a),
					(name_b, self.b)):

				w(name.ljust(16))
				for j in range(i, end):
					w(chr(genome[j]))

				w("\t{}\n".format(j+1))

			# Now where they differ
			w(''.ljust(16))
			for an, bn in zip(self.a[i:end], self.b[i:end]):
				c = '*' if an == bn else ' '
				w(c)
			w('\n')

			# And now the residues
			for name, genome, translator in (
					(name_a, self.a, ta),
					(name_b, self.b, tb)):
				w(name.ljust(16))

				for j in range(i, end):
					residue = translator.get_residue(j)
					w(residue)
				w('\n')
			w('\n')


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument('-o', "--output", type=str, default="residues.clu")
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])
	if gs.num_genomes() != 2:
		print("You're supposed to only have two genomes for this one!")
		return

	mm = MutationMap(gs.genomes[0], gs.genomes[1])
	mm.summary()

	with open(args.output, "w") as fp:
		mm.output_clu(gs.get_name(0), gs.get_name(1), fp)

	interesting = ("CGTCTC", "GAGACC", "GGTCTC", "GAGACG")

	for pat in interesting:
		print(pat, *mm.silent_mutations_in_sequences((pat,)))

	print("Controls")
	for pat in patterns():
		if pat in interesting: continue
		print(pat, *mm.silent_mutations_in_sequences((pat,)))


if __name__ == "__main__":
	main()
