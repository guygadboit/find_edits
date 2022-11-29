#!/usr/bin/env python3
import numpy as np
import math
from scipy.stats import fisher_exact
from argparse import ArgumentParser
import translate as tr
from translate import find_codons, get_residue, parse_orfs, alternatives
from textwrap import wrap
from load_genomes import GenomeSet
from utils import array_from_string, patterns
from collections import namedtuple
from pdb import set_trace as brk


class MutationMap:
	def __init__(self, a, b, a_orfs):
		"""a and b are both aligned genomes"""
		self.a = a
		self.b = b
		self.a_orfs = a_orfs

		# Sets of nt offsets for where the silent and non-silent mutations are
		self.silent = set()
		self.non_silent = set()

		# Map the offset of each silently mutated base to the number of other
		# ways that same mutation could be achieved (by modifying the base in
		# either genome)
		self.silent_alternatives = {}

		# The total number of alternatives for each base across the whole
		# genome.
		self.total_alternatives = 0

		# Just counts for these
		self.num_silent_codon_changes = 0
		self.num_non_silent_codon_changes = 0

		self.find_mutations()

	def find_mutations(self):
		for a_codon, (i, j) in tr.find_codons(self.a, self.a_orfs):
			b_codon = self.b[i:j]

			a_residue = tr.get_residue(a_codon)
			b_residue = tr.get_residue(b_codon)

			silent = a_residue == b_residue

			if silent:
				save_in = self.silent
			else:
				save_in = self.non_silent

			codon_changed = False
			for x in range(i, j):
				if self.a[x] != self.b[x]:
					save_in.add(x)

					if silent:
						ax = alternatives(a_residue, x-i) - 1
						bx = alternatives(b_residue, x-i) - 1
						self.silent_alternatives[x] = ax * bx

					codon_changed = True

			if codon_changed:
				if a_residue == b_residue:
					self.num_silent_codon_changes += 1
				else:
					self.num_non_silent_codon_changes += 1

		for v in self.silent_alternatives.values():
			self.total_alternatives += v

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

		print("For each silent mutation the number of alternatives")
		print(self.silent_alternatives)


	@staticmethod
	def _match(genomes, patterns):
		for pat in patterns:
			apat = array_from_string(pat)
			n = len(apat)
			for i in range(genomes.shape[1] - n - 1):
				for j in range(2):
					if np.array_equal(genomes[j][i:i+n], apat):
						yield i, i+n


	def silent_mutations_in_sequences(self, patterns):
		"""patterns is a list of things like "GAGACC". Return the number of
		silent mutations in those patterns, the number outside them and the
		odds ratio described in the paper"""
		genomes = np.vstack((self.a, self.b))
		num = 0
		matches = 0

		# Use for computing the odds ratio
		a, b, c, d = 0.0, 0.0, 0.0, 0.0

		for start, end in self._match(genomes, patterns):
			for k in range(start, end):
				if k in self.silent:
					a += 1
				b += 1

		# b is the number of non-silently or non-mutated nts in the sites
		b -= a

		# Silent mutations everywhere else
		c = len(self.silent) - a

		# Other nucleotides everywhere else
		d = len(self.a) - a - c - b

		contingency_table = np.array([[a, b], [c, d]], dtype=float)
		OR, p = fisher_exact(contingency_table)

		return a, c, OR, p

	def sum_alternatives(self, patterns):
		genomes = np.vstack((self.a, self.b))

		a, b, c, d = 0.0, 0.0, 0.0, 0.0

		inside = 0
		for start, end in self._match(genomes, patterns):
			for k in range(start, end):
				if k in self.silent_alternatives:
					a += 1
					b += self.silent_alternatives[k]

		c = len(self.silent) - a
		d = self.total_alternatives - b

		contingency_table = np.array([[a, b], [c, d]], dtype=float)
		OR, p = fisher_exact(contingency_table)

		return OR, p

	def output_clu(self, name_a, name_b, fp):
		"""Output in a sort of clu format but plus the residues"""
		translator = tr.Translator([self.a, self.b], self.a_orfs)

		name_a = name_a.split()[0][:12]
		name_b = name_b.split()[0][:12]

		def w(*args, **kwargs):
			print(*args, **kwargs, end='', file=fp)

		w("Sequence alignment plus residues\n\n\n")

		for i in range(0, len(self.a), 60):
			end = min(len(self.a), i+60)
			for name, genome in ((name_a, self.a), (name_b, self.b)):

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
			residues_a, residues_b = "", ""
			for j in range(i, end):
				residues = translator.get_residues(j)
				residues_a += residues[0]
				residues_b += residues[1]

			w("{}{}\n".format(name_a.ljust(16), residues_a))
			w("{}{}\n".format(name_b.ljust(16), residues_b))
			w('\n')


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)

	# These should be the ORFs for the first genome of the pair
	ap.add_argument("-r", "--orfs", type=str, default="WH1-orfs")

	ap.add_argument('-o', "--output", type=str, default="residues.clu")
	ap.add_argument('-n', "--residues-only", action="store_true",
			help="Just make the residues file")
	ap.add_argument('-e', "--exhaustive", action="store_true")
	args = ap.parse_args()

	gs = GenomeSet(args.fname[0])
	if gs.num_genomes() != 2:
		print("You're supposed to only have two genomes for this one!")
		return

	mm = MutationMap(gs.genomes[0], gs.genomes[1], parse_orfs(args.orfs))

	with open(args.output, "w") as fp:
		mm.output_clu(gs.get_name(0), gs.get_name(1), fp)

	interesting = ("CGTCTC", "GAGACC", "GGTCTC", "GAGACG")
	if args.residues_only:
		return

# 	mm.summary()

	print("<num_silent in sites> <num_silent outside> <OR> <p>")

	# Consider them all together
	print("All together")
# 	print(*mm.silent_mutations_in_sequences(interesting))
	print("Alternatives", *mm.sum_alternatives(interesting))

	# And one at a time
	for pat in interesting:
# 		print(pat, *mm.silent_mutations_in_sequences((pat,)))
		print("Alternatives", pat, *mm.sum_alternatives((pat,)))

	if not args.exhaustive:
		return

	print("Controls")

	total, count = 0.0, 0.0

	for pat in patterns():
		if pat in interesting: continue
		print(pat, *mm.sum_alternatives((pat,)))

# 		a, c, OR, p = mm.silent_mutations_in_sequences((pat,))
# 		if not math.isnan(OR):
# 			total += OR
# 			count += 1
# 		print(pat, a, c, OR, p)

	print("Average OR where defined:", total / count)


if __name__ == "__main__":
	main()
