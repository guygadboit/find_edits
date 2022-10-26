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
	return CODON_TABLE.get(codon, 'X')
