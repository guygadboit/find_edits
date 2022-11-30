from argparse import *
import re
from pdb import set_trace as brk


# Convert a residues.clu file back into a fasta file
def main():
	ap = ArgumentParser()
	ap.add_argument("-k", "--key", type=str)
	ap.add_argument("-o", "--outname", type=str)
	ap.add_argument("fname", nargs=1)

	args = ap.parse_args()
	outline = ""

	with open(args.outname, "wt") as output:
		print(">{}".format(args.key), file=output)
		with open(args.fname[0]) as fp:
			for line in fp:
				m = re.match(r'([\w\d\.]+)\s+([GACT\-]+)(\s+\d+)$', line)
				if not m or m.group(1) != args.key:
					continue

				data = m.group(2)
				data = re.sub(r'-', '', data)
				outline += data

				if len(outline) >= 70:
					print(outline[:70], file=output)
					outline = outline[70:]

		print(outline, file=output)

	print("Wrote {}".format(args.outname))


if __name__ == "__main__":
	main()
