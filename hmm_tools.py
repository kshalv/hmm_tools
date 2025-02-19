#!/usr/bin/env python

import helper_functions
import argparse

def main():
	# create main parser
	parser = argparse.ArgumentParser(description='This program takes an input folder of HMMs & loops through each genome in a directory of prokka annotated genomes \
		generating a score thresholded count of each hmm for each genome. Several intermediate outputs are produced.') 

	# create sub-parsers 
	subparsers = parser.add_subparsers(dest='functions', help='Available functions')

	# create sub parser for hmm searching 
	hmm_parser = subparsers.add_parser('hmm', help='Run for HMM searching.')
	hmm_parser.add_argument('-g', '--genomes', help='Path to genome directory')
	hmm_parser.add_argument('-hm', '--hmms', help='Path to hmm directory. All files should have a .HMM extension and named according to gene')
	hmm_parser.add_argument('-e', '--evalue', help='Specify e-value for parsing hmm output.')
	hmm_parser.add_argument('-sc', '--score', choices=['tc', 'nc'], help='Specify score for parsing hmm output.')

	# create sub parser for sequence pulling
	seq_parser = subparsers.add_parser('seq', help='Run for sequence pulling from hmm hits.')
	seq_parser.add_argument('-hi', '--hits', help='Path to hmm_hits. This is an intermediate directory produced by hmm searcher.' )
	seq_parser.add_argument('-hm', '--hmms', help='Path to hmm directory')
	seq_parser.add_argument('-g', '--genomes', help='Path to genome directory')
	seq_parser.add_argument('-s', '--sequence', choices=['nt','aa'], help='Specify sequence type. Nucleotide (-nt) or amino acid (-aa)')
	seq_parser.add_argument('-f', '--filter', choices=['all', 'tc', 'nc'], help='Specify filtering threshold. "All" removes score threshold, whereas "tc" & "nc" imposes the TC and NC cutoff from the seed HMM.')

	args = parser.parse_args()

	# run code on inputs and outputs
	# run code on inputs and outputs
	if args.functions == 'hmm':
		out_path = helper_functions.hmmer_run(args.genomes, args.hmms)
		if not args.score:
			print('Error: please specify the --score option.')
		else:
			if args.evalue is not None:
				hit_path = helper_functions.hmmer_parser(out_path, args.hmms, args.evalue)
			else:
				hit_path = helper_functions.hmmer_parser(out_path, args.hmms, 1e-03)
			helper_functions.filt_count(args.hmms, hit_path, args.score)
	elif args.functions == 'seq':
		if not args.filter:
			print('Error: please specify the --filter option.')
		helper_functions.seq_puller(args.hits, args.hmms, args.genomes, args.filter, args.sequence)
	else:
		print("Invalid function. Use 'hmm' or 'seq'")

if __name__ == '__main__':
	main()

