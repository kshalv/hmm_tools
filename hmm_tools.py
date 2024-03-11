#!/usr/bin/env python

import helper_functions
import subprocess 
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
	hmm_parser.add_argument('-f', '--filt', choices=['all', 'threshold'], help='Apply filtering to hmm counting. "All" will count all genomes (even those without hits),\
		"filt" will only record genomes with at least one hit beyond the hmm score threshold.')

	# create sub parser for sequence pulling
	seq_parser = subparsers.add_parser('seq', help='Run for sequence pulling from hmm hits.')
	seq_parser.add_argument('-hi', '--hits', help='Path to hmm_hits. This is an intermediate directory produced by hmm searcher.' )
	seq_parser.add_argument('-hm', '--hmms', help='Path to hmm directory')
	seq_parser.add_argument('-g', '--genomes', help='Path to genome directory')

	args = parser.parse_args()

	# run code on inputs and outputs
	if args.functions == 'hmm':
		if not args.filt: 
			print('Error: please specify the --filt option.')
		out_path = helper_functions.hmmer_run(args.genomes, args.hmms)
		hit_path = helper_functions.hmmer_parser(out_path, args.hmms)
		helper_functions.filt_count(args.hmms, hit_path, args.filt)
	elif args.functions == 'seq':
		helper_functions.seq_puller(args.hits, args.hmms, args.genomes)
	else: 
		print("Invalid function. Use 'hmm' or 'seq'")

if __name__ == '__main__':
	main()

