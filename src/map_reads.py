#!/usr/bin/env python

import os
import argparse
import sys
import subprocess
import time



##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));
import utils
import mappers



def setPATH(config_file=None):

	#THIS IS where the path to tophat executables should be defined
	ORIG_PATH = os.environ['PATH']
	TOPHAT_DIR="/home/apratap/softwares/tophat-2.0.11.Linux_x86_64/"
	BOWTIE_DIR="/home/apratap/softwares/bowtie2-2.2.1/"
	STAR_DIR="/home/apratap/softwares/STAR_2.3.0e/"	
	NEW_PATH=[TOPHAT_DIR,BOWTIE_DIR,STAR_DIR,ORIG_PATH]
	os.environ['PATH'] =':'.join(NEW_PATH)
	


def build_parser():
	parent_parser = argparse.ArgumentParser(add_help=False)
	parent_parser.add_argument('--gi', '--genome_index', required=True, dest='genome_index', help='path to genome index')
	parent_parser.add_argument('--ref_genes',dest='ref_genes', help='path to reference gene set')
	parent_parser.add_argument('--cores', '-c',default=1, type=str, dest='cores',help='number of cores to use')
	parent_parser.add_argument('-o', '--outdir', dest='outdir',default = None, help='output dir (defaults to fastq base dir')
	parent_parser.add_argument('-p', '--prefix',dest='prefix',default = None,help='output dir (defaults to fastq base dir named)')
	parent_parser.add_argument('--read1', '--r1', required=True, dest='read1', nargs = '+', help='read 1 fastq')
	parent_parser.add_argument('--read2', '--r2', dest='read2', default = list(), nargs = '+',help='read 2 fastq')
	parent_parser.add_argument('--optional_args', "--optargs",type=str, dest='optional_args', 
							   help='space separated optional args to the mapper')
	parent_parser.add_argument('--force', "-f", action="store_true", dest='force', default=False, help='force_run')
	parent_parser.add_argument('--debug', "-d", action="store_true", dest='debug', default=False, help='debug')

	parser = argparse.ArgumentParser(add_help=False)
	subparsers = parser.add_subparsers( title='commands', description='The following commands are available:', 
										help='For additional help: "%s <COMMAND> -h"' % os.path.basename(__file__))

	#running tophat
	parser_tophat = subparsers.add_parser('tophat', help='run tophat2 alignment', parents=[parent_parser])
	parser_tophat.set_defaults(func=mappers.run_tophat2)

	#running STAR
	parser_STAR = subparsers.add_parser('star', help='run star alignment', parents=[parent_parser])
	parser_STAR.set_defaults(func=mappers.run_STAR)

	return parser


def run(args):

	#set the paths to the dir
	setPATH()

	if 'func' in args:
		try:
			func = args.func
			func(args)
		except Exception as ex:
			if args.debug:
				raise
			else:
				sys.stderr.write(str(ex))


def __fix_dirs_N_initialize_logger(args):
	#outdir 
	if args.outdir is None:
		args.outdir = os.path.abspath(os.path.dirname(args.read1[0]))
	if args.prefix is None:
		args.prefix = '%s_%s' % (os.path.basename(__file__).replace('.py',''), time.strftime('%Y-%m-%d_%H-%M-%S'))
	else:
		args.prefix = '%s_%s' % (args.prefix,time.strftime('%Y-%m-%d_%H-%M-%S'))
	
	#get the logger
	log_file_name = args.outdir + '/' + args.prefix + '.log'
	args.logger = utils._get_logger(log_file_name)



if __name__ == "__main__":

	#get the args
	args = build_parser().parse_args()

	#internal func
	__fix_dirs_N_initialize_logger(args)

	#log the user arguments
	args.logger.info('###########################')
	args.logger.info("Command line options used ")
	[ args.logger.info('%s=%s' %(arg_name,value)) for arg_name,value in args._get_kwargs()]
	args.logger.info('###########################')

	#map the reads
	run(args)


  
