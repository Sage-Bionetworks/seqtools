#!/usr/bin/env python



import os
import sys
import argparse
import functools



##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));
import bamUtils
import utils
import parallel




def build_parser():

    description  =  """
                    ######
                    Script to process BAMs
                    ######
                    """
    
    #user command line args processing
    parser = argparse.ArgumentParser(description=description,add_help=True)
    parser.add_argument('-nc','--num_cores',dest='num_cores',default=1,type=int,metavar='',help='number of slots to use')
    parser.add_argument('-f', '--force_run',dest='force_run',default=False,action='store_true',help='force run even if expected output is present')
    parser.add_argument('-d', '--debug',dest='debug',default=False,action='store_true',help='debug mode activated')
    parser.add_argument('-o', '--output_dir',dest='output_dir',default=None,help='dir to place the output in..<default: same as input file dir>')


    subparsers = parser.add_subparsers( title='commands', description='The following commands are available:', 
                                        help='For additional help: "bam.py <COMMAND> -h"')

    #commands for running bamstats
    parser_mapstats = subparsers.add_parser('mapstats', help='generate the mapping stats for sam/bam file/s')
    parser_mapstats.add_argument('-b','--bams',dest='bam_files',nargs='+',default=None,metavar='',help='list of bam/sam files')
    parser_mapstats.set_defaults(func=bamUtils.get_bamStats)


    return parser



def run(args):
    if 'func' in args:
        try:
            func = functools.partial(args.func,force=args.force_run,output_dir=args.output_dir)
            result = parallel.parallelize_func_single_input(func,args.bam_files,num_cores=args.num_cores)
            return result
        except Exception as ex:
            if args.debug:
                raise
            else:
                print ex
                #sys.stderr.write(ex)


def main():
    SUB = 'main'
    
    args = build_parser().parse_args()
    

    #create output dir if not exists
    if args.output_dir:
        args.output_dir = os.path.abspath(args.output_dir)
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
    
    #run the selected command    
    run(args)

if __name__ == "__main__":
    main()
    








        
            
    
    
    






