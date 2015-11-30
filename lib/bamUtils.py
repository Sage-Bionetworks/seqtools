#!/usr/bin/env python


import pysam
import os
import time
import numpy as np
import math
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import itertools
import collections
import pandas
import logging
import pybedtools
import pickle
import tempfile
import time
import re
import yaml



##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));
import utils
import bio_format_parsers



def index_bam(bamFile):
    """
    Taking a bamFile as input
    produce a index file for the bam (bamFile.fai)
    """
    SUB = 'index_bam'
    bamFile_bai = bamFile + '.bai'
    if os.path.exists(bamFile_bai):
        print '[%s]: Bam file %s already indexed ' % (SUB,bamFile)
    else:
        print '[%s]: indexing %s' % (SUB,bamFile)
        args = ['samtools','index',bamFile]
        return_code = subprocess.check_call(args) 
        if return_code == 0:
            print '[%s]: Created index for %s' % (SUB,os.path.basename(bamFile))


def name_sort_bam(bamFile,parallel=False, threads = 2,out_q=None):
    """
    Taking a bamFile as input
    produce a sorted bam file for the bam (bamFile)
    """
    SUB = 'name_sort_bam'
    bamFile_name_sorted_prefix = bamFile.replace('.bam','_namesorted')
    bamFile_name_sorted        = bamFile_name_sorted_prefix + '.bam'
    if os.path.exists(bamFile_name_sorted):
        print '[%s]: Bam file %s already name sorted ' % (SUB,bamFile)
        if out_q:
            out_q.put(bamFile_name_sorted)
        return bamFile_name_sorted
    else:
        print '[%s]: Name sorting bam file %s' % (SUB,bamFile)
        if parallel:
            args = ['sambamba','sort','-n','-t',str(threads),'-o',bamFile_name_sorted,bamFile]
        else:
            args = ['samtools','sort','-n',bamFile,bamFile_name_sorted_prefix]

        return_code = subprocess.check_call(args)  
        if return_code == 0:
            print '[%s]: Created name sorted bam for %s' % (SUB,os.path.basename(bamFile))
            if out_q:
                out_q.put(bamFile_name_sorted)
            return bamFile_name_sorted
        else:
            print '[%s]: Error creating name sorted bam for %s' % (SUB,os.path.basename(bamFile))
            if out_q:
                out_q.put(False)
            
    
def coordinate_sort_bam(bamFile):
    """
    Taking a bamFile as input
    produce a coordinate sorted bam file
    """
    SUB = 'coordinate_sort_bam'
    bamFile_coordinate_sorted_preifx = bamFile.replace('.bam','_cord_sorted')
    bamFile_coordinate_sorted        = bamFile_coordinate_sorted_preifx + '.bam'
    if os.path.exists(bamFile_coordinate_sorted):
        print '[%s]: Bam file %s already coordinate sorted ' % (SUB,bamFile)
        return bamFile_coordinate_sorted
    else:
        args = ['samtools','sort',bamFile,bamFile_coordinate_sorted_preifx]
        return_code = subprocess.check_call(args) 
        if return_code == 0:
            print '[%s]: Created coordinate sorted bam for %s' % (SUB,os.path.basename(bamFile))
            return bamFile_coordinate_sorted


def remove_chr(bamFile,delete_original=False):
    """
    For bam files to be biodalliance friendly, chromosome numbers MUST only
    be numbers. (eg. chr1 ---> 1) This function removes the chr
    """
    header = bamFile.replace('.bam','_header.sam')
    bamFile_sam = bamFile.replace('.bam','_nochr.bam')
    #Get header of bam file
    args = ['samtools','view','-H',bamFile, '-o', header]
    return_code = subprocess.check_call(args) 
    if return_code == 0:
        print '[%s]: Extracted header for %s' % (header,os.path.basename(bamFile))
    #Get rid of all the chr in header file
    noChr = ['sed', '-i', 's/chr//g',header]
    nochr_return = subprocess.check_call(noChr) 
    if nochr_return == 0:
        print '[%s]: Removed "chr"' % (header)
    #Rehead bam file so no "chr"
    rehead = ['samtools', 'reheader', header, bamFile,'-o' bamFile_sam]
    rehead_return = subprocess.check_call(rehead) 
    if rehead_return == 0:
        print '[%s]: Created bam file without "chr"' % (bamFile_sam)
    os.remove(header)
    #If the original bam is not needed, delete original and rename new bam to the original name
    if delete_original:
        os.remove(bamFile)
        os.rename(bamFile_sam,bamFile)
        bamFile_sam = bamFile
        print 'Deleted original (%s) and renamed %s to %s' % (bamFile, bamFile_sam, bamFile)
    return bamFile_sam

def get_percent_duplication(bamFile, picardPath = "/opt/picard"):
    bamFile = coordinate_sort_bam(bamFile)
    metrics = bamFile.replace('.bam','_metrics.txt')
    output_duplicates = bamFile.replace('.bam','_duplicates.bam')
    args = ['java', '-jar', os.path.join(picardPath,'dist/picard.jar'), 'MarkDuplicates', 'INPUT=', bamFile,'OUTPUT=',output_duplicates,'METRICS_FILE=', metrics]
    return_code = subprocess.check_call(args) 
    if return_code == 0:
        dups = pandas.read_table(metrics,skiprows=6)
        percent = dups['PERCENT_DUPLICATION']
        #os.remove(metrics)
        os.remove(output_duplicates)
        print percent
        return metrics
    else:
        print "picard MarkDuplicates failed"



def get_library_complexity(bamFile, picardPath = "/opt/picard"):
    bamFile = coordinate_sort_bam(bamFile)
    metrics = bamFile.replace('.bam','_compmetrics.txt')
    args = ['java', '-jar', os.path.join(picardPath,'dist/picard.jar'), 'EstimateLibraryComplexity', 'INPUT=', bamFile,'OUTPUT=',metrics]
    return_code = subprocess.check_call(args) 
    if return_code == 0:
        dups = pandas.read_table(metrics,skiprows=6)
        percent = dups['PERCENT_DUPLICATION'][0]
        #os.remove(metrics)
        print percent
        return metrics
    else:
        print "picard EstimateLibraryComplexity failed"



def get_coverage_of_a_bam_file(bamFile,chr=None):
    '''
    given a bamFile and optional chr name
    the method will return a dict with chr names as keys and 
    coverage at each chr bp as a numpy one dimensional array 
    max_depth in pileup used is 1,000,000
    '''
    
    SUB = 'get_coverage_of_a_bam_file'
    
    bamFile_bai = bamFile + '.bai'
    if not os.path.exists(bamFile_bai):
        bamUtils.index_bam(bamFile)
    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    bamFile_handle = pysam.Samfile(bamFile,'rb')
    if chr is None:
        references = bamFile_handle.references
    else:
        print 'got chr from user %s' % chr 
        references = [chr]

    #storing the lengths of each chromosome
    chr_lengths = {}
    for chromosome,chr_len in itertools.izip(bamFile_handle.references,bamFile_handle.lengths):
        chr_lengths[chromosome] = chr_len
    
    #dict to store the coverage by chr
    cov_by_chr = {};
    
    for chromosome in references:
        chr_len = chr_lengths.get(chromosome,None)
        if chr_len is None:
            print '%s chromsome had no value for length ' % chromosome
        print 'Chr: %s \t %d' %(chromosome,chr_len)
        #inintialize a numpy array of len == length of chromosome with cov = 0
        chr_coverage = np.zeros(chr_len)
        print '[%s]: initialized numpy arry for storing coverage for chr %s of len %d' % (SUB,chromosome,len(chr_coverage))
        for pileup in bamFile_handle.pileup(chromosome,max_depth=1000000):
            #print '%d \t %d' %(pileup.pos,pileup.n)
            #print 'Before %d \t %d' %(pileup.pos,chr_coverage[pileup.pos])
            chr_coverage[pileup.pos] = pileup.n
            #print 'After %d \t %d' %(pileup.pos,chr_coverage[pileup.pos])
        cov_by_chr[chromosome] = chr_coverage;
        
    return cov_by_chr
    
        
def sam_to_bam(samFile,fai_file=None,delete_sam=False):
    """
    Take a samfile and convert into bam file
    """
    SUB = 'sam_to_bam'
    
    if fai_file is None:
        sam_fh  = get_mappedFile_FH(samFile)
        temp_genome_file = tempfile.NamedTemporaryFile(mode='w+t',suffix='.tmp', delete=False) #when delete is True file is deleted as soon as it closed
        [ temp_genome_file.write('%s\t%s\n' % (chr,len)) for chr,len in itertools.izip(sam_fh.references,sam_fh.lengths) ]
        temp_genome_file.close()
        fai_file = temp_genome_file.name
        
    bamFile = samFile.replace('.sam','.bam')
    if os.path.exists(bamFile):
        print '[%s]: Bam file %s already present for sam file %s  ' % (SUB,bamFile,samFile)
        return bamFile
    else:
        args = ['samtools', 'view', '-bS', '-t' ,fai_file, samFile, '-o',  bamFile ]
        #print args
        return_code = subprocess.check_call(args,) 
        if return_code == 0:
            print '[%s]: Created name sorted bam for %s' % (SUB,os.path.basename(bamFile))
            if delete_sam:
                os.remove(samFile)
            return bamFile

    
def bam_to_bedGraph(bamFile, scale=False, sorted=False, **kwargs):
    SUB='bam_to_bedGraph'
    if sorted is not True:
        #sort the bam file by coordinates
        bamFile_sorted = coordinate_sort_bam(bamFile)

    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    bam_prefix = get_mappedFile_prefix(bamFile)
    out_bedgraphFile = base_path + '/' + bamFile_name.replace('.bam','') + '.bedgraph'
    
    if os.path.exists(out_bedgraphFile):
        print '[%s]: Expected output %s already present no processing required' % (SUB,out_bedgraphFile)
        return(out_bedgraphFile)

    #scale the read counts if asked
    if scale is True:
        readcount = get_mapped_read_count(bamFile_sorted)
        factor = 1 / (readcount / 1e6)
    else:
        factor = 1
     
    #bedgraph option line
    bedgraph_line = 'track type=bedGraph name=%s color=%s altColor=%s' % (bam_prefix,'43,131,186','171,121,164') 
    #find the genome coverage
    bam_bedtool = pybedtools.BedTool(bamFile_sorted)
    bam_bedtool.genome_coverage(bga=True, scale=factor,trackopts=bedgraph_line, output=out_bedgraphFile)
    
    print '[%s]: Converted bam %s to bedgraph file'  % (SUB,bamFile)
    return out_bedgraphFile



def bam_to_bigWig(bamFile,scale=False, sorted=False, **kwargs):
    """
    For a given bam file
    create a bigWig File
    """
    SUB = "bam_to_bigWig"
    #user supplied args if any
    force = kwargs.get('force',False)
    output_dir = kwargs.get('output_dir',False)

    #output bigWig filename
    out_bigWig_file = get_mappedFile_prefix(bamFile,output_dir=output_dir)  + '.bw'    
    if os.path.exists(out_bigWig_file):
        print '[%s]: Expected output %s already present no processing required' % (SUB,out_bigWig_file)
        return(out_bigWig_file)

    #constructing genome size file on the fly from bam file
    bam_fh = get_mappedFile_FH(bamFile)
    bam_prefix = get_mappedFile_prefix(bamFile)
    temp_genome_file = tempfile.NamedTemporaryFile(mode='w+t',suffix='.tmp', delete=False) #when delete is True file is deleted as soon as it closed
    [ temp_genome_file.write('%s\t%s\n' % (chr,len)) for chr,len in itertools.izip(bam_fh.references,bam_fh.lengths) ]
    temp_genome_file.close()


    #1. convert bam to bedgraph first 
    bedGraph_file = bam_to_bedGraph(bamFile, scale=scale, sorted=False, **kwargs )

    #2. convert to bigwig using the UCSC script   #should be in the path
    cmds = ['bedGraphToBigWig',bedGraph_file, str(temp_genome_file.name),out_bigWig_file]
    p = subprocess.Popen(cmds,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    
    if not stderr:
        print 'Created a bigWig file %s for %s' % (os.path.basename(out_bigWig_file),
                                                   os.path.basename(bamFile))
    else:
        print stderr
    return(out_bigWig_file)



def create_detailed_insert_size_plot(bamFile):
    SUB='create_insert_size_plot'
    
    namesorted_bamFile =  name_sort_bam(bamFile)
    base_path = os.path.dirname(namesorted_bamFile) or '.'
    bamFile_name = os.path.basename(namesorted_bamFile)
    out_insertSize_plotFile = base_path + '/' + bamFile_name.replace('.bam','') + '_insert_size.png'
    out_insertSize_logFile = base_path + '/' + bamFile_name.replace('.bam','') + '_insert_size.log'
    
    #get the dual logger
    logging = utils._get_logger(out_insertSize_logFile,logger_name = __name__)
    bamfile_handle = get_mappedFile_FH(namesorted_bamFile)
    bamfile_prefix = get_mappedFile_prefix(namesorted_bamFile)

    bag = collections.Counter()
    insert_size_bag = collections.defaultdict(list)
    
    for counter,read1,read2 in itertools.izip(itertools.count(1),bamfile_handle,bamfile_handle):
        pair_type = get_readPair_type(read1,read2)
        bag[pair_type] +=  1
        insert_size_bag[pair_type].append(abs(read1.isize))  #ref http://docs.python.org/2/library/collections.html, defaultdict examples
    
    #plotting the insert size
    alphas = {1:1,2:.7,3:.5,4:.3,5:.2,6:.1}
    for count,(pair_type,insert_sizes) in enumerate(insert_size_bag.items()):
        plt.hist([ math.log10(x+1) for x in insert_sizes] ,bins=100,label=pair_type,histtype='step')
        
        
        
    plt.xlabel('Insert Size(log 10) based on mapping')
    plt.ylabel('frequency')
    plt.legend()
    plt.savefig(out_insertSize_plotFile,dpi=175)
    
    overall_insert_size = [] 
    
#    print 'pair_type \t count \t mean \t std \t median \n'
    logging.info( '{0:#>70}'.format(''))
    logging.info('{0:30} {1:>10} {2:>10} {3:>10}'.format('pair_type','mean','std','median'))
    logging.info('{0:#>70}'.format('')) 
    for pair_type,insert_sizes in insert_size_bag.items():
        temp_array  = np.array(insert_sizes)
        count = bag[pair_type]
        overall_insert_size.extend(insert_sizes)
        logging.info('{0:30} {1:10.2f} {2:10.2f} {3:10.2f}'.format( pair_type, np.mean(temp_array),np.median(temp_array), np.std(temp_array)))
    
    #overall insert size
    temp_array = np.array(overall_insert_size)
    logging.info('{0:30} {1:10.2f} {2:10.2f} {3:10.2f}'.format( 'overall', np.mean(temp_array),np.median(temp_array), np.std(temp_array)))

    logging.info('[%s]: Created insert size plot for bam: %s'  % (SUB,bamFile_name))
    return out_insertSize_plotFile




def remove_low_coverage_reads(bamFile,min_coverage=4):
    SUB = 'remove_low_coverage_reads'
    
    bamFile_bai = bamFile + '.bai'
    if not os.path.exists(bamFile_bai):
        bamUtils.index_bam(bamFile)
    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    bamFile_handle = pysam.Samfile(bamFile,'rb')
    references = bamFile_handle.references
    out_bamFile = base_path + '/' + bamFile_name.replace('.bam','') + '_mincov_'+str(min_coverage)+'.bam'
    if os.path.exists(out_bamFile):
        print '[%s]: Expected output %s already present no processing required' % (SUB,out_bamFile)
        return out_bamFile
    out_bamFile_handle = pysam.Samfile(out_bamFile,'wb',template=bamFile_handle)
    count_num_reads_low_coverage = 0;
    for chromosome in references:
        low_coverage_read_headers = {}
        for pileupcol in bamFile_handle.pileup(chromosome):
           if pileupcol.n < min_coverage and pileupcol.n > 0:
                for read in pileupcol.pileups:
                    low_coverage_read_headers[read.alignment.qname] = 1
        for read in bamFile_handle.fetch(chromosome):
            if ( low_coverage_read_headers.get(read.qname,'None') == 'None'):
                out_bamFile_handle.write(read)
            else:
                count_num_reads_low_coverage += 1
    
    print '[%s]: Wrote file %s, found and skipped %d low coverage reads' % (SUB,os.path.basename(out_bamFile), count_num_reads_low_coverage)
    return out_bamFile


def generate_read_count_per_chr_per_bam(bam_files):
    '''
    given a list of bamfiles generate a pandas dataframe of 
    read count per reference in the bamfile and the strand
    '''
    SUB = 'generate_read_count_per_chr_per_bam'
    
    class AutoVivification(dict):
        """Implementation of perl's autovivification feature."""
        def __getitem__(self, item):
            try:
                return dict.__getitem__(self, item)
            except KeyError:
                value = self[item] = type(self)()
                return value
        
    
    for bam_file in bam_files:
        base_dir = os.path.dirname(bam_file)
        file_name = os.path.basename(bam_file)
        out_file = base_dir + '/' + file_name +'.refCounts.txt'
        
        
        if os.path.exists(out_file):
            print '[%s]: Expected output %s already present no processing required' % (SUB,file_name)
            continue
        else:
            #else open the file
            out_fh = open(out_file,'w')
            print '[%s]: Processing %s' % (SUB,file_name)
            bam_file_fh = pysam.Samfile(bam_file,'rb')
            hash = AutoVivification()
            for read in bam_file_fh:
                if read.is_unmapped:
                    continue
                else:
                    reference = read.tid
                    strand    = get_strand_for_firstStrand_RNA_Seq_read(read)
                if hash[reference][strand]:
                     hash[reference][strand] += 1
                else:
                    hash[reference][strand] = 1
            references = bam_file_fh.references
            lengths =   bam_file_fh.lengths
        
            ref_to_len_dict = {}
            #creating a reference to length dict
            for (ref,len) in itertools.izip(references,lengths):
                ref_to_len_dict[ref.strip()] = len
            
            #writing to disk
            out_fh.write('%s \t %s \t %s \t %s \n' %('reference','reference_length','count_positive_strand','count_negative_strand'))
            for ref_id in hash:
                ref_name    = bam_file_fh.getrname(ref_id).strip()
                count_pos   = hash[ref_id]['+'] or 0
                count_neg   = hash[ref_id]['-'] or 0
                if ref_to_len_dict.get(ref_name, None) is None:
                    ref_len='NA'
                    print '[%s]: No reference length found for %s' % (SUB,ref_name)
                else:
                    ref_len = ref_to_len_dict[ref_name]
                out_fh.write('%s \t %s \t %s \t %s\n' %(ref_name,ref_len,count_pos,count_neg))
            out_fh.close()
            bam_file_fh.close()

  
def split_bam_by_chromosomes(bamFile,dir_name=None):
    
    SUB ='split_bam_by_chromosomes'
    #print 'pysam version available %s' % pysam.__version__
    
    bamFile_bai = bamFile + '.bai'
    if not os.path.exists(bamFile_bai):
        bamUtils.index_bam(bamFile)
    
    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    file_prefix=os.path.basename(bamFile).replace('.bam','')
    
    ##check and create a folder to keep the split bam files
    if not dir_name:
        dir_name = 'bams_split_by_chr'
    base_path_split_files = base_path + '/' + dir_name
    if not os.path.exists(base_path_split_files):
        os.makedirs(base_path_split_files)
    
    #open the bam file
    in_bamFile_handle = pysam.Samfile(bamFile,'rb')
    references = in_bamFile_handle.references
    
    # initialize the container to hold the bam file names splitted by chromosome
    list_of_split_bam_files = []
    list_of_chromosomes = []
    for chr in references:
        
        #forming the bam file name for this chr
        chr_bam = base_path_split_files + '/' + file_prefix + '_' + chr + '.bam'
        
        #check if chr_bam exists already, if so just store the name and continue
        if os.path.exists(chr_bam):
            print '[%s]: for %s bam %s exists' % (SUB, chr, chr_bam)
            list_of_split_bam_files.append(chr_bam)
            list_of_chromosomes.append(chr) 
        else:
            chr_bam_handle = pysam.Samfile(chr_bam,'wb', template=in_bamFile_handle)
            
            for read in in_bamFile_handle.fetch(chr):
                if read:
                    chr_bam_handle.write(read)
                    
            print '[%s]: created bam %s ' % (SUB,chr_bam)
            list_of_split_bam_files.append(chr_bam)
            list_of_chromosomes.append(chr)
        
    return (list_of_split_bam_files,list_of_chromosomes)


##
## Taken from https://github.com/daler/pybedtools/blob/master/pybedtools/contrib/bigwig.py
def get_mapped_read_count(bam, force=False):
    """
    Scale is cached in a bam.scale file containing the number of mapped reads.
    Use force=True to override caching.
    """
    scale_fn = bam + '.scale'
    if os.path.exists(scale_fn) and not force:
        for line in open(scale_fn):
            if line.startswith('#'):
                continue
            readcount = float(line.strip())
            return readcount

    cmds = ['samtools',
            'view',
            '-c',
            '-F', '0x4',
            bam]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise ValueError('samtools says: %s' % stderr)

    readcount = float(stdout)

    # write to file so the next time you need the lib size you can access
    # it quickly
    if not os.path.exists(scale_fn):
        fout = open(scale_fn, 'w')
        fout.write(str(readcount) + '\n')
        fout.close()
    return readcount


def plot_bamCoverage(bamFile,min_coverage=2):
    SUB ='plot_bamCoverage'
    
    bamFile_bai = bamFile + '.bai'
    if not os.path.exists(bamFile_bai):
        index_bam(bamFile)
    
    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    bamFile_handle = pysam.Samfile(bamFile,'rb')
    references = bamFile_handle.references
    out_genome_cov_plot = base_path + '/' + bamFile_name.replace('.bam','') + '_mincov_'+str(min_coverage)+'.png'
    out_cov_hist = base_path + '/' + bamFile_name.replace('.bam','') + '_histogram_mincov_'+str(min_coverage)+'.png'
    
    if os.path.exists(out_genome_cov_plot):
        print 'Expected output %s already present ...! No processed required' % (out_genome_cov_plot)
        return out_genome_cov_plot
    
    nuc_coverage = []
    nuc_coverage_gt_min_coverage = []
    for chromosome in references:
        for pileupcol in bamFile_handle.pileup(chromosome,max_depth=50000):
           if pileupcol.n > min_coverage:
                nuc_coverage.append(pileupcol.n)
                nuc_coverage_gt_min_coverage.append(pileupcol.n)
           else:
                nuc_coverage.append(0)
    
    #plot1
    plt.figure()
    plt.plot(nuc_coverage,'-')
    plt.xlabel('Genome')
    plt.ylabel('coverage')
    plt.savefig(out_genome_cov_plot)
    
    #plot2
    plt.figure()
    plt.hist([ math.log10(x+1) for x in nuc_coverage_gt_min_coverage],bins=50)
    plt.grid(True)
    plt.xlabel('Coverage(log 10)')
    plt.ylabel('frequency')
    plt.savefig(out_cov_hist,dpi=175)
    
    temp_array = np.array(nuc_coverage_gt_min_coverage)
    print '{0:40} {1:10.2f} {2:10.2f} {3:10.2f} '.format( 'Coverage(mean, std, median)',
                                                          np.mean(temp_array),
                                                          np.std(temp_array),
                                                          np.median(temp_array))
        
    print '[%s]: created Coverage plots: \n 1. %s \n 2. %s ' % (SUB,out_genome_cov_plot,out_cov_hist)
    



def bam_to_fastq(bamFile, paired=False,**kwargs):

    bamPrefix = get_mappedFile_prefix(bamFile)

    if paired is True:
        sorted_bam = name_sort_bam(bamFile)
        read1 = bamPrefix + '_read1.fastq'
        read2 = bamPrefix + '_read2.fastq'
        pybedtools.BedTool.bam_to_fastq(sorted_bam,fq=read1,fq2=read2,**kwargs)
        return(read1,read2)
    else:
        fastq = bamPrefix + '.fastq'
        pybedtools.BedTool.bam_to_fastq(bamFile,fq=fastq,**kwargs)
        return(fastq)


def bam2Fastq(bamFile,paired=False):
    """
    Bam to Fastq using bedtools:
    bedtools bamtofastq is a conversion utility for extracting FASTQ records from sequence alignments in BAM format. 
    Can also create two FASTQ files for paired-end sequences, but must be sorted first.
    """
    bamPrefix = get_mappedFile_prefix(bamFile)

    if paired is True:
        read1 = bamPrefix + '.read1.fastq'
        read2 = bamPrefix + '.read2.fastq'
        sorted_bam = name_sort_bam(bamFile)

        p=subprocess.Popen('bedtools bamtofastq -i %s -fq %s -fq2 %s' % (sorted_bam,read1,read2),stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        if stderr:
            raise ValueError('bedtools bamtofastq -i says: %s' % stderr)
        return(read1, read2)
    else: 
        fastq = bamPrefix + '.fastq'
        p=subprocess.Popen("bedtools bamtofastq -i %s -fq %s" % (bamfile, fastq),stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        if stderr:
            raise ValueError('bedtools bamtofastq says: %s' % stderr)
        return(fastq)


def get_strand_for_RNA_Seq_read(read,lib_protocol):
    '''
    given that the input RNA-Seq data is based on  first strand sequencing lab protocol
    for a given pair, 
    a.) read pair is coming from sense strand if
        read 1 is reverse and read 2 is not reversed
        
    b.) read pair is coming from antisense strand if 
        read 1 is not reversed and read 2 is reversed
    

    ##############
    Background theory for the first strand protocol
    ##############
    
    Stranded RNA-Seq protocol is stranded with respect to mRNA and not with repsect to genome
    when we are trying to split the reads to sense and antisense strand we are basically going to split them wrt to mRNA
    
        in first strand RNA-Seq the first strand is actually the rev comp of RNA is being sequenced
        so what is being sequenced is actually the reverse strand of the mRNA,
        
        : so read 1 will map to the opposite strand of RNA and read 2 on the same
        
    
    Specific Example:
        --------------------------------> (genome 5' to 3')
            (r2)--->        <---- (r1)
            
        since read 1 is coming from the first strand of mRNA and in this case r1 maps on the opp strand of the reference,
        this actually is reflecting the read coming from a sense mRNA and thus shud be a sense read pair, going to sense strand coverage bin

        
    !!!!!!!!!!
    ***just the reverse for second strand protocol****
    !!!!!!!!!
    
    Assumption:
        Input data is JGI first strand normally
        for second strand libs : we have to swap the strand
        paired end, facing inwards
    '''

    
    if read.is_unmapped:
        return 'unmapped'
    
    strand = ''
    
    #paried reads
    if read.is_paired:
        if ( ( read.is_read1 and not read.is_reverse ) or (read.is_read2 and read.is_reverse) ):
            strand = '-'
        elif ( ( read.is_read1 and read.is_reverse ) or ( read.is_read2 and not read.is_reverse ) ):
            strand = '+'
    
    #single reads
    #WARNING: if the paired end reads are mapped as
    #         single end reads then wierd results will be produced
    #        for strandedness check
    else:
        if (read.is_reverse):
            strand = '+'
        elif (not read.is_reverse):
            strand = '-'
    
    
    #logic to swap the strand if the lib protocol is second stranded
    if lib_protocol == 'first-strand':
        return strand
    elif lib_protocol == 'second-strand':
        if strand == '+':
            return '-'
        elif strand == '-':
            return '+'
    else:
        print 'Cant understand the library protocol: %s' % lib_protocol
        sys.exit(2)
    
    '''
     #logic for read 1
    if read.is_read1:
        if read.is_reverse:
        # read is read1 and reversed : so coming from same strand of first Strand cDNA which is revComp of mRNA, so positive strand 
        else:
        # read is read1 and NOT reversed : so coming from opp strand of first Strand cDNA which is same as mRNA, so negative strand

    #opp logic for read 2
    if read.is_read2:
        if read.is_reverse:
        # read is read2 and reversed : so coming from opp strand of first Strand cDNA which is same as mRNA, so negative strand 
        else:
        # read is read2 and NOT reversed : so coming from same strand of first Strand cDNA which is revComp of mRNA, so positive strand
    '''
    #END get_strand_for_firstStrand_RNA_Seq_read
    
        
    
def get_mappedFile_FH(mapped_file,mode='r'):
    """
    get apt file handle for reading sam / bam file
    """
    SUB='get_mappedFile_FH'
    
    if mapped_file.endswith('.bam'):
        file_handle = pysam.Samfile(mapped_file,('%sb' % mode))
        return (file_handle)
    elif mapped_file.endswith('.sam'):
        file_handle = pysam.Samfile(mapped_file,('%s' % mode))
        return (file_handle)
    else:
        print '[%s]: cant understand the file extension..tried looking for .sam  and .bam' % SUB

#bam summary




def get_bamStats(bamFile,**kwargs):
    '''
    generate bamStats from a sam/bam file and save the same 
    secondary mapped reads are counted separately
    '''
    #stats that will be recorded
    bamStats = { 'bamFile': bamFile,
                 'reads' : 0,
                 'read1' : 0,
                 'read2' : 0,
                 'read2' : 0,
                 'mapped' : 0,
                 'mapped_pairs' : 0,
                 'chimeric_mapped' : 0,
                 'secondary_mapped' : 0,
                 'unmapped' : 0,
                 'mapped_read1' : 0,
                 'mapped_read2' : 0,
                 'secondary_mapped_read1' : 0,
                 'secondary_mapped_read2' : 0,
                 'unmapped_read1': 0,
                 'unmapped_read2': 0,
                 'QC_failed' : 0,
                 'PCR_optical_dups' : 0
                 }
    
    
    #user supplied args if any
    force = kwargs.get('force',False)
    output_dir = kwargs.get('output_dir',False)
    
    #stats to be saved with the following file name
    bamStats_file = get_mappedFile_prefix(bamFile,output_dir=output_dir) + '.bamStats'
    
    #check is pre-computed stats exist
    if os.path.exists(bamStats_file) and force == False:
        #read the stats from the file and return
        bamStats = yaml.load(open(bamStats_file,'r'))
        return bamStats
    else:
        #do the read counting and save the stats
        #open the bam file for reading
        bam_fh = get_mappedFile_FH(bamFile)
        #start the timer
        start = time.clock()
        loop_start = time.clock()
        for count,read in enumerate(bam_fh):
            #counter
            if count % 1000000 == 0 and count > 0:
                time_taken = time.clock() - start
                total_time_taken = time.clock() - loop_start
                sys.stdout.write('Processed 1,000,000 reads in %0.2f seconds total_reads: %d in %0.2f seconds \n' % (time_taken,count,total_time_taken))
                start = time.clock()
            
            #check for multi mapped reads first up and keep a separate counter for them
            if read.is_secondary:
                bamStats['secondary_mapped'] += 1
                if read.is_read1:
                    bamStats['secondary_mapped_read1'] += 1
                elif read.is_read2:
                    bamStats['secondary_mapped_read2'] += 1
                # skip counting this read further
                continue
            #get the raw read count
            else:
                # count raw reads
                bamStats['reads'] += 1
                if read.is_read1:
                    bamStats['read1'] += 1
                if read.is_read2:
                    bamStats['read2'] += 1
    
            #check QC failed reads
            if read.is_qcfail:
                bamStats['QC_failed'] += 1
                continue
            
            #PCR / Optical dups
            if read.is_duplicate:
                bamStats['PCR_optical_dups'] += 1
                
            #check if unmapped reads
            if read.is_unmapped:
                bamStats['unmapped'] += 1
                if read.is_read1:
                    bamStats['unmapped_read1'] += 1
                elif read.is_read2:
                    bamStats['unmapped_read2'] += 1
            #if read is mapped
            else:
               bamStats['mapped'] += 1
               #chimeric read pairs
               if read.tid != read.rnext and read.rnext != '*':
                   bamStats['chimeric_mapped'] += 1
               if read.is_read1:
                   bamStats['mapped_read1'] += 1
               elif read.is_read2:
                   bamStats['mapped_read2'] += 1

            #check if read is part of a properly mapped pair
            if read.is_proper_pair:
                bamStats['mapped_pairs'] += 1
                    
        #save the stats in a YAML format
        with open(bamStats_file,'w') as fh:
            fh.write(yaml.dump(bamStats, default_flow_style=False))
        print 'Wrote %s file' % bamStats_file
        return (bamStats)
              
'''
#check if the the pair is mapped on the same chromosome of other for one of the reads        
                    if read.is_read1:
                        if read.is_proper_pair:
                            bamStats['mapped_pairs'] += 1
                        else:
                            bamStats['mapped_pairs'] += 1
                            bamStats['chimeric_mapped_pairs'] += 1
            
    
                    if read.is_read1 and not read.is_unmapped:
                        if read.tid == read.rnext:
                            bamStats['mapped_pairs'] += 1
                        else:
                            bamStats['mapped_pairs'] += 1
                            bamStats['chimeric_mapped_pairs'] += 1
'''

def get_mappedFile_prefix(mapped_file,mode='r',output_dir=None):
    """
    get prefix for a sam/bam file
    output_dir : is used to append that dir as a dir name behind the file prefix
               : useful to direct the bam/sam file to a diff direc
    """
    SUB='get_mappedFile_prefix'
    
    if output_dir:
        dirname = output_dir
    else:
        dirname = os.path.dirname(mapped_file)
    if dirname == "": ##fix when script is run from the same dir where the file is
        dirname = '.'
    if mapped_file.endswith('.bam'):
        file_prefix = dirname + '/' +  os.path.basename(mapped_file).replace('.bam','')
        return (file_prefix)
    elif mapped_file.endswith('.sam'):
        file_prefix = dirname + '/' + os.path.basename(mapped_file).replace('.sam','')
        return (file_prefix)
    else:
        print '[%s]: cant understand the file extension..tried looking for .sam  and .bam' % SUB
    
    
        
def split_RNASeq_reads_bam_by_strand(mappedFile, library_protocol=None, create_split_files=True,force_run=False, **kwargs):
    '''
    given an input bam file RNA-Seq reads (with firststrand/ secondstrand protocol) 
    split it into 
       1. reads capturing cDNA ends on the sense strand
       2. reads capturing cDNA ends on the antisense strand
       
       PS: must specify which protocol of library construction is used
       
       library_protocol can take the following values:
           a.) first-strand
           b.) second-strand
    '''
    
    SUB = 'split_RNASeq_reads_bam_by_strand'
    
    base_dir = os.path.dirname(mappedFile) or '.'
    input_mappedFile_handle  = get_mappedFile_FH(mappedFile)
    input_mappedFile_prefix  = get_mappedFile_prefix(mappedFile)
    log_file_name =  base_dir+'/'+input_mappedFile_prefix+'_strand_split.log'
    
    
    if kwargs.get('logger',None) is None:
        #get a new dual logger
        logger = utils._get_logger(log_file_name,logger_name = __name__)
    else:
        #used the passed one
        logger = kwargs['logger']
    
    
    if library_protocol is None:
        logger.error('Error[%s]: library protocol is not mentioned \n\n exting with error code: 2')
        sys.exit(2)
    
    
    sense_bamFile     = base_dir+'/'+input_mappedFile_prefix+'_sense.bam'
    antisense_bamFile = base_dir+'/'+input_mappedFile_prefix+'_antisense.bam';
    if create_split_files:    
        if ( os.path.exists(sense_bamFile) and os.path.exists(antisense_bamFile) and ( force_run == False) ):
            logger.info( '[%s]: Expected output  \n %s \n %s  already present .....no processing required' % (SUB,os.path.basename(sense_bamFile), os.path.basename(antisense_bamFile)))
            return (sense_bamFile,antisense_bamFile)
        else:
            sense_bamFile_handle = pysam.Samfile(sense_bamFile,'wb',template=input_mappedFile_handle)
            antisense_bamFile_handle = pysam.Samfile(antisense_bamFile,'wb',template=input_mappedFile_handle)
    

    logger.info('splitting %s based on lib protocol: %s' % (mappedFile,library_protocol))
    
    
    #initialize counters
    count_num_sense_strand_reads = 0;
    count_num_antisense_strand_reads = 0;
    count_num_unmapped_reads = 0;
    count_num_mapped_reads = 0
     
    prog_start = time.clock()
    start_loop = prog_start
    for (count_read,read) in itertools.izip(itertools.count(1),input_mappedFile_handle):
        if (count_read % 1000000 == 0):
            end_loop = time.clock()
            time_taken = end_loop - start_loop
            logger.info('[%s]: Processed %d reads in %d secs' % (SUB,count_read,time_taken))
    
        strand = get_strand_for_RNA_Seq_read(read,library_protocol)
        if strand == '+':
            count_num_sense_strand_reads += 1
            if create_split_files:
                sense_bamFile_handle.write(read)
            count_num_mapped_reads += 1
        elif strand == '-':
            count_num_antisense_strand_reads += 1
            if create_split_files:
                antisense_bamFile_handle.write(read)
            count_num_mapped_reads += 1
        elif strand == 'unmapped':
            count_num_unmapped_reads += 1
            continue
        else:
            logger.info("[%s]: Warn-> can't determine strand for read %s" % (SUB,read.qname))
            continue
        
    #closefile handles
    input_mappedFile_handle.close()
    if create_split_files:
        antisense_bamFile_handle.close()
        sense_bamFile_handle.close()
        #index split bams
        index_bam(sense_bamFile)
        index_bam(antisense_bamFile)
    
    
    prog_end = time.clock()
    prog_time_taken = prog_end - prog_start
    
    logger.info('#########################')
    logger.info('Summary')
    logger.info('#########################')
    logger.info('#reads input bam : %d' % (count_read))
    logger.info('#reads sense bam : %d' % (count_num_sense_strand_reads))
    logger.info('#reads antisense bam : %d' % (count_num_antisense_strand_reads))
    logger.info('#mapped reads : %d' % (count_num_mapped_reads))
    if create_split_files:
        logger.info('Created sense bam : %s' % (os.path.basename(sense_bamFile)))
        logger.info('Create antisense bam : %s' % (os.path.basename(antisense_bamFile)))
    logger.info('#total time : %d'    % (prog_time_taken))
            
    return (sense_bamFile,antisense_bamFile,
            count_read,count_num_mapped_reads,
            count_num_unmapped_reads,
            count_num_sense_strand_reads,count_num_antisense_strand_reads)
    
    
#END split_first_Strand_reads_bam_by_strand



def debug(read1,read2,message):
    """
    print the debug information
    for a read pair
    """
    SUB='debug'
    print '[%s]: Message: %s' % (SUB,message)
    print '%s \n %s \n type: %s \n' % ( read1, read2, type(read1))


    
def is_chimeric_pair(read1,read2):
    ##both reads should map to diff chromosome
    if read1.tid != read2.tid:
        return True
    else:
        return False

def is_a_sequencing_readPair(read1,read2):
    
    if read1.qname != read2.qname:
        return False
    else:
        return True


def is_intraChrom_pair(read1,read2):
    ##both reads should map to same chromosome
    if read1.tid == read2.tid:
        return True
    else:
        return False

def get_insert_size_readPair(read1,read2):
    '''
    for the pairs where read1 and read2 are mapped independently
    '''
    
    if is_chimeric_pair(read1,read2):
        return '*'

    #finding the absolute insert size of the fragment
    abs_insert_size = 0
    if (read1.pos < read2.pos):
        abs_insert_size = read2.pos + read2.qlen - read1.pos
        return abs_insert_size
    elif (read2.pos < read1.pos):
        abs_insert_size = read1.pos + read1.qlen - read2.pos
        return abs_insert_size
    elif (read1.pos == read2.pos):
        abs_insert_size = 0
        return abs_insert_size
    else:
        debug_messsage= "Can't determine the abs insert size"
        debug(read1,read2,debug_messsage)
        return 0


def get_readPair_type(read1,read2):
    """
    given a read pair, read1 and read2
    the function checks exhaustively for all possibilities 
    
    1. chimera
    2. singleton
    3. outwards
    4. inwards
    5. both-forward
    6. both-reverse
    
    and return the applicable one

    """
    
    
    if not is_a_sequencing_readPair(read1,read2):
        return 'not-a-sequencing-pair'
    
    if read1.is_duplicate or read2.is_duplicate:
        return 'duplicate'
    
    if is_chimeric_pair(read1,read2):
        return 'chimera'
    
    if read1.is_unmapped or read2.is_unmapped:
        return 'singleton'
    
    #case 1 : read 1 maps on the opp strand 
        #two sub cases possible
    if (read1.is_reverse and not read2.is_reverse):
        # <------(r1) ----->(r2)
        if (read1.pos <= read2.pos):
            return 'outwards'
        # ----->(r2)  <------(r1) 
        elif (read2.pos <= read1.pos):
            return 'inwards'
        
        
    #case 2 : read 1 is on the same strand
        #two sub cases possible
    elif (not read1.is_reverse and read2.is_reverse):
        # <---------(r2)  -------->(r1)
        if (read2.pos <= read1.pos):
            return 'outwards'
        # ------>(r1)    <-------(r2)
        elif (read1.pos <= read2.pos):
            return 'inwards'
        
        
    #case 3 : read 1 and 2 on the same strand
        #two sub cases possible
    elif (not read1.is_reverse and not read2.is_reverse):         
        # ------>(r2) ------>(r1)
        if (read2.pos <= read1.pos):
            return 'both-forward'
        # ------>(r1) ------>(r2)
        elif (read1.pos <= read2.pos):
            return 'both-forward'
        
    
    #case 4 : read 1 and 2 map on the opp strand
        #two sub cases possible
    elif ( read1.is_reverse and read2.is_reverse):
        # <------(r1) <------(r2)
        if (read2.pos <= read1.pos):
            return 'both-reverse'
        # <------(r2) <------(r1))
        elif (read1.pos <= read2.pos):
            return 'both-reverse'
    
    
    #Error Case
    else:
        message="Can't determine read pair type in SUB: %s" % (SUB)
        debug(read1,read2,message)
        return (False)
    
    return True





def feature_cov_in_bam(features_bed_file,bams,stranded=None,
                       min_map_qual=1,lib_protocol=None,
                       sample_names=None,count_dup=True,
                       normalize=True):
    
    """
        For a given feature list in BED format
        the method produces a list of read count/bam(sample) across all features
        
        //TODO
        implement normalization process
        
    
    """

    curr_dir = os.getcwd()
    log_file = curr_dir + '/features_cov.log'
    feature_cov_file = curr_dir + '/features_cov.tsv'
    
    logger = utils._get_logger(log_file,logger_name = __name__ )
    
    temp_list_feature_cov = []

    #standed counting for a RNA-Seq stranded lib    
    if (stranded):
        if lib_protocol is None:
            logger.info("No library type defined ..required for stranded counting")
        

        #split the bam into sense and antisense
        split_bam_info = [ split_RNASeq_reads_bam_by_strand(bam,library_protocol=lib_protocol,logger=logger) for bam in bams ]
        
        #get the list of sense and antisense bam
        sense_bams = [ info[0] for info in split_bam_info ]
        antisense_bams = [ info[1] for info in split_bam_info ]
        
        #split the features into sense and antisense for proper counting
        (sense_features,antisense_features)=bio_format_parsers.split_bed_by_strand(features_bed_file)
        
        sense_features_bedtool = pybedtools.BedTool(sense_features)
        antisense_features_bedtool = pybedtools.BedTool(antisense_features)
        
        #calculate the cov (#reads/feature) for multiple bam files
        sense_features_cov = sense_features_bedtool.multi_bam_coverage(bams=sense_bams,q=min_map_qual,D=count_dup)
        antisense_features_cov = antisense_features_bedtool.multi_bam_coverage(bams=antisense_bams,q=min_map_qual,D=count_dup)
        
        #chain the sense and antisense features
        feature_cov = pybedtools.BedTool(itertools.chain(sense_features_cov,antisense_features_cov))
        
        #create a temp list for conversion into pandas df
        temp_list_feature_cov = [x[:] for x in feature_cov]
        
    #non stranded counting
    else:
        features_bedtool = pybedtools.BedTool(features_bed_file)
        feature_cov = features_bedtool.multi_bam_coverage(bams=bams,q=min_map_qual,D=count_dup).saveas(feature_cov)
        
        #create a temp list for conversion into pandas df
        temp_list_feature_cov = [x[:] for x in feature_cov]
        
    
    #columns naming
    column_names = ['chr','start','end','feature_name','score','strand']
    if sample_names is not None:
        column_names.extend(sample_names)    
    else:
        #extract the first 5 characters from the bam file
        #assuming the lib name is prefixed to the bam file
        [ column_names.append(os.path.basename(bam)[:5].replace('.','')) for bam in bams ]
        
    #convert the data into pandasdf
    df = pandas.DataFrame.from_records(temp_list_feature_cov,columns=column_names,coerce_float=True)

    df.to_csv(feature_cov_file,sep="\t",na_rep='NA',index=False)
    
    print "the read count / feature / bam is saved in %s file" % (feature_cov_file)
    
    if normalize:
        print 'Normalization of counts not ready in the method..TODO'
    
    #pickle.dump(df,open('df.pickle','wb'))
    



class mappingSummary(object):
    """
    given a bam file generate the alignment quality stats
    
    #########
    example:
    ##########
    
    #instantiate
    bamSummary = mappingSummary(bamFile,max_readlen=150,quality_offset=33,minQual=30)
    
    #evaluate the indels, mismatches, N%
    # returns a pandas df
    df = bamSummary.eval_alignQ()


    """
    
    #initializer
    def __init__(self,bam,max_readlen,quality_offset=33,minQual=20):
        self.bam = bam
        self.quality_offset = quality_offset
        self.max_readlen = max_readlen
        self.minQual = minQual
        
        #stats with all the reads
        self.mismatches = np.zeros(max_readlen)
        self.indels = np.zeros(max_readlen)
        self.qual = np.zeros(max_readlen)
        
        #stats with reads where no N is found
        self.mismatches_rN = np.zeros(max_readlen)
        self.indels_rN = np.zeros(max_readlen)
        self.qual_rN = np.zeros(max_readlen)
        
        self.Ns     =  np.zeros(max_readlen)    
        self.num_reads = 0
        self.num_reads_rN = 0
        self.unmapped_reads = 0
        self.Ns_unmapped_reads = np.zeros(max_readlen)
        
        #compile the patternt for matching the MD
        self.MD_pattern = '([0-9]+)([A-Z]|\^[A-Z]+)*'
        self.MD_pattern_compiled = re.compile(self.MD_pattern)
        
    
    
    def eval_alignQ(self):
        """
            #main looper
            #which iterates over bam file
        """
        bam_fh=pysam.Samfile(self.bam,'rb')
        
        for count,read in enumerate(bam_fh):
            self.num_reads += 1
            read_N_seen = self._check_Ns_across_read(read)
            if read_N_seen is False:
                self.num_reads_rN += 1
            self._check_quality_across_read(read,read_N_seen)
            self._check_mismatch_indels_across_read(read,read_N_seen)
            
            if self.num_reads % 1000000 == 0 and self.num_reads > 0:
                print '%s : processed %s reads' % (os.path.basename(self.bam),self.num_reads) 
        print '%s : processed %d reads' % (os.path.basename(self.bam),self.num_reads) 
        
        #create a summary df to be returned
        df = pandas.DataFrame( { 'N_percent'                     :     (self.Ns/self.num_reads)*100,
                                 'N_percent_unmapped_reads'      :     (self.Ns_unmapped_reads/self.unmapped_reads)*100, 
                                 'mismatches'     :     (self.mismatches/self.num_reads)*100,
                                 'indels'         :     (self.indels/self.num_reads)*100,
                                 'bases_<_Q30'    :     (self.qual/self.num_reads)*100,
                                 'mismatches_rN'  :     (self.mismatches_rN/self.num_reads_rN)*100,
                                 'indels_rN'      :     (self.indels_rN/self.num_reads_rN)*100,
                                 'bases_<_Q30_rN' :     (self.qual_rN/self.num_reads_rN)*100
                               }
                              )
        #store the db for later access
        self.results = df
        return df
        
    
    
    
    def plot_alignment_quality(self):
        """
            plot the results
        """
        
    
        df = self.results
        fig = plt.figure(figsize=(15,10))
        fig.subplots_adjust(hspace=.2,wspace=.3)
        fig.suptitle('Alignment Quality Summary:  %s ' % os.path.basename(self.bam) ,size=14)
        
        ax1 = fig.add_subplot(221)
        df[df.index <= self.max_readlen - 30][['N_percent','mismatches','mismatches_rN']].plot(ax=ax1)
        ax1.set_xlabel('read position')
        ax1.set_ylabel('Percent')
        ax1.set_title('Mismatches v/s Ns percent per read bp position')
        
        ax2 = fig.add_subplot(222)
        df[df.index <= self.max_readlen - 30][['indels','indels_rN']].plot(ax=ax2)
        ax2.set_xlabel('read position')
        ax2.set_ylabel('Percent')
        ax2.set_title('Indels v/s Ns percent per read bp position')
            
        ax3 = fig.add_subplot(223)
        df[df.index <= self.max_readlen - 30][['bases_<_Q30','bases_<_Q30_rN']].plot(ax=ax3)
        ax3.set_xlabel('read position')
        ax3.set_ylabel('Percentage  Bases < Q30')
        ax3.set_title('percent bases < Q30 per cycle')
        
        ax4 = fig.add_subplot(224)
        df[df.index <= self.max_readlen - 30][['N_percent_unmapped_reads']].plot(ax=ax4)
        ax4.set_xlabel('read position')
        ax4.set_ylabel('Percent')
        ax4.set_title('UNMAPPED reads N percent per read bp position')
        
        outplot_name = self.bam + '_alignQ.png'
        plt.savefig(outplot_name)
        
        
    
    
    def _check_mismatch_indels_across_read(self,read,read_N_seen):
        if read.is_unmapped:
            return None
        
        #get the cigar and split MD tag
        cigar,split_MD_tag = self._parse_MD_tag(read)
           
        #jumper for MD tag
        CIAGR_jumper_for_MD = 0
    
        #make sense of CIGAR string
        position_cigar = 0
        for count,(cigar_code,matches) in enumerate(cigar):
            """    
            CIGAR CODE in pysam parsed info
            M    BAM_CMATCH        0
            I    BAM_CINS          1
            D    BAM_CDEL          2
            N    BAM_CREF_SKIP     3
            S    BAM_CSOFT_CLIP    4
            H    BAM_CHARD_CLIP    5
            P    BAM_CPAD          6
            =    BAM_CEQUAL        7
            X    BAM_CDIFF         8
            """
        
            #matches or mismatches
            if cigar_code == 0:
                position_cigar += matches
            
            #hard or soft clipping
            elif cigar_code == 4 or cigar_code == 5:
                self.mismatches[position_cigar:position_cigar+matches] += 1
                if read_N_seen is False:
                    self.mismatches_rN[position_cigar:position_cigar+matches] += 1
                position_cigar += matches
                if count == 0: #if the first item in the CIGAR is soft or hard clipping
                    #calculating the jump ahead for MD tag
                    #find any soft or hard clippings in the beginning of CIGAR
                    #unfortunately they are not reported in MD and better be calculated here
                    CIAGR_jumper_for_MD = matches
            
            #insertions in the read wrt to the reference genome
            elif cigar_code == 1:
                self.indels[position_cigar] += 1
                if read_N_seen is False:
                    self.indels_rN[position_cigar] += 1
                position_cigar += matches
            
            #deletion in the read wrt to the reference genome :
            elif cigar_code == 2:
                self.indels[position_cigar] += 1
                if read_N_seen is False:
                    self.indels_rN[position_cigar] += 1
                #cigar position doesnt change in a deletion as the read length stays same
                #the deletion is in the read wrt to the reference
                #position_cigar += matches
                
                #Example
                #   AGTGATGGG---------GGGGTTCCAGGTGGAGACGAGGACTCC   ( seqeunced read)
                #     --     ---------
                #   AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC  (reference genome)
            else:
                print 'New CIGAR code seen: %d' % cigar_code
                print 'PROBS :%s, MD:%s , CIGAR:%s, position:%d' % (read.qname,MD_tag,read.cigar,position_cigar)
                sys.exit(1)
        
        
        #make sense of MD tag
        #mainly for mismatches
        position_MD = 0
        #jump ahead if there were soft/hard clipping in CIGAR string
        position_MD += CIAGR_jumper_for_MD
        
        for x in split_MD_tag:
            if x == '' or x is None:
                continue
            #mismatches
            if x in ['A','G','C','T','N'] and not x.startswith('^'):  #^ indicate deletion in MD string
                self.mismatches[position_MD:position_MD+len(x)] += 1
                if read_N_seen is False:
                    self.mismatches_rN[position_MD:position_MD+len(x)] += 1
                position_MD += len(x) #now increase the position of MD
            
            #see if string is a integer AND move ahead the position in string
            try:
                x=int(x)
                position_MD += x
            except ValueError:
                pass
    
    
    def _check_Ns_across_read(self,read):
        if read.is_reverse:
            seq = read.seq[::-1]
        else:
            seq = read.seq
        read_positions_with_N =[matches.start() for matches in re.finditer('N',seq)]
        if read_positions_with_N:
            self.Ns[read_positions_with_N] += 1
            #unmapped reads stats
            if read.is_unmapped:
                self.unmapped_reads += 1
                self.Ns_unmapped_reads[read_positions_with_N] += 1
            return True
        else:
            return False
        
    def _check_quality_across_read(self,read,read_N_seen):
        #get the indices of the read position with quality score < min_quality
        if read.is_reverse:
            qual = read.qual[::-1]
        else:
            qual = read.qual
        indices= [ count for count,q in enumerate(qual) if (ord(q)- self.quality_offset) < self.minQual ]
        self.qual[indices] += 1
        
        if read_N_seen is False:
            self.qual_rN[indices] += 1
    
    
    def _parse_MD_tag(self,read):
        #create a dict of SAM tags
        tags_dict = dict(read.tags)
        #get the MD tag value
        MD_tag    = tags_dict.get('MD',None)
        if MD_tag is None:
            return None,None
        
        #MD tag regex matching
        MD_tag_matches = self.MD_pattern_compiled.split(MD_tag)
        if MD_tag_matches is None:
            print 'No match found %s' % read.qname
            return None,None
        
        #reverse the CIGAR and MD tag for correct intepretation
        if read.is_reverse:
            cigar = read.cigar[::-1]
            split_MD_tag = MD_tag_matches[::-1]
        else:
            cigar = read.cigar
            split_MD_tag = MD_tag_matches
        return cigar,split_MD_tag
        


    
class mappingSummary(object):
    """
    given a bam file generate the alignment quality stats
    
    #########
    example:
    ##########
    
    #instantiate
    bamSummary = mappingSummary(bamFile,max_readlen=150,quality_offset=33,minQual=30)
    
    #evaluate the indels, mismatches, N%
    # returns a pandas df
    df = bamSummary.eval_alignQ()


    """
    
    #initializer
    def __init__(self,bam,max_readlen,quality_offset=33,minQual=20):
        self.bam = bam
        self.quality_offset = quality_offset
        self.max_readlen = max_readlen
        self.minQual = minQual
        
        #stats with all the reads
        self.mismatches = np.zeros(max_readlen)
        self.indels = np.zeros(max_readlen)
        self.qual = np.zeros(max_readlen)
        
        #stats with reads where no N is found
        self.mismatches_rN = np.zeros(max_readlen)
        self.indels_rN = np.zeros(max_readlen)
        self.qual_rN = np.zeros(max_readlen)
        
        self.Ns     =  np.zeros(max_readlen)    
        self.num_reads = 0
        self.num_reads_rN = 0
        self.unmapped_reads = 0
        self.Ns_unmapped_reads = np.zeros(max_readlen)
        
        #compile the patternt for matching the MD
        self.MD_pattern = '([0-9]+)([A-Z]|\^[A-Z]+)*'
        self.MD_pattern_compiled = re.compile(self.MD_pattern)
        
    
    
    def eval_alignQ(self):
        """
            #main looper
            #which iterates over bam file
        """
        bam_fh=pysam.Samfile(self.bam,'rb')
        
        for count,read in enumerate(bam_fh):
            
            
            if count == 5:
                pass
                #sys.exit()
                
                
            self.num_reads += 1
            read_N_seen = self._check_Ns_across_read(read)
            if read_N_seen is False:
                self.num_reads_rN += 1
            self._check_quality_across_read(read,read_N_seen)
            self._check_mismatch_indels_across_read(read,read_N_seen)
            
            if self.num_reads % 100000 == 0 and self.num_reads > 0:
                print 'processed %d reads' % self.num_reads
        print 'In total processed %d reads' % self.num_reads
        
        #create a summary df to be returned
        df = pandas.DataFrame( { 'N_percent'                     :     (self.Ns/self.num_reads)*100,
                                 'N_percent_unmapped_reads'      :     (self.Ns_unmapped_reads/self.unmapped_reads)*100, 
                                 'mismatches'     :     (self.mismatches/self.num_reads)*100,
                                 'indels'         :     (self.indels/self.num_reads)*100,
                                 'bases_<_Q30'    :     (self.qual/self.num_reads)*100,
                                 'mismatches_rN'  :     (self.mismatches_rN/self.num_reads_rN)*100,
                                 'indels_rN'      :     (self.indels_rN/self.num_reads_rN)*100,
                                 'bases_<_Q30_rN' :     (self.qual_rN/self.num_reads_rN)*100
                               }
                              )
        #store the db for later access
        self.results = df
        return df
        
    
    
    
    def plot_alignment_quality(self):
        """
            plot the results
        """
        
    
        df = self.results
        fig = plt.figure(figsize=(15,10))
        fig.subplots_adjust(hspace=.2,wspace=.3)
        fig.suptitle('Alignment Quality Summary:  %s ' % os.path.basename(self.bam) ,size=14)
        
        ax1 = fig.add_subplot(221)
        df[df.index <= self.max_readlen - 20][['N_percent','mismatches','mismatches_rN']].plot(ax=ax1)
        ax1.set_xlabel('read position')
        ax1.set_ylabel('Percent')
        ax1.set_title('Mismatches v/s Ns percent per read bp position')
        
        ax2 = fig.add_subplot(222)
        df[df.index <= self.max_readlen - 20][['indels','indels_rN']].plot(ax=ax2)
        ax2.set_xlabel('read position')
        ax2.set_ylabel('Percent')
        ax2.set_title('Indels v/s Ns percent per read bp position')
            
        ax3 = fig.add_subplot(223)
        df[df.index <= self.max_readlen - 20][['bases_<_Q30','bases_<_Q30_rN']].plot(ax=ax3)
        ax3.set_xlabel('read position')
        ax3.set_ylabel('Percentage  Bases < Q30')
        ax3.set_title('percent bases < Q30 per cycle')
        
        ax4 = fig.add_subplot(224)
        df[df.index <= self.max_readlen - 20][['N_percent_unmapped_reads']].plot(ax=ax4)
        ax4.set_xlabel('read position')
        ax4.set_ylabel('Percent')
        ax4.set_title('UNMAPPED reads N percent per read bp position')
        
        outplot_name = self.bam + '_alignQ.png'
        plt.savefig(outplot_name)
        
        
    
    
    def _check_mismatch_indels_across_read(self,read,read_N_seen):
        if read.is_unmapped:
            return None
        
        #get the cigar and split MD tag
        cigar,split_MD_tag = self._parse_MD_tag(read)
           
        #jumper for MD tag
        CIAGR_jumper_for_MD = 0
    
        #make sense of CIGAR string
        position_cigar = 0
        for count,(cigar_code,matches) in enumerate(cigar):
            """    
            CIGAR CODE in pysam parsed info
            M    BAM_CMATCH        0
            I    BAM_CINS          1
            D    BAM_CDEL          2
            N    BAM_CREF_SKIP     3
            S    BAM_CSOFT_CLIP    4
            H    BAM_CHARD_CLIP    5
            P    BAM_CPAD          6
            =    BAM_CEQUAL        7
            X    BAM_CDIFF         8
            """
        
            #matches or mismatches
            if cigar_code == 0:
                position_cigar += matches
            
            #hard or soft clipping
            elif cigar_code == 4 or cigar_code == 5:
                self.mismatches[position_cigar:position_cigar+matches] += 1
                if read_N_seen is False:
                    self.mismatches_rN[position_cigar:position_cigar+matches] += 1
                position_cigar += matches
                if count == 0: #if the first item in the CIGAR is soft or hard clipping
                    #calculating the jump ahead for MD tag
                    #find any soft or hard clippings in the beginning of CIGAR
                    #unfortunately they are not reported in MD and better be calculated here
                    CIAGR_jumper_for_MD = matches
            
            #insertions in the read wrt to the reference genome
            elif cigar_code == 1:
                self.indels[position_cigar] += 1
                if read_N_seen is False:
                    self.indels_rN[position_cigar] += 1
                position_cigar += matches
            
            #deletion in the read wrt to the reference genome :
            elif cigar_code == 2:
                self.indels[position_cigar] += 1
                if read_N_seen is False:
                    self.indels_rN[position_cigar] += 1
                #cigar position doesnt change in a deletion as the read length stays same
                #the deletion is in the read wrt to the reference
                #position_cigar += matches
                
                #Example
                #   AGTGATGGG---------GGGGTTCCAGGTGGAGACGAGGACTCC   ( seqeunced read)
                #     --     ---------
                #   AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC  (reference genome)
            else:
                print 'New CIGAR code seen: %d' % cigar_code
                print 'PROBS :%s, MD:%s , CIGAR:%s, position:%d' % (read.qname,MD_tag,read.cigar,position_cigar)
                sys.exit(1)
        
        
        #make sense of MD tag
        #mainly for mismatches
        position_MD = 0
        #jump ahead if there were soft/hard clipping in CIGAR string
        position_MD += CIAGR_jumper_for_MD
        
        for x in split_MD_tag:
            if x == '' or x is None:
                continue
            #mismatches
            if x in ['A','G','C','T','N'] and not x.startswith('^'):  #^ indicate deletion in MD string
                self.mismatches[position_MD:position_MD+len(x)] += 1
                if read_N_seen is False:
                    self.mismatches_rN[position_MD:position_MD+len(x)] += 1
                position_MD += len(x) #now increase the position of MD
            
            #see if string is a integer AND move ahead the position in string
            try:
                x=int(x)
                position_MD += x
            except ValueError:
                pass
    
    
    def _check_Ns_across_read(self,read):
        if read.is_reverse:
            seq = read.seq[::-1]
        else:
            seq = read.seq
        read_positions_with_N =[matches.start() for matches in re.finditer('N',seq)]
        if read_positions_with_N:
            self.Ns[read_positions_with_N] += 1
            #unmapped reads stats
            if read.is_unmapped:
                self.unmapped_reads += 1
                self.Ns_unmapped_reads[read_positions_with_N] += 1
            return True
        else:
            return False
        
    def _check_quality_across_read(self,read,read_N_seen):
        #get the indices of the read position with quality score < min_quality
        if read.is_reverse:
            qual = read.qual[::-1]
        else:
            qual = read.qual
            
        indices= [ count for count,q in enumerate(qual) if (ord(q)- self.quality_offset) < self.minQual ]
        self.qual[indices] += 1
        
        if read_N_seen is False:
            self.qual_rN[indices] += 1
    
    
    def _parse_MD_tag(self,read):
        #create a dict of SAM tags
        tags_dict = dict(read.tags)
        #get the MD tag value
        MD_tag    = tags_dict.get('MD',None)
        if MD_tag is None:
            return None,None
        
        #MD tag regex matching
        MD_tag_matches = self.MD_pattern_compiled.split(MD_tag)
        if MD_tag_matches is None:
            print 'No match found %s' % read.qname
            return None,None
        
        #reverse the CIGAR and MD tag for correct intepretation
        if read.is_reverse:
            cigar = read.cigar[::-1]
            split_MD_tag = MD_tag_matches[::-1]
        else:
            cigar = read.cigar
            split_MD_tag = MD_tag_matches
        return cigar,split_MD_tag
        
    

    
        
        
        
        
        
        
        
        
        
        
        
            
        
    
    
    











"""




# 
#read1.flag =  1 + 2 + 16 + 64         # ( read is paired, properly aligned, rev_comp, first in pair);
#read2.flag =  1 + 2 + 32 + 128        # ( read is paired, properly aligned, mate_is_rev_comp, second in pair);
#
#read1.flag =  1 + 2 + 16 + 32 + 64         # ( read is paired, properly aligned, rev_comp, mate_is_rev_comp, first in pair);
#read2.flag =  1 + 2 + 16 + 32 + 128        # ( read is paired, properly aligned, rev_comp, mate_is_rev_comp, second in pair);
#
#read1.flag =  1 + 2 + 32 + 64         # ( read is paired, properly aligned, first in pair);
#read2.flag =  1 + 2 + 16 + 128        # ( read is paired, properly aligned, second in pair);
#        
#read1.flag =  1 + 2 + 32 + 64         # ( read is paired, properly aligned, mate_is_rev_comp, first in pair);
#read2.flag =  1 + 2 + 16 + 128        # ( read is paired, properly aligned, rev_comp , second in pair);










def split_bam_into_inter_N_intra_chromosomal_bam(bamFile):
    '''
    given a bam file the method will split it into two bins:
        1. read pairs that are mapping to two different chromosomes : chimeras
        2. read pairs that map on the same chromosome
        3. Not sure what needs to be done about the singletons and unmapped pairs : for now they will be removed
    '''
    
    #check if bam is indexed else index it
    bamUtils.index_bam(bamFile)
    
    #bamFile handle
    bamFile_handle = pysam.Samfile(bamFile,'rb')
    
    chimeric_pairs = 0
    proper_mapped_pairs = 0
    unmapped_pairs = 0
    singleton_pairs = 0
    
    #name sort the file(done just once)
    nameSorted_bam_file = bamUtils.name_sort_bam(bamFile)
    nameSorted_bamFile_handle = pysam.Samfile(nameSorted_bam_file,'rb'
    
    
    for count,read1,read2 in itertools.izip(itertools.count(1),nameSorted_bamFile_handle,nameSorted_bamFile_handle):
        
        if count % 1000000 == 0:
            print 'Processed %d pairs ' % count
            
#        if count >= 10000000:
#            break


    
"""
    
    



"""
def bam_to_fastq(bamFile,type=None):
    '''
    given a bam file and type argument(mapped,unmapped,all)
    create a paired fastq
    '''
       
    SUB ='split_bam_by_chromosomes'
    
    print 'pysam version available %s' % pysam.__version__
        
    if type is None:
        print "[%s] : Need to know whether to create fastq for unmapped/mapped/all reads"
        sys.exit()
    
    bamFile_bai = bamFile + '.bai'
    if not os.path.exists(bamFile_bai):
        bamUtils.index_bam(bamFile)
    
    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    file_prefix=os.path.basename(bamFile).replace('.bam','')
    
     
    #open the bam file
    in_bamFile_handle = pysam.Samfile(bamFile,'rb')
    references = in_bamFile_handle.references
    
    # initialize the container to hold the bam file names splitted by chromosome
    list_of_split_bam_files = []
    list_of_chromosomes = []
    for chr in references:
        
        #forming the bam file name for this chr
        chr_bam = base_path_split_files + '/' + file_prefix + '_' + chr + '.bam'
        
        #check if chr_bam exists already, if so just store the name and continue
        if os.path.exists(chr_bam):
            print '[%s]: for %s bam %s exists' % (SUB, chr, chr_bam)
            list_of_split_bam_files.append(chr_bam)
            list_of_chromosomes.append(chr) 
        else:
            chr_bam_handle = pysam.Samfile(chr_bam,'wb', template=in_bamFile_handle)
            
            for read in in_bamFile_handle.fetch(chr):
                if read:
                    chr_bam_handle.write(read)
                    
            print '[%s]: created bam %s ' % (SUB,chr_bam)
            list_of_split_bam_files.append(chr_bam)
            list_of_chromosomes.append(chr)
        
    return (list_of_split_bam_files,list_of_chromosomes)

"""
