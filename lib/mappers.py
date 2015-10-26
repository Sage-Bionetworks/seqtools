#!/usr/bin/env python


import pandas


#local modules
import utils
import bamUtils


def run_tophat2(args):

	read1 = ','.join(args.read1)
	read2 = ','.join(args.read2)

	tophat2_args = ['tophat2', '-p', args.cores, '--output-dir', args.outdir+'/'+args.prefix]
	if args.outdir is None:
		'tophatmapped_%s_out' % time.strftime('%Y-%m-%d:%H-%M-%S')
	if args.ref_genes is not None:
		tophat2_args.extend(['-G',args.ref_genes])
	if args.optional_args is not None:
		tophat2_args.extend([args.optional_args])

	tophat2_args.extend([args.genome_index,read1, read2 ])

	#run the external command
	utils.run_external_command(tophat2_args, logger=args.logger)



def read_tophat_alignSummary(align_summary_files):
	
	overall_mapping_percent = []
	summary_files = []
	for summary_file in align_summary_files:
		with open(summary_file) as f:
			for line in f:
				m = re.search('^(.+?)% overall.*$', line)
				if m:
					overall_mapping_percent.append(float(m.group(1)))
					summary_files.append(summary_file)
	df = pd.DataFrame({'sample':summary_files, 'map_percent':overall_mapping_percent})
	df = df.set_index(['sample'])
	return df


def get_star_alignSummary(alignment_summary_files):
	
	star_stats = []
	for alignment_summary_file in alignment_summary_files:
		df = pandas.read_csv(alignment_summary_file,sep="\t",names=['stats',alignment_summary_file])
		df.stats = df.stats.map(lambda x: str(x).strip().replace(' |',''))
		df = df.set_index(['stats'])
		star_stats.append(df)
	star_stats = pandas.concat(star_stats, axis=1)

	uniq_reads_percent = df.ix['Uniquely mapped reads %'].map(lambda x: x.replace('%','')).astype(float)
	multimapped_reads_percent = df.ix['% of reads mapped to multiple loci'].map(lambda x: x.replace('%','')).astype(float)
	total_mapped_reads = uniq_reads_percent + multimapped_reads_percent

	return star_stats
		


def run_STAR(args):

	read1 = ','.join(args.read1)
	read2 = ','.join(args.read2)

	output_dir =  args.outdir + '/' + args.prefix + '/'
	utils.create_dir(output_dir)

	STAR_args = ['STAR', '--runThreadN', args.cores, '--genomeDir', args.genome_index,
				 '--outFileNamePrefix ',output_dir , '--outSAMunmapped', 'Within']
	if args.outdir is None:
		args.outdir = 'STARmapped_%s_out' % time.strftime('%Y-%m-%d:%H-%M-%S')
	if args.optional_args is not None:
		STAR_args.extend([args.optional_args])
	STAR_args.extend(['--readFilesIn',read1, read2 ])
	
	#to make sure we convert all the args to string
	STAR_args = map(str,STAR_args)
	
	#run STAR
	utils.run_external_command(STAR_args, logger=args.logger)

	#expected output sam file
	expected_sam_file =  output_dir + 'Aligned.out.sam'

	#convert sam to bam
	bamFile = bamUtils.sam_to_bam(expected_sam_file, delete_sam=True)
