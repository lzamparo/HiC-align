import os
import logging
import argparse
import subprocess
import glob
import gzip
import numpy as np

from joblib import Parallel, delayed

try:
	from hiclib import mapping
	from mirnylib import h5dict, genome
except ImportError:
	print "Missing either hiclib or mirnylib imports, exiting..."
	os._exit()


# calculate the alignment step-size increment
def calculate_step(length, minlen, approxStep=10, maxSteps=4):
    """returns minimum length and step based on the
    length of sequence and proposed minimum length"""

    actualDif = length - minlen
    if actualDif < approxStep * 0.6:
        return length, 100

    numIter = np.array(np.around(actualDif / float(approxStep)), dtype=int)
    if numIter == 0:
        numIter = 1
    if numIter > maxSteps:
        numIter = maxSteps
    actualStep = actualDif / numIter

    minlen = length - actualStep * numIter

    return minlen, actualStep

# quick sanity check on a given fastq file.  The second line should contain 
# a sequence of length greater than 10
def check_len(fastq_file):
	f = gzip.open(fastq_file, 'r')
	f.readline()
	length = len(f.readline()) - 1
	if length < 10:
		raise ValueError("Length of your sequence is {0}. Something is wrong ".format(length))
	else:
		return length



# map the reads from each paired end reads fastq file (first,second) to 
# corresponding .sam files, then store both mapped sam files in an hdf5 file (outfile)
def map_reads(first_fq,second_fq,outfile):
	first_sam = first_fq.split(".fastq.gz")[0] + ".sam"
	second_sam = second_fq.split(".fastq.gz")[0] + ".sam"

	# map the first fastq file -> sam file
	length = check_len(first_fq)
	min_len, step_size = calculate_step(length - seq_skip_start, min_map_len)
	mapping.iterative_mapping(bowtie_path = bowtie_path,
		bowtie_index_path = bowtie_index,
		fastq_path=first_fq,
		out_sam_path=os.path.join(args.samdir,first_sam),
		min_seq_len=min_len,
		len_step=step_size,
		seq_start=seq_skip_start,
		nthreads=threads,
		bowtie_flags=bowtie_flags)

	# map the second fastq file -> sam file
	length = check_len(second_fq)
	min_len, step_size = calculate_step(length - seq_skip_start, min_map_len)
	mapping.iterative_mapping(
	bowtie_path = bowtie_path,
	bowtie_index_path = bowtie_index,
	fastq_path=second_fq,
	out_sam_path=os.path.join(args.samdir,second_sam),
	min_seq_len=min_len,
	len_step=step_size,
	seq_start=seq_skip_start,
	nthreads=threads,
	bowtie_flags=bowtie_flags)

	# parse the mapped sequences into a the hdf5 dict structure,
	# assign the ultra-sonic fragments to restriction fragments. <- what the hell does this even mean?
	out_dict = os.path.join(args.samdir,outfile)
	mapped_reads = h5dict.h5dict(out_dict)
	sf1, sf2 = [os.path.join(args.samdir,first_sam), os.path.join(args.samdir,second_sam)]
	mapping.parse_sam(sam_basename1=sf1, sam_basename2=sf2,
	out_dict=mapped_reads, genome_db=genome_db, save_seqs=False, maxReads=10000000, IDLen=50, enzyme_name='HindIII')


if __name__ == "__main__":

	logging.basicConfig(level=logging.DEBUG)

	parser = argparse.ArgumentParser()

	parser.add_argument("fastqdir", help="parse the fastq files in this dir")
	parser.add_argument("hg19", default="/home/zamparol/data/hg19/", help="reference genome is located in this dir")
	parser.add_argument("index", help="the bowtie2 index for hg19 is here")
	parser.add_argument("samdir", help="put the sam output in this dir")
	args = parser.parse_args()

	# find bowtie2, and genome index
	bowtie_path = subprocess.check_output("which bowtie2", shell=True)
	bowtie_path = bowtie_path.strip()
	
	# Their bowtieIndex seems finicky
	# e.g bowtieIndex = "../bin/bowtie2/index/{0}".format(genomeName)  # change this if your index is named differently from the genome

	bowtie_index = args.index # change this if your index is named differently from the genome
	bowtie_flags = "--very-sensitive"

	# number of threads for bowtie alignment
	threads = 4

	# scratch space for temporary files during mapping, binning
	tmp_dir = "/home/zamparol/scratch"  

	# they need some object to represent hg19
	genome_db = genome.Genome(args.hg19, readChrms=["#", "X", "Y"], cacheDir=tmp_dir)
	genome_name = genome_db.folderName  # automatically infer genome name from the folder name. Name is used for naming folders ("mapped-hg19", etc)

	seq_skip_start = 2  # skip first 2 bp of the read, if you want
	min_map_len = 25  # start mapping at this length

	# go to the fastq dir, read both pair files and compose the sam file output names
	os.chdir(args.fastqdir)
	first_ends = glob.glob('*_R1_*.fastq.gz')
	second_ends = glob.glob('*_R2_*.fastq.gz')
	first_ends.sort()
	second_ends.sort()
	
	outfiles = [f.split("_R1_")[0] + (f.split("_R1")[1]).split(".fastq.gz")[0] + ".hdf5" for f in first_ends]

	# map the reads in parallel:
	Parallel(n_jobs=4)(delayed(map_reads)(p,q,outf) for p,q,outf in zip(first_ends,second_ends,outfiles))
