import os
import loging
import argparse
import subprocess
import glob
from joblib import Parallel, delayed
import pdb

try:
	from hiclib import mapping
	from mirnylib import h5dict, genome
except ImportError:
	print "Missing either hiclib or mirnylib imports, exiting..."
	os._exit()


# with the params established ^, calculate the alignment step size increment
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

# quick sanity check on a given fastq file.  The second line should contain a non-zero length
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
	min_len, step_size = calculate_step(length - seq_skip_start, min_map_len)
	first_sam = first_fq.split(".fastq.gz")[0] + ".sam"
	second_sam = second_fq.split(".fastq.gz")[0] + ".sam"

	# map the fastq files -> sam files
	mapping.iterative_mapping(
	bowtie_path = bowtie_path,
	bowtie_index_path = bowtie_index,
	fastq_path=first_fq,
	out_sam_path=os.path.join(args.samdir,first_sam),
	min_seq_len=min_len,
	len_steps=step_size,
	seq_start=seq_skip_start,
	nthreads=threads,
	bowtie_flags=bowtie_flags)

	mapping.iterative_mapping(
	bowtie_path = bowtie_path,
	bowtie_index_path = bowtie_index,
	fastq_path=second_fq,
	out_sam_path=os.path.join(args.samdir,second_sam),
	min_seq_len=min_len,
	len_steps=step_size,
	seq_start=seq_skip_start,
	nthreads=threads,
	bowtie_flags=bowtie_flags)

	# parse the mapped sequences into a Python data structure,
	# assign the ultra-sonic fragments to restriction fragments. <- what the hell does this even mean?
	mapped_reads = h5dict.h5dict(outfile)
	sf1, sf2 = [os.path.join(args.samdir,first_sam), os.path.join(args.samdir,second_sam)]
	mapping.parse_sam(sam_basename1=sf1, sam_basename2=sf2,
	out_dict=mapped_reads, genome_db=genome_db, save_seqs=False, maxReads=10000000, IDLen=50)


if __name__ == "__main__":

	logging.basicConfig(level=logging.DEBUG)

	parser = argparse.ArgumentParser()

	parser.add_argument("fastqdir", help="parse the fastq files in this dir")
	parser.add_argument("hg19", default="/home/zamparol/data/hg19/" help="reference genome is located in this dir")
	parser.add_argument("index", help="the bowtie2 index for hg19 is here")
	parser.add_argument("samdir", help="put the sam output in this dir")
	args = parser.parse_args()

	# find bowtie2, and genome index
	bowtie_path = subprocess.check_output("which bowtie2", shell=True)
	bowtie_path = bowtie_path.strip()
	bowtie_index = args.index # change this if your index is named differently from the genome
	bowtie_flags = "--very-sensitive"

	# number of threads for bowtie alignment
	threads = 4

	pdb.set_trace()

	# scratch space for temporary files during mapping, binning
	tmp_dir = "/home/zamparol/scratch"  

	# they need some object to represent hg19
	genome_db = genome.Genome(args.hg19, readChrms=["#", "X", "Y"], cacheDir="tmp_dir")
	genome_name = genome_db.folderName  # automatically infer genome name from the folder name. Name is used for naming folders ("mapped-hg19", etc)

	seq_skip_start = 2  # skip first 2 bp of the read, if you want
	min_map_len = 25  # start mapping at this length

	# go to the fastq dir, read both pair files and compose the sam file output names
	os.chdir(args.fastqdir)
	first_ends = glob.glob('*_R1_*.fastq.gz')
	second_ends = glob.glob('*_R2_*.fastq.gz')
	first_ends.sort()
	second_ends.sort()
	
	outfiles = [f.split("_R1_")[0] + (f.split("_R1")[1]).split(".fastq.gz")[0] + ".hdf5" for f in r_one]

	# map the reads in parallel:
	Parallel(n_jobs=8)(delayed(map_reads)(p,q,outf) for p,q,outf in zip(first_ends,second_ends,outfiles))
	


# This is a current version of the pipeline which we use for processing many Hi-C datasets. This version is for .sra files or fastq.gz files split by two sides.

# Usage:

# 1. place .sra files in the fastq folder
# 2. edit 00_mapData.py. Specify genome and paths to bowtie, and other arguments
# 3. provide runs.tsv file as described in 01_makeDatasetsFile.py
# 4. possibly edit 02_mergeDatasets.py to specify resolutions (especially when working with smaller genomes)
# 5. run 02_mergeDatasets.py   -  it will do several things.
#     -- merge (if necessary) different .hdf5 files corresponding to one replica of one datasets
#     -- do fragment-level filtering of merged files
#     -- save heatmaps at different resolutions
#     -- merge files from different replicas of the same experiment
#     -- save heatmaps