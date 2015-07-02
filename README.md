# HiC-align
Suite of scripts to align and count HiC reads

- Data is on HAL in shared directory for our lab
- Mark sent me his filter, count scripts.  Added them to the repo.

# Trying all this out using [hiclib](http://mirnylab.bitbucket.org/hiclib/)

Steps from start to finish using the hiclib based python code:

1. Run `align_reads_hiclib.py` to align the reads in fastq files to the genome using bowtie.  The command I used is :

`python align_reads_hiclib.py <fastqdir> <hg19 dir> <bowtie index file prefix> <sam file output dir>`  

2. Run `generate_datasets_file.py` to generate a tab separated values list of .hdf5 chunks (containing bam files) for mapping
reads to bins for all the files in the experiments.  The command I used is:

 `python generate_datasets_file.py -b ~/projects/HiC-align/data/ -r runs.tsv -d datasets.tsv`

3. Filtering some reads TODO

4. Mapping reads, assembling heatmaps TODO

5. Save heatmaps, one per chromosome per replicate TODO
