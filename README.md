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

3. Run `filter_and_assemble_heatmaps.py` to filter some reads due to poor mapability or length (<10bp, >10kb),
groupp the mapped reads into bins (at different resolutions), and output heatmaps. See the script for details.  
This generates whole genome heatmaps at 200kb, 100kb resolutions, and per chromosome heatmaps at 100kb resolution.  
'Resolution' here means mean bin size, I think.  The command I used is:

`python filter_and_assemble_heatmaps.py ../data/dataset.tsv`

4. Run `plot_diagonal_correlation(200000)` from `plot_heatmaps.py` to visualize the Spearman correlation between both replicates 
at different genomic distances.  
