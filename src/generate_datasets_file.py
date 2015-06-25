from mirnylib.h5dict import h5dict
#import pandas as pd
#from mirnylib.systemutils import setExceptionHook
import os
from argparse import ArgumentParser

"""
Adapted from the hiclib/examples/pipeline2015/01_makeDatasetsFile.py

The goal of this file is to create a datasets.tsv file.

The file has the following structure:
Filename                                Experiment  Replicate   Genome  RestrictionEnzyme

mapped-hg19/SRR3141592/chunk0001.hdf5   K562        R1          hg19    HindIII
mapped-hg19/SRR3141592/chunk0002.hdf5   K562        R1          hg19    HindIII
mapped-hg19/SRR3141593/chunk0001.hdf5   K562        R1          hg19    HindIII
mapped-hg19/SRR3141593/chunk0002.hdf5   K562        R1          hg19    HindIII
mapped-hg19/SRR3141594/chunk0001.hdf5   K565        R2          hg19    HindIII

"""

parser = ArgumentParser()
parser.add_argument("-b", "--basedir", help="the base directory for the data in the project")
parser.add_argument("-r", "--runsfile", help="the runs.tsv file")
parser.add_argument("-d", "--datafile", help="the dataset file to be output (default name is datasets.tsv, same dir as runs file)")
args = parser.parse_args()

# open and parse the runs file
runs_file = open(os.path.join(args.basedir,args.runsfile),"r")
runs = [run.split() for run in runs_file.readlines() if not run.startswith("#")]
runs_file.close()

# process each record in the runs file, write out to the data sets file
datasets_file = open(os.path.join(args.basedir,args.datafile),"w")

# print header for datasets file
datasets_file.write("# The file has the following structure:\n")
datasets_file.write("# Filename\tExperiment\tReplicate\tGenome\tRestrictionEnzyme\n")

for run in runs:
    input_dir, experiment, replicate, genome, restriction_enzyme = run
    filenames = [j for j in os.listdir(os.path.join(args.basedir,input_dir)) if j.endswith(".hdf5") ]
    for fname in filenames:
        try:
            mydict = h5dict(os.path.join(args.basedir,input_dir,fname),'r')
        except:
            pass
        if "strands1" not in mydict:
            raise
        if len(mydict.get_dataset("strands1")) < 10000:
            raise
        datasets_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(os.path.join(input_dir,fname),experiment,
                                                    replicate, genome, restriction_enzyme))

datasets_file.close()