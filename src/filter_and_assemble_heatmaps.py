"""
This example script does the following:

-Loads individual files from the location defined in datasets.tsv file
-Parses individual files in memory  (inMemory = True in "for onename in in_files))
   -If your system does not have enough memory, you might need to switch to hdf5 here.
-Merges files corresponding to the same experiment together, on the HDD.
-Filters datasets, builds heatmaps

I'd like to do all steps save for the last one:

- Combines multiple replicas of the same experiment together, builds heatmaps

--Datasets are defined in the datasets.tsv file

--genome is defined by genomeFolder function, and workingGenome identifyer

--output files are arranged to folders named by their workingGenome IDs


Warnings:
    Running this over NFS might cause unexpected slow-downs because NFS is
    unhappy with repeated read/write cycles to the same file

    You could do entire thing in memory, if you have RAM or your datasets are small.
    Actually, using HDF5 is then equivalent to storing compressed data in RAM,
    and might be in fact pretty fast.
    
    General recommendation: if you have 16+GB of RAM, and .sra (.fastq.gz) files were less than 30GB, then you should be fine with parsing things in memory. 

"""

from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap, setExceptionHook
import numpy as np
import os
from argparse import ArgumentParser

from joblib import Parallel, delayed

setExceptionHook()


def ensure(f):
    d = os.path.dirname(f)
    if os.path.isdir(d):
        return f
    else:
        try:
            os.makedirs(d)
        except:
            raise ValueError("Cannot create directory")
    return f


def genomeFolder(name):
    return os.path.join("/home/zamparol/data", name)  # Fetch genome folder by genome name

### LZ: original script parameters, I'll set my own
whole_genome_resolutions_Kb = [200]  # [2000,1000,500,200]
by_chromosome_resolutions_Kb = [100] # [100,40]
hi_res_with_overlap_resolutions_Kb = []  # [20,10]        #add 10 here if you have more than 16GB RAM
super_hi_res_with_overlap_resolutions_Kb = []  # [5]
skip = 1  # how many to skip for single replica datasets


# wholeGenomeResolutionsKb = [2000,1000,500,200]
# byChromosomeResolutionsKb = [100, 40]
# HiResWithOverlapResolutionsKb = [20,10]     #add 10 here if you have more than 16GB RAM
# SuperHiResWithOverlapResolutionsKb = [5]
# skip = 1                                    #how many to skip for single replica datasets


def parse_mapped_reads(onename, niceval, working_genome, enzyme, stat_folder, in_memory):
    """
    Parse the given h5 mapped reads, output to a partial fragment file.

    :param onename: (string) the name of this fragment
    :param niceval: (int) positive int to rank the priority of this job
    :param working_genome: (string) name of the genome aganist which we mapped the reads
    :param enzyme: (string) name of the restriction enzyme used to cut fragments
    :param stat_folder: (string) folder in which to write fragment mapping stats
    :return: none explicit, h5 file is saved to disk
    """

    # set the niceness of this sub-process:
    os.nice(niceval)

    np.random.seed()
    # Parsing individual files, either in memory or on disk
    if in_memory:
        finalname = onename + "_parsed.frag"
        TR = HiCdataset("bla" + str(np.random.randint(100000000000)), genome=genomeFolder(working_genome),
                        maximumMoleculeLength=500, enzymeName=enzyme, tmpFolder="tmp",
                        inMemory=True)  # remove inMemory if you don't have enough RAM

        TR.parseInputData(dictLike=onename)
        print onename
        TR.save(ensure(finalname))
        folder, fname = os.path.split(onename)
        statSubFolder = os.path.join(stat_folder, folder)

        TR.printMetadata(saveTo=ensure(os.path.join(statSubFolder, fname + ".stat")))
    else:
        # Create dataset at destination, parse on HDD, then no need to save.
        TR = HiCdataset(ensure(onename + "_parsed.frag"),
                        genome=genomeFolder(working_genome), enzymeName=enzyme, tmpFolder="tmp",
                        maximumMoleculeLength=500, mode='w')
        TR.parseInputData(dictLike=onename, enzymeToFillRsites=enzyme)
        TR.printMetadata(saveTo=ensure(os.path.join(stat_folder, onename + ".stat")))




def refine_dataset(filenames, niceval, delete=True, parse_in_memory=True):
    """
    Map the fragments from each replicate to chromosomes (in parallel)

    Parameters
    ----------
    filenames[0] is a list of filenames of incoming files
    filenames[1] is a folder for outgoing file
    filenames[2] is a working genome, that is output directory
    filenames[3] is an enzyme for a given experiment

    create : bool, optional
        If True, parse each file.
        If False, assume that files were already parsed
        (e.g. if you are just playing around with filtering parameters)
    delete : bool, optional
        If True, delete parsed files after merging.
        Man, these files may be huge... if you don't have a 10TB RAID, this may be useful.
    parseInMemory : bool, optional
        Perform parsing input files in memory.
    """

    in_files = filenames[0]
    out_file = filenames[1]

    stat_folder = os.path.join("statistics", out_file)

    working_genome = filenames[2]
    enzyme = filenames[3]

    nice_list = [niceval for i in in_files]
    parse_list = [parse_in_memory for i in in_files]
    genome_list = [working_genome for i in in_files]
    enzyme_list = [enzyme for i in in_files]
    stat_folder_list = [stat_folder for i in in_files]

    #map(parse_onename, in_files)
    Parallel(n_jobs=20)(delayed(parse_mapped_reads)(infile, nice_val, genome, enzyme, stat_folder, parse_val) for infile, nice_val, genome, enzyme, stat_folder, parse_val in
               zip(in_files, nice_list, genome_list, enzyme_list, stat_folder_list, parse_list))

    # Merge in all parsed files from one experiment
    print("Merging files all together, applying filters...")
    TR = HiCdataset(ensure(out_file + "_merged.frag"),
                    genome=genomeFolder(working_genome), enzymeName=enzyme, tmpFolder="tmp", dictToStoreIDs="h5dict",
                    mode="w")
    TR.merge([i + "_parsed.frag" for i in in_files])

    if delete:  # cleaning up parsed files
        for delFile in [i + "_parsed.frag" for i in in_files]:
            os.remove(delFile)
    print("done!")
    print("Filtering merged data...")
    TR = HiCdataset(out_file + "_refined.frag", enzymeName=enzyme,
                    genome=genomeFolder(working_genome), tmpFolder="tmp", dictToStoreIDs="h5dict",
                    mode='w')
    TR.load(out_file + "_merged.frag")


    # ----------------------------Set of filters applied -------------
    TR.filterDuplicates()
    TR.filterLarge(10000, 10)
    TR.filterExtreme(cutH=0.001, cutL=0)
    TR.writeFilteringStats()
    TR.printMetadata(saveTo=stat_folder + ".stat")
    print("done!")
    # ------------------------End set of filters applied----------

    print("Building heatmaps at specified resolutions...")
    TR.printStats()
    for res in whole_genome_resolutions_Kb:
        TR.saveHeatmap(out_file + "-{0}k.hm".format(res), res * 1000)

    for res in by_chromosome_resolutions_Kb:
        TR.saveByChromosomeHeatmap(out_file + "-{0}k.byChr".format(res), res * 1000)

    for res in hi_res_with_overlap_resolutions_Kb:
        TR.saveHiResHeatmapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res * 1000)

    for res in super_hi_res_with_overlap_resolutions_Kb[:-skip]:
        TR.saveSuperHighResMapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res * 1000)
    print("done!")



### Script begins here ###
parser = ArgumentParser()
parser.add_argument("datasets", help="use this datasets .tsv file")
parser.add_argument("-n", "--niceness", default=10, help="nice value for subprocesses that use")
args = parser.parse_args()

fsplit = os.path.split(args.datasets)
if len(fsplit[0]) > 0:
    os.chdir(fsplit[0])
filename = fsplit[1]

# Catalogue all .hdf5 bam file constructs to be binned and filtered
data_file = open(filename, 'r')
data_files = data_file.readlines()
data_files = [line.strip() for line in data_files if not line.startswith("#")]
data_file.close()

# Ensure each line has 5 fields
for line in data_files:
    if len(line.split()) != 5:
        print "incomplete line", line
        raise

# Each experiment is defined by (experiment name, replicate name, genome, enzyme)
experiment_names = {tuple(elem.split()[1:5]) for elem in data_files}
by_experiment = []
combined_experiment_names = []

for experiment in experiment_names:
    # experiment is HeLa,R1,hg19,HindIII
    this_experiment_name, this_replicate, this_genome, this_enzyme = experiment
    # out_name is experiment-replicate-restriction-enzyme
    filenames = [i.split()[0] for i in data_files if tuple(i.split()[1:5]) == experiment]
    out_name = "-".join([this_experiment_name, this_replicate, this_enzyme])

    # filenames, path to save, genome, enzyme
    by_experiment.append((filenames, os.path.join(this_genome, out_name), this_genome, this_enzyme))
    # merged files belonging to one expriment
    combined_experiment_names.append((this_experiment_name, os.path.join(this_genome, out_name), this_genome, this_enzyme))

# map the reads in parallel:
for experiment in by_experiment:
    refine_dataset(experiment, args.niceness)

# TODO: cleaned up to here.  The refineDatasets function *should* produce chromosome by chromosome heatmaps
# that can be examined later (somehow).

# now merge different experiments all together
# note that the first column is not here, as it is a replica

# def extract_combined_experiment_fields(experiment):
#     return ((experiment[0], experiment[2], experiment[3]))
#
# experiments = set([extract_combined_experiment_fields(i) for i in combined_experiment_names])
#
# for experiment in experiments:
#     this_genome = experiment[1]
#     myExperimentNames = [i[1] + "_refined.frag" for i in combined_experiment_names if (i[0], i[2], i[3]) == (experiment[0], experiment[1],experiment[2])]
#     assert len(myExperimentNames) > 0
#     if len(myExperimentNames) > 0:
#         #If we have more than one experiment (replica) for the same data, we can combine.
#         TR = HiCdataset(os.path.join(this_genome, "%s-all-%s_refined.frag" %
#                                      (experiment[0],experiment[2])), genome=genomeFolder(this_genome),
#                                      enzymeName = experiment[2],tmpFolder = "tmp",dictToStoreIDs="h5dict")
#         statSaveName = os.path.join("statistics", this_genome, "%s-all-%s_refined.stat" % (experiment[0], experiment[2]))
#
#         TR.merge(myExperimentNames)
#
#         TR.printMetadata(saveTo=statSaveName)
#
#         for res in whole_genome_resolutions_Kb:
#             TR.saveHeatmap(os.path.join(this_genome, "%s-all-%s-{0}k.hm" % (experiment[0], experiment[2])).format(res), res*1000)
#
#         for res in by_chromosome_resolutions_Kb:
#             TR.saveByChromosomeHeatmap(os.path.join(this_genome, "%s-all-%s-{0}k.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
#
#         for res in hi_res_with_overlap_resolutions_Kb:
#             TR.saveHiResHeatmapWithOverlaps(os.path.join(this_genome, "%s-all-%s-{0}k_HighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
#
#         for res in super_hi_res_with_overlap_resolutions_Kb:
#             TR.saveSuperHighResMapWithOverlaps(os.path.join(this_genome, "%s-all-%s-{0}k_SuperHighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
