import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from mirnylib import plotting
from mirnylib.plotting import removeAxes

from hiclib import fragmentHiC
from hiclib.binnedData import binnedData, binnedDataAnalysis

from scipy.stats import spearmanr as corr


fragments_dir = "/Users/zamparol/projects/HiC-align/results/fragments"
heatmaps_dir = "/Users/zamparol/projects/HiC-align/results/heatmaps/"
hg19_dir = "/Users/zamparol/projects/HiC-align/data/hg19"

'''Plots figure with correlation at different binning.
Note the caching and creating of binned heatmaps flags below.
Suppplementary paper figure
'''

# load up both reps
rep1_frags = fragmentHiC.HiCdataset("throw_away", hg19_dir, enzymeName="HindIII",
                              override=False, inMemory=True)
rep1_frags.load(os.path.join(fragments_dir, "Martin-R1-HindIII_refined.frag"))

rep2_frags = fragmentHiC.HiCdataset("throw_away", hg19_dir, enzymeName="HindIII",
                              override=False, inMemory=True)
rep2_frags.load(os.path.join(fragments_dir, "Martin-R2-HindIII_refined.frag"))

heatmap_dict = {100000: ['Martin-R1-HindIII-100k.byChr','Martin-R2-HindIII-100k.byChr'],
                200000: ['Martin-R1-HindIII-200k.hm','Martin-R2-HindIII-200k.hm']}

def plot_diagonal_correlation(resolution):
    '''
    Plot spearman correlation of replicates across different off-diagonal bands
    :param: resolution (int) the resolution for bin construction
    '''
    "Correlation of diagonal bins - paper figure"

    up_to = 50
    bands = np.arange(2, up_to)
    binned_data = binnedData(resolution, hg19_dir)

    # try to load the data
    pair = tuple(heatmap_dict[resolution])
    for hm in heatmap_dict[resolution]:
        binned_data.simpleLoad(os.path.join(heatmaps_dir,hm),hm)

    # remove main diagonal, filter 'poor' regions (??)
    binned_data.removeDiagonal(1)
    binned_data.removePoorRegions()
    binned_data.removeZeros()
    cors = []
    for i in bands:
        cors.append(corr(
                       np.diagonal(binned_data.dataDict[pair[0]], i),
                       np.diagonal(binned_data.dataDict[pair[1]], i)
                       )[0])


    # binned_data.iterativeCorrectWithoutSS(M=1)
    # cors2 = [[] for _ in pairs]
    # for i in bands:
    #     for j, pair in enumerate(pairs):
    #         cors2[j].append(cr(
    #                         numpy.diagonal(binned_data.dataDict[pair[0]], i),
    #                         numpy.diagonal(binned_data.dataDict[pair[1]], i)
    #                         )[0])
    #
    # binned_data.iterativeCorrectWithoutSS(M=20)
    # cors3 = [[] for _ in pairs]
    # for i in bands:
    #     for j, pair in enumerate(pairs):
    #         cors3[j].append(cr(
    #                         numpy.diagonal(binned_data.dataDict[pair[0]], i),
    #                         numpy.diagonal(binned_data.dataDict[pair[1]], i)
    #                         )[0])

    matplotlib.rcParams['font.sans-serif'] = 'Arial'

    #plt.figure(figsize = (2.3,1.8))
    print cors
    # print cors2
    # print cors3
    plt.figure(figsize=(10, 3))
    ax = plt.gca()

    #plt.subplot(1, len(pairs), j)
    fs = 8
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    plt.title("%s vs %s" % pair)
    # plt.plot(bands / 5., cors3[j], color="#E5A826", label="Iterative")
    # plt.plot(bands / 5., cors2[j], color="#28459A", label="Single")
    plt.plot(bands / 10., cors, color="#E55726", label="Raw")
    plt.xlabel("Genomic Separation, MB", fontsize=8)
    plt.ylabel("Spearman correlation", fontsize=8)
    plt.legend()

    legend = plt.legend(prop={"size": 6}, loc=9, handlelength=2)
    legend.draw_frame(False)
    plt.ylim((0, 1))
    removeAxes(shift=0)

    plt.show()


# # --------Filter only trans DS reads-----------------
# rep1.maskFilter(rep1.DS * (rep1.chrms1 != rep1.chrms2))
# rep2.maskFilter(rep2.DS * (rep2.chrms1 != rep2.chrms2))
#
# # Now create two halfs of one dataset and down-sample second dataset
# # ----------------------standard version code--------
# fraction = 0.5 * len(rep1.DS) / float(len(FR2.DS))
#
# rarray = numpy.random.random(len(rep1.DS))
# mask1 = rarray < 0.5
# mask3 = rarray >= 0.5
#
# rep1.maskFilter(mask1)
# rep2.maskFilter(mask3)
# rep1.save("../../../tcc/working/cache1")
# rep2.save("../../../tcc/working/cache3")

# p1 = []
# p2 = []
# p3 = []
# p4 = []
# evs = []
#
# for size in resolutions:
#
#     BD = binnedDataAnalysis(size, hg19_dir)
#     BD.simpleLoad("../../../tcc/working/HindIII_%d.hm" % size, "HindIII")
#     BD.simpleLoad("../../../tcc/working/NcoI_%d.hm" % size, "NcoI")
#     BD.simpleLoad("../../../tcc/working/control_%d.hm" % size, "control")
#     BD.removeDiagonal()
#     BD.removePoorRegions(cutoff=2)
#     BD.removeCis()
#
#     data1 = BD.dataDict["HindIII"]
#     data2 = BD.dataDict["NcoI"]
#     data3 = BD.dataDict["control"]
#
#     mask = (numpy.sum(
#         data1, axis=0) > 0) * (numpy.sum(data2, axis=0) > 0)
#     validMask = mask[:, None] * mask[None, :]
#     transmask = BD.chromosomeIndex[:, None] != BD.chromosomeIndex[None, :]
#     cormask = transmask * validMask
#
#     c1 = scipy.stats.spearmanr(data1[cormask], data2[cormask])[0]
#     c4 = scipy.stats.spearmanr(data1[cormask], data3[cormask])[0]
#
#     if size == 1:
#         evs.append(BD.interchromosomalValues("HindIII"))
#         evs.append(BD.interchromosomalValues("NcoI"))
#         evs.append(BD.interchromosomalValues("control"))
#     p4.append(c4)
#     p1.append(c1)
#
#     print "size\t%d\traw:" % size, c1,
#     BD.removeZeros()
#     BD.fakeCis()  # does iterative correction as well
#     BD.restoreZeros(value=0)
#
#     data1 = BD.dataDict["HindIII"]
#     data2 = BD.dataDict["NcoI"]
#     data3 = BD.dataDict["control"]
#     c2 = scipy.stats.spearmanr(data1[cormask], data2[cormask])[0]
#     c3 = scipy.stats.spearmanr(data1[cormask], data3[cormask])[0]
#
#     if size == 1:
#         evs.append(BD.interchromosomalValues("HindIII"))
#         evs.append(BD.interchromosomalValues("NcoI"))
#         evs.append(BD.interchromosomalValues("control"))
#         print evs
#
#     p3.append(c3)
#     p2.append(c2)
#
#     print "\tcorrected:", c2, "\tcontrol", c3
#
# plt.plot(resolutions, p1, label="Raw data, between enzymes")
# plt.plot(resolutions, p2, label="Iteratively corrected, between")
# plt.plot(resolutions, p3, label="IC, within")
# plt.xlabel("Bin size, MB")
# plt.xticks(range(1, 11))
# plt.ylabel("Spearman correlation coefficient")
# plt.legend()
# plotting.niceShow()

