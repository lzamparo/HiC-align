library(GenomicAlignments)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)

readfile <- args[1]
numsites <- as.integer(args[2])
outfile <- gsub('SRR','Counts_GM12878_SRR',readfile)
chrom <- sprintf('chr%s',c(1:22,"X","Y"))

hiCreads <- readRDS(readfile)

# Arrange Hi-C paired-end reads by chromosome

R1 <- split(hiCreads[['R1']],seqnames(hiCreads[['R1']]))
R2 <- split(hiCreads[['R2']],seqnames(hiCreads[['R2']]))
ReadsPerChrom <- list()
for(chr in chrom){
  ReadsPerChrom[[chr]] <- list("R1" = R1[[chr]],"R2" = R2[[chr]])
}
rm (R1,R2,hiCreads)
gc()
gc(reset=TRUE)

# We constructive human genome based on enzyme cutting sites

hg19chromSizes <- data.table(chr=chrom,size=seqlengths(Hsapiens)[chrom])
hg19chromGR <- GRanges(seqnames = hg19chromSizes$chr,ranges = IRanges(start = 1,end = hg19chromSizes$size))
enzymeCuts <- lapply(chrom,function(x){GRanges(seqnames=x,ranges=ranges(matchPattern('GATC',Hsapiens[[x]])))})
enzymeCuts <- do.call(c,enzymeCuts)
saveRDS(enzymeCuts,file=sprintf('/home/mark/LiebermanHiC/data/enzymeCuts.rds'))
enzymeBins <- setdiff(hg19chromGR,enzymeCuts)
enzymeBinsTab <- data.table(chr = as.character(seqnames(enzymeBins)),
                            start = as.integer(start(enzymeBins))-2,
                            end = as.integer(end(enzymeBins))+2)
setkey(hg19chromSizes,chr)
enzymeBinsTab$chrSize <- hg19chromSizes[enzymeBinsTab$chr]$size
enzymeBinsTab[start < 1]$start <- enzymeBinsTab[start < 1]$start+2
enzymeBinsTab[end > chrSize]$end <- enzymeBinsTab[end > chrSize]$end-2
enzymeBins <- GRanges(seqnames = enzymeBinsTab$chr,IRanges(start = enzymeBinsTab$start,end = enzymeBinsTab$end))
enzymeBinsTab$name <- sprintf('%s-%s-%s',enzymeBinsTab$chr,enzymeBinsTab$start,enzymeBinsTab$end)
rm (hg19chromSizes,hg19chromGR,enzymeBins)
gc()
gc(reset=TRUE)

# We cluster enzyme cutting sites every ten consecutive bins

enzymeBinsGrouping <- list()
for(i in chrom){
  chrEnzymeBins <- enzymeBinsTab[chr == i][order(start)]
  binNumbers <- 1:ceiling(nrow(chrEnzymeBins)/numsites)
  chrEnzymeBins$binNumb <- rep(binNumbers,each=numsites)[1:nrow(chrEnzymeBins)]
  chrEnzymeBins$binName <- sprintf('%s-%s',i,chrEnzymeBins$binNumb)
  enzymeBinsGrouping[[i]] <- chrEnzymeBins
}
enzymeBinsTab <- do.call(rbind,enzymeBinsGrouping)
binsTab <- enzymeBinsTab[,list(chr=chr[1],start=min(start),end=max(end)),by=binName]
binsTab$name <- sprintf('%s-%s-%s',binsTab$chr,binsTab$start,binsTab$end)
binsGR <- GRanges(seqnames = binsTab$chr,IRanges(start=binsTab$start,end=binsTab$end))
names(binsGR) <- binsTab$name
rm(chrEnzymeBins,enzymeBinsTab)
gc()
gc(reset=TRUE)

quantile(width(binsGR))
0%      25%      50%      75%     100% 
102     2969     3816     4845 30016886 

> sum(as.numeric(width(binsGR)))
[1] 3095676316
> sum(as.numeric(width(reduce(binsGR))))
[1] 3095676316
> 
saveRDS(binsGR,file = sprintf('/home/mark/LiebermanHiC/data/binsGR.rds'))

enzymeCuts <- readRDS('/home/mark/LiebermanHiC/data/enzymeCuts.rds')
binsGR <- readRDS('/home/mark/LiebermanHiC/data/binsGR.rds')

# For each read, we find the distance to the closest enzyme cutting site

ReadsPerChromTab <- list()
for (chr in chrom){
  
  t0 <- proc.time()
  
  names(ReadsPerChrom[[chr]][["R1"]]) <- 1:length(ReadsPerChrom[[chr]][["R1"]])
  names(ReadsPerChrom[[chr]][["R2"]]) <- 1:length(ReadsPerChrom[[chr]][["R2"]])
  
  ReadsPerChromTab[[chr]] <- data.table(R1 = names(ReadsPerChrom[[chr]][["R1"]]),R2 = names(ReadsPerChrom[[chr]][["R2"]]))
  
  nearestEnzcut <- nearest(x = ReadsPerChrom[[chr]][["R1"]],subject=enzymeCuts)
  ReadsPerChromTab[[chr]]$d1 <- distance(x = ReadsPerChrom[[chr]][["R1"]],y = enzymeCuts[nearestEnzcut])
  
  nearestEnzcut <- nearest(x = ReadsPerChrom[[chr]][["R2"]],subject=enzymeCuts)
  ReadsPerChromTab[[chr]]$d2 <- distance(x = ReadsPerChrom[[chr]][["R2"]],y = enzymeCuts[nearestEnzcut])
  
  eT <- (proc.time() - t0)["elapsed"]
  show(sprintf("%s took %3.2f seconds",chr,eT))
}
rm (nearestEnzcut,t0,eT)
gc()
gc(reset=TRUE)


# The max time to complete per chromosome ~ 52 secs
# nrow(ReadsPerChromTab[['chr1']])
# [1] 3041711
# nrow(ReadsPerChromTab[['chr1']][d1 <= 500 & d2 <= 500])
# [1] 3026972
# nrow(ReadsPerChromTab[['chr16']])
# [1] 1031398
# nrow(ReadsPerChromTab[['chr16']][d1 <= 500 & d2 <= 500])
# [1] 1026552

# Filtered Hi-C paired-end reads

Reads <- list()
for(chr in chrom){
  reads <- ReadsPerChromTab[[chr]]
  ind <- which(reads$d1 <= 500 & reads$d2 <= 500)
  Reads[[chr]] <- lapply(ReadsPerChrom[[chr]],function(x){x[ind]})
}
rm (ReadsPerChromTab)
gc()
gc(reset=TRUE)

# Count the number of Hi-C paired reads mapping to bins

registerDoMC(24)
t0 <- proc.time()
 AlignToBins <- foreach(chr = sprintf('%s',chrom)) %dopar% {

  bins <- binsGR[which(seqnames(binsGR)==chr)]
  ovl1 <- findOverlaps(bins,Reads[[chr]][["R1"]])
  ovl1Tab <- data.table(bins = names(bins[queryHits(ovl1)]),
                        reads = names(Reads[[chr]][["R1"]][subjectHits(ovl1)]))
  intersectRanges <- pintersect(bins[queryHits(ovl1)],Reads[[chr]][['R1']][subjectHits(ovl1)])
  ovl1Tab$intersectwidth <- width(intersectRanges)
  ovl1Tab <- ovl1Tab[,list(bins = bins[which.max(intersectwidth)]),by = reads]
  
  ovl2 <- findOverlaps(bins,Reads[[chr]][["R2"]])
  ovl2Tab <- data.table(bins = names(bins[queryHits(ovl2)]),
                        reads = names(Reads[[chr]][["R2"]][subjectHits(ovl2)]))
  intersectRanges <- pintersect(bins[queryHits(ovl2)],Reads[[chr]][['R2']][subjectHits(ovl2)])
  ovl2Tab$intersectwidth <- width(intersectRanges)
  ovl2Tab <- ovl2Tab[,list(bins = bins[which.max(intersectwidth)]),by = reads]
  res <- list("R1" = ovl1Tab, "R2" = ovl2Tab) 
}
eT <- (proc.time() - t0)["elapsed"]
show(sprintf('Time %3.2f seconds',eT))
names(AlignToBins) <- chrom

#offDiag <- 2e6
registerDoMC(24)

t0 <- proc.time()
BinPairCountList <- foreach(chr = sprintf('%s',chrom)) %dopar% {

  setkey(AlignToBins[[chr]][["R1"]],reads)
  PairedBins <- AlignToBins[[chr]][["R1"]][AlignToBins[[chr]][["R2"]]$read]
  PairedBins$bins2 <- AlignToBins[[chr]][["R2"]]$bins
  PairedBins <- PairedBins[,list(reads,binI = bins, binJ = bins2)]
  PairedBins$startI <- as.integer(start(binsGR[PairedBins$binI]))
  PairedBins$startJ <- as.integer(start(binsGR[PairedBins$binJ]))
  PairedBins.Inorder <- PairedBins[startI <= startJ]
  PairedBins.Notordered <- PairedBins[startI > startJ]
  PairedBins.ordered <- rbind(PairedBins.Inorder,PairedBins.Notordered[,
                                                                         list(reads,binI = binJ, binJ = binI, startI = startJ, startJ = startI)])
  PairedBins.ordered <- subset(PairedBins.ordered, startI < startJ)
  countPairedBins <- PairedBins.ordered[,list(counts = .N),by=list(binI,binJ)]
  rm (PairedBins,PairedBins.Inorder,PairedBins.Notordered)
  gc()
  gc(reset=TRUE)
  
  binCenters <- GRanges(seqnames=chr,IRanges(start=start(resize(binsGR[countPairedBins$binI],fix='center',width=1)),
                                             end = start(resize(binsGR[countPairedBins$binJ],fix='center',width=1))))
  countPairedBins$D <- width(binCenters)
  countPairedBins
  #countPairedBins <- subset(countPairedBins,D <= offDiag)
  
  # This section of code that is commented out was to get zero counts paired bins with linear genomic distance <= 2 Mb
  # append it to non-zero count paired bins.
  
  #   regions <- binsGR[which(seqnames(binsGR) == chr)]
  #   regions <- resize(regions,fix="center",width=1)
  #   ovl <- findOverlaps(regions,
  #                            resize(regions,width = offDiag*2,fix = 'center'))
  #   DT <- data.table(binI = queryHits(ovl),binJ = subjectHits(ovl))
  #   DT <- subset(DT,binI < binJ)
  #   binI <- names(regions[DT$binI])
  #   binJ <- names(regions[DT$binJ])
  #   AllPairedBins <- data.table(binI = binI, binJ = binJ)
  #   rm (ovl,DT,regions,binI,binJ)
  #   gc()
  #   gc(reset = TRUE)
  
  #   AllPairedBins$startI <- as.integer(start(binsGR[AllPairedBins$binI]))
  #   AllPairedBins$startJ <- as.integer(start(binsGR[AllPairedBins$binJ]))
  #   AllPairedBins.Inorder <- AllPairedBins[startI <= startJ]
  #   AllPairedBins.Notordered <- AllPairedBins[startI > startJ]
  #   AllPairedBins.ordered <- rbind(AllPairedBins.Inorder,AllPairedBins.Notordered[,
  #                      list(binI = binJ, binJ = binI, startI = startJ, startJ = startI)])
  #   binCenters <- GRanges(seqnames=chr,IRanges(start=start(resize(binsGR[AllPairedBins$binI],fix='center',width=1)),
  #                                              end = start(resize(binsGR[AllPairedBins$binJ],fix='center',width=1))))
  #   AllPairedBins.ordered$D <- width(binCenters)
  #   AllPairedBins <- AllPairedBins.ordered[,list(binI=binI,binJ=binJ,D=D)]
  #   rm(AllPairedBins.ordered,AllPairedBins.Inorder,AllPairedBins.Notordered,binCenters)
  #   gc()
  #   gc(reset=TRUE)
  #   
  #   temp <- countPairedBins[,list(binI=binI,binJ=binJ,D=D)]
  #   setkey(AllPairedBins,key="binI","binJ","D")
  #   zeroPairedBins <- AllPairedBins[!temp]
  #   zeroPairedBins$counts <- 0
  #   allpairedBins <- rbind(countPairedBins,zeroPairedBins)
  #   allpairedBins <- subset(allpairedBins,D <= offDiag)
  #   allpairedBins <- allpairedBins[,list(binI=binI,binJ=binJ,counts=counts,D=D)]
  
  
  # If you have your data in chunks, then you"ll need to run the following lines of code to aggregrate the counts
  #   files <- list.files(path="/Volumes/G-RAID Studio/LiebermanHiC/Contact-Maps/",pattern="Counts_GM12878*",full.names=T,recursive=F)
  #   listDT <- lapply(files,function(x){readRDS(x)[[chrom[chr]]]})
  #   countPairedBins <- do.call(rbind,listDT)
  #   countPairedBins <- countPairedBins[,list(counts=sum(counts),D=D[1]),by=list(binI,binJ)]
  
}
eT <- (proc.time() - t0)["elapsed"]
show(sprintf('Time require to get counts %3.2f seconds',eT))
names(BinPairCountList) <- chrom
saveRDS(BinPairCountList,outfile)
rm (list=ls())
gc()
gc(reset=TRUE)

