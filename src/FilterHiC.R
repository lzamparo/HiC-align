library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(doMC)
registerDoMC(24)
args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)
n  <- as.integer(args[1])

  files <- sprintf('/home/mark/LiebermanHiC/BAM/SRR16585%d_%d_aln.sorted.bam',n,1:2)
  file <- sprintf('/home/mark/LiebermanHiC/Reads/SRR16585%d.rds',n)
  chrom <- sprintf('chr%s',c(1:22,"X","Y"))
  stime <- system.time({
          res <- foreach(i=1:24,.combine=c) %dopar% {
                  chr <- chrom[i]
              what <- c("rname","qname","strand","pos","mapq")
              which <- GRanges(chr,IRanges(start=1,end=seqlengths(Hsapiens)[chr]))
              param <- ScanBamParam(what = what,which=which,
                flag = scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE),
                tag=c("AS","XS"))
              rdL <- readGAlignments(files[1], use.names=TRUE, param=param)
              rdR <- readGAlignments(files[2], use.names=TRUE, param=param)
              rdL.filter1 <- rdL[which(elementMetadata(rdL)$mapq!=0)]
              rdR.filter1 <- rdR[which(elementMetadata(rdR)$mapq!=0)]
                  rm(rdL,rdR)
              gc()

              gc(reset=TRUE)
              dL <- elementMetadata(rdL.filter1)$AS - elementMetadata(rdL.filter1)$XS
              dR <- elementMetadata(rdR.filter1)$AS - elementMetadata(rdR.filter1)$XS
              rdL.filter2 <- rdL.filter1[dL > quantile(dL,probs=0.1)]
              rdR.filter2 <- rdR.filter1[dR > quantile(dR,probs=0.1)]
              rm(rdL.filter1,rdR.filter1,dL,dR)
              gc()

              gc(reset=TRUE)
              RL <- data.table(id=names(rdL.filter2),chrom1=as.character(seqnames(rdL.filter2)),pos1=start(rdL.filter2),width1=width(rdL.filter2),
                strand1=as.character(strand(rdL.filter2)))
              RR <- data.table(id=names(rdR.filter2),chrom2=as.character(seqnames(rdR.filter2)),pos2=start(rdR.filter2),width2=width(rdR.filter2),
                strand2=as.character(strand(rdR.filter2)))
              rm(rdL.filter2,rdR.filter2)
              gc()

              gc(reset=TRUE)
              setkey(RL,"id")
              setkey(RR,"id")
              R <- merge(RL,RR,by="id")
              G1 <- GRanges(R$chrom1,IRanges(start=R$pos1,width=R$width1),strand=R$strand1)
              seqlengths(G1) <- seqlengths(Hsapiens)[chrom[i]]
                  G2 <- GRanges(R$chrom2,IRanges(start=R$pos2,width=R$width2),strand=R$strand2)
                  seqlengths(G1) <- seqlengths(Hsapiens)[chrom[i]]
              common <- !duplicated(G1) & !duplicated(G2)
              H1 <- G1[common]
              H2 <- G2[common]
              H <- list(H1,H2)
                }
          })[3]/60
  show(sprintf('Time required to get reads: %3.2f mins',stime))
  rm(dL,dR,R,G1,G2,H1,H2,RL,RR,common)
  gc()

  gc(reset=TRUE)
  R1 <- do.call(c,res[2*(1:(length(res)/2))-1])
  R2 <- do.call(c,res[2*(1:(length(res)/2))])
  reads <- GRangesList(R1=R1,R2=R2)
  rm(R1,R2,param,what,which)
  gc()
  gc(reset=TRUE)
  saveRDS(reads,file=file)