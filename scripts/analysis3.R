#################################################################################
#################################################################################
setwd("./output")
library("reshape2")
library("grid")
library("stringr")
library("ggplot2")
library("ggrepel")
library("lattice")
library(KernSmooth)
library("plyr")
library(scales)


mntn=read.delim("plot4.txt", header=T)#In my version of R, but NOT on the bioross server, R adds an extra NA column
#mntn=subset(mntn, select=-X)
cnds=read.delim("cands4.txt", header =T, sep = "\t")

mntn1=droplevels(mntn[((mntn$ref=="C"&mntn$alt=="T")|(mntn$ref=="G"&mntn$alt=="A"))&(mntn$chr!="Mt")&(mntn$chr!="Pt"),])

wt=mntn1$wt.ref/(mntn1$wt.ref+mntn1$wt.alt)
mut=mntn1$mut.ref/(mntn1$mut.ref+mntn1$mut.alt)
ratio=wt-mut
head(mntn1)

tbl1=data.frame(At_num=mntn1$At_num, gene=mntn1$gene, chr=mntn1$chr, pos=mntn1$pos, mut.ref=mntn1$mut.ref, mut.alt=mntn1$mut.alt, wt.ref=mntn1$wt.ref, wt.alt=mntn1$wt.alt, wt.ratio=wt, mut.ratio=mut, ratio)
head(tbl1)

 tbl1.cands=tbl1[(tbl1[,3] %in% cnds[,1]) & (tbl1[,4] %in% cnds[,2]),]
breaks=seq(0, max(tbl1$pos), round(max(tbl1$pos)/3, digits=-7))


x11()
#separated chromosomes original data
ggplot(tbl1, aes(x=pos, y=ratio))+geom_point(data=tbl1, aes(x=pos, y=ratio, color=as.factor(tbl1[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl1.cands, aes(x=pos, y=ratio), shape=5)+geom_text_repel(data=tbl1.cands, aes(x=pos, y=ratio, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)

ggsave("Rplot.pdf")











