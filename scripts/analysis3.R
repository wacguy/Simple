#################################################################################
#################################################################################
# setwd("/Users/guywachsman/Guy/EMS/pipeline/M3_300_7/")
# mntn=read.delim("M3_300_7.plot4.txt", header=T)
# cnds=read.delim("M3_300_7.cands4.txt", header =T, sep = "\t")

setwd("./output")
library("ggplot2")
library("ggrepel")

# library("reshape2")
# library("grid")
# library("stringr")
# library("lattice")
# library("plyr")
# library(scales)
# library(KernSmooth)

mntn=read.delim("plot4.txt", header=T)#In my previous version of R, but NOT on the bioross server, R adds an extra NA column
#mntn=subset(mntn, select=-X)
cnds=read.delim("cands4.txt", header =T, sep = "\t")

mntn1=droplevels(mntn[((mntn$ref=="C"&mntn$alt=="T")|(mntn$ref=="G"&mntn$alt=="A"))&(mntn$chr!="Mt")&(mntn$chr!="Pt"),])

# tbl1[tbl1[,1]=="AT3G28550",]
# head(mntn)

wt=mntn1$wt.ref/(mntn1$wt.ref+mntn1$wt.alt)
mut=mntn1$mut.ref/(mntn1$mut.ref+mntn1$mut.alt)
ratio=wt-mut

tbl1=data.frame(At_num=mntn1$At_num, gene=mntn1$gene, chr=mntn1$chr, pos=mntn1$pos, mut.ref=mntn1$mut.ref, mut.alt=mntn1$mut.alt, wt.ref=mntn1$wt.ref, wt.alt=mntn1$wt.alt, wt.ratio=wt, mut.ratio=mut, ratio)
tbl1=tbl1[complete.cases(tbl1),]

#tbl1.cands=tbl1[(tbl1[,3] %in% cnds[,1]) & (tbl1[,4] %in% cnds[,2]),]
breaks=seq(0, max(tbl1$pos), round(max(tbl1$pos)/3, digits=-7))
tbl2=tbl1[(tbl1[,11]>0.1),]

#########################################################################
#separated chromosomes original data; after filtering (tbl2)-LOESS fitted
#########################################################################
tbl2=tbl1[(tbl1[,11]>0.1),]
t2_s=split(tbl2, tbl2$chr)
lll=lapply(t2_s, function(x) {loess(x$ratio~x$pos, degree=2, span=0.3)})
mmm=lapply(lll, '[[', 'fitted')
fitted=Reduce(c, mmm)
tbl3=data.frame(tbl2, fitted)
tbl3.cands=tbl3[(tbl3[,3] %in% cnds[,1]) & (tbl3[,4] %in% cnds[,2]),]


x11()
ggplot(tbl3, aes(x=pos, y=fitted))+geom_point(data=tbl3, aes(x=pos, y=fitted, color=as.factor(tbl3[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl3.cands, aes(x=pos, y=fitted), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=fitted, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)+labs(x="position", y="ratio")

ggsave("Rplot.pdf")
#######################################
#######################################

#separated chromosomes original data; after filtering (tbl2)
# ggplot(tbl2, aes(x=pos, y=ratio))+geom_point(data=tbl2, aes(x=pos, y=ratio, color=as.factor(tbl2[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl1.cands, aes(x=pos, y=ratio), shape=5)+geom_text_repel(data=tbl1.cands, aes(x=pos, y=ratio, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)

#x11()
#separated chromosomes original data
#ggplot(tbl1, aes(x=pos, y=ratio))+geom_point(data=tbl1, aes(x=pos, y=ratio, color=as.factor(tbl1[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl1.cands, aes(x=pos, y=ratio), shape=5)+geom_text_repel(data=tbl1.cands, aes(x=pos, y=ratio, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)








