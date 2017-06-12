#################################################################################
#################################################################################
# setwd("/Users/guywachsman/Guy/EMS1/pipeline/M2_194-proof/")
# mntn=read.delim("/Users/guywachsman/Guy/EMS1/pipeline/M2_194-proof/EMS.allSNPs.txt", header=T)
# cnds=read.delim("/Users/guywachsman/Guy/EMS1/pipeline/M2_194-proof/EMS.candidates.txt", header =T, sep = "\t")
# setwd("/Users/guywachsman/Guy/EMS1/pipeline/SIMPLE/Simple")


setwd("./output")
library("ggplot2")
library("ggrepel")
library("reshape2")

# library("grid")
# library("stringr")
# library("lattice")
# library("plyr")
# library(scales)
# library(KernSmooth)

#JEN inputting the line information and setting filenames
args <- commandArgs(trailingOnly=TRUE)
line <- args[1]
plotfilename <- paste(line, ".allSNPs.txt", sep="")
candidatefilename <- paste(line, ".candidates.txt", sep="")


mntn=read.delim(plotfilename, header=T)#In my previous version of R, but NOT on the bioross server, R adds an extra NA column
cnds=read.delim(candidatefilename, header =T, sep = "\t")

#mntn1=droplevels(mntn[((mntn$chr!="Mt")&(mntn$chr!="Pt")),])
mntn1= mntn[!(mntn$chr=="Mt" | mntn$chr=="Pt"),]


#mntn1=droplevels(mntn[((mntn$ref=="C"&mntn$alt=="T")|(mntn$ref=="G"&mntn$alt=="A"))&(mntn$chr!="Mt")&(mntn$chr!="Pt"),])

# tbl1[tbl1[,1]=="AT3G28550",]
# tail(mntn, n=10)
# dim(mntn1)


wt=mntn1$wt.ref/(mntn1$wt.ref+mntn1$wt.alt)
mut=mntn1$mut.ref/(mntn1$mut.ref+mntn1$mut.alt)
ratio=wt-mut

tbl1=data.frame(At_num=mntn1$At_num, gene=mntn1$gene, chr=mntn1$chr, pos=mntn1$pos, mut.ref=mntn1$mut.ref, mut.alt=mntn1$mut.alt, wt.ref=mntn1$wt.ref, wt.alt=mntn1$wt.alt, mut.ratio=mut, wt.ratio=wt, ratio)
tbl1=tbl1[complete.cases(tbl1),]

#tbl1.cands=tbl1[(tbl1[,3] %in% cnds[,1]) & (tbl1[,4] %in% cnds[,2]),]
breaks=seq(0, max(tbl1$pos), round(max(tbl1$pos)/3, digits=-7))

#########################################################################
#separated chromosomes original data; after filtering (tbl2)-LOESS fitted
#removing ratios below 0.1
#########################################################################
tbl2=tbl1[(tbl1[,11]>0.1),]
t2_s=split(tbl2, tbl2$chr)
lll=lapply(t2_s, function(x) {loess(x$ratio~x$pos, degree=2, span=0.3)})
mmm=lapply(lll, '[[', 'fitted')
fitted=Reduce(c, mmm)
tbl3=data.frame(tbl2, fitted)
tbl3.cands=tbl3[(tbl3[,3] %in% cnds[,1]) & (tbl3[,4] %in% cnds[,2]),]

#making the x-axes labels less messy
fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     l <- gsub("0e\\+00","0",l)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e\\+","e",l)
     l <- gsub("e", "%*%10^", l)
     l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
     l <- gsub("\\.0", "", l)
     #l <- gsub("\\'[\\.0]'", "", l)
     # return this as an expression
     parse(text=l)
}

#getting the loess fitted plot
x11()
ggplot(tbl3, aes(pos, fitted)) + geom_point(aes(color=chr),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands, aes(x=pos, y=fitted), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=fitted, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="ratio")

#JEN changed file name
Rplot_loess1_file <- paste(line, ".Rplot.loess.1.pdf", sep="")
ggsave(Rplot_loess1_file)

#########################################################################
#separated chromosomes original data; after filtering (tbl2)-LOESS fitted
#removing ratios below 0.3
#########################################################################
tbl2=tbl1[(tbl1[,11]>0.3),]
t2_s=split(tbl2, tbl2$chr)
lll=lapply(t2_s, function(x) {loess(x$ratio~x$pos, degree=2, span=0.3)})
mmm=lapply(lll, '[[', 'fitted')
fitted=Reduce(c, mmm)
tbl3=data.frame(tbl2, fitted)
tbl3.cands=tbl3[(tbl3[,3] %in% cnds[,1]) & (tbl3[,4] %in% cnds[,2]),]

#making the x-axes labels less messy
fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     l <- gsub("0e\\+00","0",l)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e\\+","e",l)
     l <- gsub("e", "%*%10^", l)
     l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
     l <- gsub("\\.0", "", l)
     #l <- gsub("\\'[\\.0]'", "", l)
     # return this as an expression
     parse(text=l)
}

#getting the loess fitted plot
x11()
ggplot(tbl3, aes(pos, fitted)) + geom_point(aes(color=chr),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands, aes(x=pos, y=fitted), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=fitted, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="ratio")

# JEN changed file name
Rplot_loess3_file <- paste(line, ".Rplot.loess.3.pdf", sep="")
ggsave(Rplot_loess3_file)

#and getting the allele frequency
tbl3.m=melt(tbl3, id.vars=c('At_num', 'gene', 'chr', 'pos', 'mut.ref', 'mut.alt', 'wt.ref', 'wt.alt', 'ratio', 'fitted'))

tbl3.cands.m=melt(tbl3.cands, id.vars=c('At_num', 'gene', 'chr', 'pos', 'mut.ref', 'mut.alt', 'wt.ref', 'wt.alt', 'ratio', 'fitted'))


x11()
ggplot(tbl3.m, aes(pos, value)) + geom_point(aes(color=variable),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands.m, aes(x=pos, y=value), shape=5)+geom_text_repel(data=tbl3.cands.m, aes(x=pos, y=value, label=gene), size=3)+theme()+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="allele frequency")

#JEN changed file name
Rplot_allele_file <- paste(line, ".Rplot_allele.pdf", sep="")
ggsave(Rplot_allele_file)
#######################################
#######################################

#JEN - comment - best to remove this stuff below, don't you think.

#separated chromosomes original data; after filtering (tbl2)
# ggplot(tbl3, aes(x=pos, y=ratio))+geom_point(data=tbl3, aes(x=pos, y=ratio, color=as.factor(tbl3[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl3.cands, aes(x=pos, y=ratio), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=ratio, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)

#x11()
#separated chromosomes original data
#ggplot(tbl3, aes(x=pos, y=ratio))+geom_point(data=tbl3, aes(x=pos, y=ratio, color=as.factor(tbl3[,3])),size=0.3)+facet_wrap(~chr, scales='free_x')+geom_point(data=tbl3.cands, aes(x=pos, y=ratio), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=ratio, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks)

















