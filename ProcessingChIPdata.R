#####LOADING DATA#####
#Importing home made functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)
library(genomation)
library(Repitools)
library(regioneR)
setwd("C:/Users/ValentinFC/Desktop/Article/IGV")

########################
### DATA IMPORTATION ###
########################

### Import all TSS, Enhancer and GeneBody ###
UniverseTSS<-toGRanges("Genes9_2TSS5Kb_chrosomose.bed")
UniverseGenes<-toGRanges("Genes9_2_chrosomose.bed")
UniverseEnhancers<-toGRanges("Genes9_2Enhancers_chrosomose.bed")

### Import TSS, Enhancer and GeneBody linked to USP21 sensitive genes ###
BackgroundTSS<-toGRanges("USP21TSS_chr.bed") 
BackgroundGenes<-toGRanges("USP21Genes.bed")
BackgroundEnhancers<-toGRanges("USP21Enhancers.bed") 

### Set of peaks to test ###
PeakSet<-toGRanges("H2AK119ub1SpermPeaks.bed")

########################
### PERMUTATION TEST ###
########################

pdf(file = "EnrichmentUSP21ReplicatedU21.pdf",width=10,height=10)

ptTSS<-permTest(A=BackgroundTSS, ntimes=10000,count.once = TRUE, randomize.function=resampleRegions,evaluate.function=numOverlaps, B=PeakSet,allow.overlaps = TRUE, verbose=T,universe = UniverseTSS)
summary(ptTSS)
plot(ptTSS)

ptGenes<-permTest(A=BackgroundGenes, ntimes=10000,count.once = TRUE, randomize.function=resampleRegions,evaluate.function=numOverlaps, B=PeakSet,allow.overlaps = TRUE, verbose=T,universe = UniverseGenes)
summary(ptGenes)
plot(ptGenes)

ptEnhancers<-permTest(A=BackgroundEnhancers, ntimes=10000,count.once = TRUE, randomize.function=resampleRegions,evaluate.function=numOverlaps, B=PeakSet,allow.overlaps = TRUE, verbose=T,universe = UniverseEnhancers)
summary(ptEnhancers)
plot(ptEnhancers)

dev.off()

##########################
### RANDOMISATION TEST ###
##########################

gen<-toGRanges("~/These/Analyses/GFF/XL9_Chromosomes.bed")

pdf(file = "EnrichmentTestFig_CpG_repeats_CO.pdf",width=10,height=10)

PeakSet<-toGRanges("Egg-U21_only.bed")
CpG<-toGRanges("cpgIslandExt.bed")
FigCpG<-permTest(A=PeakSet, B=CpG, ntimes=1000, randomize.function=randomizeRegions, evaluate.function=numOverlaps, count.once=TRUE, genome=gen, allow.overlaps = TRUE,verbose=T)
summary(FigCpG)
plot(FigCpG)

repeats<-toGRanges("repeats.bed")
Figrepeats<-permTest(A=PeakSet, B=repeats, ntimes=1000, randomize.function=randomizeRegions, evaluate.function=numOverlaps, count.once=TRUE, genome=gen, allow.overlaps = TRUE,verbose=T)
summary(Figrepeats)
plot(Figrepeats)

IntergenicRegions<-toGRanges("intergenic.bed") 
FigIntergenic<-permTest(A=IntergenicRegions, B=PeakSet, ntimes=1000, randomize.function=randomizeRegions, evaluate.function=numOverlaps, count.once=FALSE, genome=gen, allow.overlaps = TRUE,verbose=T)
summary(FigIntergenic)
plot(FigIntergenic)

dev.off()

############################
### TF MOTIF DOT HEATMAP ###
############################

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")

library(patchwork) 
library(cowplot)
library(ggtree)
library(theme_cowplot)
library(viridis)
library(ggrepel)
setwd("C:/Users/ValentinFC/Desktop/Article/IGV")

### Process findMotif.pl before
pdf("DotHeatmapHOMERFilter.pdf",width = 100,height=100)

gene_cluster[which(!row.names(gene_cluster)%in%as.numeric(row.names(gene_cluster)[which(gene_cluster$MinlogP.Value.2 > -5 & gene_cluster$MinlogP.Value.3 > -5 & gene_cluster$MinlogP.Value.4 > -5)])),] -> DataLogPVal5

gene_cluster <- read.tsv("../HomerLogPValueInf5.tsv")
ggplot(gene_cluster, aes(x=Cluster, y = MOTIF_NAME, size = log10(1-MinlogP.Value), color=as.numeric(Percentage_motif))) + 
  geom_point() +  
  ylab("GeneMotif") + 
  scale_colour_gradientn(colours = c("red","blue","green"),values = c(0,0.1,1), limits=c(0, 100))


dev.off()

#######################
### PEAK ANNOTATION ###
#######################

library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("~/These/Analyses/GFF/XENLA_9.2_Xenbase.gtf")
setwd("C:/Users/ValentinFC/Desktop/Article/IGV")

promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000)
cntdir <- "./"
pat <- "_0.001.bed"
myfiles <- list.files(path = cntdir,
                      pattern = pat,
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)
pdf(file="ChipSeekerHMMGroups.pdf",width=10,height=10)
for (i in 1:length(myfiles)) {
  peakAnno <- annotatePeak(myfiles[[i]], tssRegion=c(-5000, 5000),TxDb=txdb)
  plotAnnoPie(peakAnno)
}
dev.off()
