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

### Import all TSS, Enhancer and GeneBody coordinates ###
UniverseTSS<-toGRanges("Genes9_2TSS5Kb_chrosomose.bed")
UniverseGenes<-toGRanges("Genes9_2_chrosomose.bed")
UniverseEnhancers<-toGRanges("Genes9_2Enhancers_chrosomose.bed")

### Import TSS, Enhancer and GeneBody coordinates linked to USP21 sensitive genes ###
BackgroundTSS<-toGRanges("USP21TSS_chr.bed") 
BackgroundGenes<-toGRanges("USP21Genes.bed")
BackgroundEnhancers<-toGRanges("USP21Enhancers.bed") 

### Set of peaks to test ###
PeakSet<-toGRanges("H2AK119ub1SpermPeaks.bed")

########################
### PERMUTATION TEST ###
########################

pdf(file = "EnrichmentUSP21ReplicatedU21.pdf",width=10,height=10)

# Overlap of 10000 resampling of TSS from USP21 sensitive genes over all TSS with H2AK119ub1 sperm peaks
ptTSS<-permTest(A=BackgroundTSS, ntimes=10000,count.once = TRUE, randomize.function=resampleRegions,evaluate.function=numOverlaps, B=PeakSet,allow.overlaps = TRUE, verbose=T,universe = UniverseTSS)
summary(ptTSS)
plot(ptTSS)

# Overlap of 10000 resampling of Gene body from USP21 sensitive genes over all TSS with H2AK119ub1 sperm peaks
ptGenes<-permTest(A=BackgroundGenes, ntimes=10000,count.once = TRUE, randomize.function=resampleRegions,evaluate.function=numOverlaps, B=PeakSet,allow.overlaps = TRUE, verbose=T,universe = UniverseGenes)
summary(ptGenes)
plot(ptGenes)

# Overlap of 10000 resampling of Enhancers associated with USP21 sensitive genes over all TSS with H2AK119ub1 sperm peaks
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
# 1000 randomization of peaks over the genome computing overlap over CpG islands
FigCpG<-permTest(A=PeakSet, B=CpG, ntimes=1000, randomize.function=randomizeRegions, evaluate.function=numOverlaps, count.once=TRUE, genome=gen, allow.overlaps = TRUE,verbose=T)
summary(FigCpG)
plot(FigCpG)

repeats<-toGRanges("repeats.bed")
# 1000 randomization of peaks over the genome computing overlap over repeats elements
Figrepeats<-permTest(A=PeakSet, B=repeats, ntimes=1000, randomize.function=randomizeRegions, evaluate.function=numOverlaps, count.once=TRUE, genome=gen, allow.overlaps = TRUE,verbose=T)
summary(Figrepeats)
plot(Figrepeats)

# 1000 randomization of peaks over the genome computing overlap over intergenic regions (same protocol with Exon, Intron, TSS, Genebody...)
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

### Process findMotif.pl and setup file input
pdf("DotHeatmapHOMERFilter.pdf",width = 100,height=100)

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

############################################
##### R Script TSS Profile P-Value #########
############################################

library(dplyr)
setwd("C:/Users/ValentinFC/Desktop/Article/IGV")
table <- lire("Droso_USP21Sensitive_LadderH2Aub.tsv")
# Select Kb refion around TSS (here 1kb around TSS, 20% bins before and after TSS)
BorneDown <- 81
BorneUp <- 120
# Select sample to compare
ToCompare <- 37
BaseLine <-42

# Select and convert value into numeric
NormalState <- as.numeric(table[BaseLine,c(BorneDown:BorneUp)])
ComparedState <- as.numeric(table[ToCompare,c(BorneDown:BorneUp)])

# Dataframe formating for test processing
my_data <- data.frame(group = rep(c("AllGenes", "MZ"), each = length(NormalState)), values = c(NormalState, ComparedState))

# Paired wilcoxon Mann-Whitney test
wilcox.test(values ~ group, data = my_data, paired = TRUE,alternative="less")

# Extract number of value, median and IQR for each sample
group_by(my_data, group) %>%
  summarise(
    count = n(),
    median = median(values, na.rm = TRUE),
    IQR = IQR(values, na.rm = TRUE)
  )

###############################################################
##### R Script PValues : Output a .probabilities file #########
###############################################################

## PART1 : PValue calculating

## H2A vs NoH2A
test150_110 = seq(1,length(FragmentLen$V1))
test70 = seq(1,length(FragmentLen$V1))
f150_110 <- 0.5*(colSums(FragmentLen[,c(4:7)])[1]/colSums(FragmentLen[,c(4:7)])[4]+colSums(FragmentLen[,c(4:7)])[2]/colSums(FragmentLen[,c(4:7)])[4])
f70 <- colSums(FragmentLen[,c(4:7)])[3]/colSums(FragmentLen[,c(4:7)])[4] 

for (i in 1:length(FragmentLen$V1)){
  if(i%%50000 == 0){
    print(i)}
  if(FragmentLen[i,7]>0){
    test150_110[i] <- prop.test((FragmentLen[i,4]+FragmentLen[i,5]),FragmentLen[i,7],p=f150_110,alternative = "greater")$p.value
    test70[i] <- prop.test(FragmentLen[i,6],FragmentLen[i,7],p=f70,alternative = "greater")$p.value
  } else {
    test150_110[i] <- 1
    test70[i] <- 1
  }
}
pval = cbind(test150_110, test70)
write.table(cbind(FragmentLen,((FragmentLen$V4+FragmentLen$V5)/FragmentLen$V7),(FragmentLen$V6/FragmentLen$V7),pval,1),file=paste0(FragmentLenFile,"_H2A_vs_NonH2A.probabilities",sep=""),sep = '\t',quote = F,col.names = F,row.names = F)

## Nucl vs seminucl
test150 = seq(1,length(FragmentLen$V1))
test110_70 = seq(1,length(FragmentLen$V1))
f150 <- colSums(FragmentLen[,c(4:7)])[1]/colSums(FragmentLen[,c(4:7)])[4] 
f110_70 <- 0.5*(colSums(FragmentLen[,c(4:7)])[2]/colSums(FragmentLen[,c(4:7)])[4]+colSums(FragmentLen[,c(4:7)])[3]/colSums(FragmentLen[,c(4:7)])[4])

for (i in 1:length(FragmentLen$V1)){
  if(i%%50000 == 0){ print(i)}
  if(FragmentLen[i,7]>0){
    test110_70[i] <- prop.test((FragmentLen[i,5]+FragmentLen[i,6]),FragmentLen[i,7],p=f110_70,alternative = "greater")$p.value
    test150[i] <- prop.test(FragmentLen[i,4],FragmentLen[i,7],p=f150,alternative = "greater")$p.value
  }else{
    test110_70[i] <- 1
    test150[i] <- 1
  }}
pval = cbind(test150, test110_70)
write.table(cbind(FragmentLen,(FragmentLen$V4/FragmentLen$V7),((FragmentLen$V5+FragmentLen$V6)/FragmentLen$V7),pval,1),file=paste0(FragmentLenFile,"_Nucl_vs_SemiNucl.probabilities",sep=""),sep = '\t',quote = F,col.names = F,row.names = F)


#############################################################

##################################
##### R Script HMM Model #########
##################################

##PART2 : HMM##
library("RHmm")
HMMFile<-"../AlreadyDone/Input.Fragm150_110_70_sum_nodup_above5_slid50.full.bedgraph.probabilities.discrete"
HMM <- read.table(HMMFile, stringsAsFactors = F, sep = "\t", header = F)
colnames(HMM)<-c("Chr","Begin","End","Count150","Count110","Count70","TotalCount","%150","%110","%70","p150","p110","p70","obs150","obs110","obs70")
range<-4000
test_obs <- HMM
test_obs[is.na(test_obs)] <- 0
hmmfit<- HMMFit(test_obs, nStates= 2, dis= 'DISCRETE')
lim<-seq(1,len(test_obs),by = range)
lim[len(lim)]<-len(test_obs)+1
states<-c()
for(i in 1:(len(lim)-1)){
  vit<-viterbi(hmmfit, test_obs[lim[i]:(lim[i+1]-1)])
  states<-c(states,vit$states)
}
write.table(states,file="Input.Fragm150_above5_hmm_state_50bpw.full",sep = '\t',quote = F,col.names = F,row.names = F)

