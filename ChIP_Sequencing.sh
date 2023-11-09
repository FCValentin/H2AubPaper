#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e  ./log/
#$ -o ./log/

bwa index -a bwtsw XL9-2_dm6-40_601DNA.fa -p XL-DM-601DNA_Hybrid
conda activate ChIP_SeqPipeline
mkdir rawfastq
mkdir fastq
mkdir FingerPrint
mkdir FragmentSize
mkdir bwa
mkdir cutadapt
mkdir macs2
mkdir log


## Put your FASTQ in rawfastq folder
## RENAME your sample like : <SAMPLE_NAME><HISTONEMARK_ChIP or INPUT_Input><SEQUENCINGLANEparameter><fastq.gz>
#Example : Sperm1_H2AK119ub1_ChIP_L001_R1_001.fastq.gz and Sperm1_INPUT_Input_L001_R1_001.fastq.gz

############################
### FOR PAIRED-END DATA ####
############################

#Merge Lanes
for f in `ls -1 rawfastq/*_L001_R1_001.fastq.gz | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'`
	do 
	cat rawfastq/${f}_L001_R1_001.fastq.gz rawfastq/${f}_L002_R1_001.fastq.gz > fastq/${f}_1.fq.gz
	cat rawfastq/${f}_L001_R2_001.fastq.gz rawfastq/${f}_L002_R2_001.fastq.gz > fastq/${f}_2.fq.gz
done

#low quality reads trimming
for f in `ls -1 fastq/*_1.fq.gz | cut -d "/" -f 2 | sed 's/_1.fq.gz//'`
	do cutadapt -q 20 -m 20 --pair-filter=any -o cutadapt/cutadapt_${f}_1.fastq.gz -p cutadapt/cutadapt_${f}_2.fastq.gz fastq/${f}_1.fq.gz fastq/${f}_2.fq.gz
	fastqc cutadapt/cutadapt_${f}_1.fastq.gz cutadapt/cutadapt_${f}_2.fastq.gz
done

#Alignment PE
for f in `ls -1 cutadapt/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do bwa mem -M -t 8 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/XL-DM-601DNA_Hybrid cutadapt/${f}_1.fastq.gz cutadapt/${f}_2.fastq.gz | samtools view -hbS | samtools sort > bwa/${f}.sort.bam
done

############################
### FOR SINGLE-END DATA ####
############################

#Merge Lanes
for f in `ls -1 rawfastq/*_L001_R1_001.fastq.gz | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'`
	do 
	cat rawfastq/${f}_L001_R1_001.fastq.gz rawfastq/${f}_L002_R1_001.fastq.gz > fastq/${f}_1.fq.gz
done

#low quality reads trimming
for f in `ls -1 fastq/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do cutadapt -q 20 -m 20 -o cutadapt/cutadapt_${f}_1.fastq.gz fastq/${f}_1.fastq.gz 
	fastqc cutadapt/cutadapt_${f}_1.fastq.gz
done

for f in `ls -1 cutadapt/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do bwa mem -M -t 8 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/XL-DM-601DNA_Hybrid cutadapt/${f}_1.fastq.gz | samtools view -hbS | samtools sort > bwa/${f}.sort.bam
done

###################################
### READ EXTRACTION PER GENOME ####
###################################

# Index BAM file, and extract subgenomes chromosomes (DNA ladder, Drosophila melanogaster, Xenopus laevis)
for f in `ls -1 bwa/*.bam  | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools index bwa/${f}.sort.bam 
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam chr1L chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S | samtools view -bS | samtools sort > bwa/${f}.XL92.sort.bam 
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam chr_2L_fly chr_2R_fly chr_3L_fly chr_3R_fly chr_4_fly chr_X_fly chr_Y_fly | samtools view -bS | samtools sort > bwa/${f}.DM6.sort.bam
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam 601_DNA | samtools view -bS | samtools sort > bwa/${f}.DNA.sort.bam
	   samtools index bwa/${f}.XL92.sort.bam 
	   samtools index bwa/${f}.DM6.sort.bam
done

################
### DATA QC ####
################

#Sample quality
for f in `ls -1 bwa/*sort.bam | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools flagstat bwa/${f}.bam > log/${f}.txt
done

#Replicates quality
multiBamSummary bins -p 4 --verbose --bamfiles <ALL_BAMFILES> -o Correlation/ChipCorrelationChIP_Input.npz
plotCorrelation --removeOutliers --corData Correlation/ChipCorrelationChIP_Input.npz -c spearman -p heatmap --colorMap bwr --plotNumbers --labels <ALL_SamplesLABEL> -o Correlation/HeatmapCorChIP_Inputspearman.svg

###############################
### PCR duplicates removal ####
###############################

# mark duplicates with Picard && index result bam for DM6 and XL92 genome, on each <HISTONEMARK> and <INPUT>
for f in `ls -1 bwa/*_INPUT_Input.sort.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.sort.bam//'`
	do 
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H3K4me3_ChIP.XL92.sort.bam O=bwa/${f}_H3K4me3_ChIP.XL92.duplicates.bam M=bwa/${f}_H3K4me3_ChIP.XL92.dup_metrics.txt
	   samtools index bwa/${f}_H3K4me3_chIP.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H2Aub_ChIP.XL92.sort.bam O=bwa/${f}_H2Aub_ChIP.XL92.duplicates.bam M=bwa/${f}_H2Aub_ChIP.XL92.dup_metrics.txt
	   samtools index bwa/${f}_H2Aub_chIP.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_INPUT_Input.XL92.sort.bam O=bwa/${f}_INPUT_Input.XL92.duplicates.bam M=bwa/${f}_INPUT_Input.XL92.dup_metrics.txt
	   samtools index bwa/${f}_INPUT_Input.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H3K4me3_ChIP.DM6.sort.bam O=bwa/${f}_H3K4me3_ChIP.DM6.duplicates.bam M=bwa/${f}_H3K4me3_ChIP.DM6.dup_metrics.txt
	   samtools index bwa/${f}_H3K4me3_chIP.DM6.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H2Aub_ChIP.DM6.sort.bam O=bwa/${f}_H2Aub_ChIP.DM6.duplicates.bam M=bwa/${f}_H2Aub_ChIP.DM6.dup_metrics.txt
	   samtools index bwa/${f}_H2Aub_chIP.DM6.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_INPUT_Input.DM6.sort.bam O=bwa/${f}_INPUT_Input.DM6.duplicates.bam M=bwa/${f}_INPUT_Input.DM6.dup_metrics.txt
	   samtools index bwa/${f}_INPUT_Input.DM6.duplicates.bam
done	   

########################
### Histone mark QC ####
########################

# Fragmentsize only on PAIRED-END DATA
for f in `ls -1 bwa/*_INPUT_Input.sort.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.sort.bam//'`
	do 
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.XL92.sort.bam bwa/${f}_H3K4me3_ChIP.XL92.sort.bam bwa/${f}_INPUT_Input.XL92.sort.bam -plot FingerPrint/${f}_fingerprint.XL92.svg 
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.DM6.sort.bam bwa/${f}_H3K4me3_ChIP.DM6.sort.bam bwa/${f}_INPUT_Input.DM6.sort.bam -plot FingerPrint/${f}_fingerprint.DM6.svg
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.sort.bam bwa/${f}_H3K4me3_ChIP.sort.bam bwa/${f}_INPUT_Input.sort.bam -plot FingerPrint/${f}_fingerprint.svg
	   
	   ### ONLY ON PAIRED-END DATA
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.XL92.sort.bam bwa/${f}_H3K4me3_ChIP.XL92.sort.bam bwa/${f}_INPUT_Input.XL92.sort.bam -o FragmentSize/${f}_FragmentSize.XL92.svg
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.DM6.sort.bam bwa/${f}_H3K4me3_ChIP.DM6.sort.bam bwa/${f}_INPUT_Input.DM6.sort.bam -o FragmentSize/${f}_FragmentSize.DM6.svg
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.sort.bam bwa/${f}_H3K4me3_ChIP.sort.bam bwa/${f}_INPUT_Input.sort.bam -o FragmentSize/${f}_FragmentSize.svg
done

# filter out duplicates
for f in `ls -1 bwa/*.duplicates.bam | cut -d "/" -f 2 | sed 's/.duplicates.bam//'`
	do samtools view -hb -F 0x400 -q 20 bwa/${f}.duplicates.bam | samtools sort -@ 6 > bwa/${f}.sort.filter.bam
	   samtools index bwa/${f}.sort.filter.bam 
done

#############################################
### DATA QC after PCR duplicates removal ####
#############################################


for f in `ls -1 bwa/*.sort.filter.bam   | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools flagstat bwa/${f}.bam > log/${f}.txt
done

#############################################
### PEAKS calling and DATA visualisation ####
#############################################

for f in `ls -1 bwa/*_INPUT_Input.XL92.sort.filter.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.XL92.sort.filter.bam//'`
	do 
	mkdir macs2/${f}
	macs2 callpeak -t bwa/${f}_H3K4me3_ChIP.DM6.sort.filter.bam -c bwa/${f}_INPUT_Input.DM6.sort.filter.bam -f BAMPE --gsize=1.2e8 -n macs2/${f}/${f}.DM6 -B -q 0.01 --nomodel --extsize 147
	macs2 callpeak -t bwa/${f}_H3K4me3_ChIP.XL92.sort.filter.bam -c bwa/${f}_INPUT_Input.XL92.sort.filter.bam -f BAMPE --gsize=2.6e9 -n macs2/${f}/${f}.XL92 -B -q 0.01 --nomodel --extsize 147
	
	bamToBed -i bwa/${f}_H2Aub_ChIP.DM6.sort.filter.bam > macs2/${f}_H2Aub_chIP.DM6.sort.filter.bed
	bamToBed -i bwa/${f}_INPUT_Input.DM6.sort.filter.bam > macs2/${f}_INPUT_Input.DM6.sort.filter.bed
	bamToBed -i bwa/${f}_H2Aub_ChIP.XL92.sort.filter.bam > macs2/${f}_H2Aub_chIP.XL92.sort.filter.bed
	bamToBed -i bwa/${f}_INPUT_Input.XL92.sort.filter.bam > macs2/${f}_INPUT_Input.XL92.sort.filter.bed

	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H2Aub_chIP.XL92.sort.filter.bam -b2 bwa/${f}_INPUT_Input.XL92.sort.filter.bam -o bigwig/${f}_H2Aub_log2ratio.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H3K4me3_ChIP.XL92.sort.filter.bam -b2 bwa/${f}_INPUT_Input.XL92.sort.filter.bam -o bigwig/${f}_H3K4me3_log2ratio.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H2Aub_chIP.DM6.sort.filter.bam -b2 bwa/${f}_INPUT_Input.DM6.sort.filter.bam -o bigwig/${f}_H2Aub_log2ratio_DM6.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H3K4me3_ChIP.DM6.sort.filter.bam -b2 bwa/${f}_INPUT_Input.DM6.sort.filter.bam -o bigwig/${f}_H3K4me3_log2ratio_DM6.bw
	
	#activate py2 env for RECOGNICER
	conda activate py2
	cd macs2/${f}
	sh ../../RECOGNICER_Droso.sh ../${f}_H2Aub_chIP.DM6.sort.filter.bed ../${f}_INPUT_Input.DM6.sort.filter.bed 0.001
	sh ../../RECOGNICER.sh ../${f}_H2Aub_chIP.XL92.sort.filter.bed ../${f}_INPUT_Input.XL92.sort.filter.bed 0.001
	cd ../../
	conda activate ChIP_SeqPipeline
done

for f in `ls -1 bwa/*_INPUT_Input.sort.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.sort.bam//'`
	do 
	bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_INPUT_Input.XL92.sort.bam -o FragmentSize/${f}_FragmentSize_INPUT.XL92.svg
done


#Fig 2A example
computeMatrix reference-point -p 4 --referencePoint center --verbose -S cutadapt_oo-WT_H2Aub_log2ratio.bw cutadapt_oo-U21_H2Aub_log2ratio.bw cutadapt_Egg-oo-WT_H2Aub_log2ratio.bw cutadapt_Egg-oo-U21_H2Aub_log2ratio.bw cutadapt_oo-WT_H3K4me3_log2ratio.bw cutadapt_oo-U21_H3K4me3_log2ratio.bw Merged_cutadapt_Egg-oo-WT_H3K4me3_log2ratio.bw Merged_cutadapt_Egg-oo-U21_H3K4me3_log2ratio.bw -R Egg-WT_only.bed Egg-commonU21-WT_only.bed Egg-U21_only.bed Egg_WT_K4Only.bed Egg_commonK4_WT-U21.bed Egg_U21_K4Only.bed --binSize 50 --missingDataAsZero --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -out DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz --outFileNameMatrix DeeptoolsMatrixDynamicsReplication_U21_WT.tab --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT.bed
plotHeatmap --samplesLabel WTH2Aub USP21H2Aub RepWTH2Aub RepUSP21H2Aub WTK4 USP21K4 RepWTK4 RepUSP21K4 --regionsLabel WT common U21 WT Common U21 -m DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz -out DeeptoolsMatrixDynamicsReplication_U21_WT.svg --colorMap bwr --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT.bed

#######################################
### HOMER MOTIF DISCOVERY ON PEAKS ####
#######################################

bedtools getfasta -fi XL9_2.fa -bed Cluster3_USP21.sort.bed -tab > MEME_sequences_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_ClusterH2Aub.bed | sed 's/\t/\n/g' > Sequences_H2Aub.fasta

bedtools getfasta -fi XL9_2.fa -bed Background_Cluster3.bed -tab > MEME_sequences_Background_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_Background_ClusterH2Aub.bed | sed 's/\t/\n/g' > MEME_sequences_Background_ClusterH2Aub.fasta

#H2Aub Cluster
findMotifs.pl Sequences_H2Aub.fasta fasta Cluster_H2Aub -fasta MEME_sequences_Background_ClusterH2Aub.fasta > H2AubCluster.txt
