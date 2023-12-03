#!/bin/bash

### Aim: WGS accurate read alignment with BWA (http://bio-bwa.sourceforge.net/bwa.shtml) and Mark Duplicates with gatk MarkDuplicatesSpark
### Author: Dario Galanti Nov 2022
### Run: bash BWA_multi_accurate_align_BinAC.sh
### Dependencies: bwa and gatk4 v4.1.8.1

###: To increase mapping accuracy we use a seed length of 23 and remove reads with MAPQ < 20. We also remove reads not mapping to reduce space.

## NB: This script is made for samples sequenced on a single lane. For samples sequenced on multiple lanes the "Extract read group" should be double checked

### PRE-STEPS:
### 1) Prepare conda environment
## conda create -n bwa
## conda install -c bioconda bwa
## conda install -c bioconda gatk4 (if issues with sometools check https://github.com/bioconda/bioconda-recipes/issues/12100 and use user's "Leipzig" solution)
### 2) Index reference genome (-p prefix)
## ~/miniconda3/envs/bwa/bin/bwa index -p Ta_genome /beegfs/work/bbmdg01/Tarvense_genome/v4_thlaspi_final/final.fasta
### 3) Check number of reads per sample
## mkdir fastqc
## echo -e Sample"\t"Mate"\t"Num_reads > Total_reads.txt
## for f in ../WGS_trimmed/trim.TA*;do spl=$(echo $f | rev | cut -c12-34 | rev); mate=$(echo $f | rev | cut -c10); reads=$(zcat $f | grep -c @);echo -e $spl"\t"$mate"\t"$reads >> Total_reads.txt;done

## Define home dir (work), tools, input and output
work=/beegfs/work/bbmdg01
bwa=~/miniconda3/envs/bwa/bin/bwa
gatk=~/miniconda3/envs/bwa/bin/gatk
samtools=~/miniconda3/envs/bwa/bin/samtools					# Only necessary if willing to remove unmapped reads (line 37)
genome=${work}/genomes/Mizus_persicae/Mizus_persicae_genome	# Index genome (see PRE-STEPS).
inDir=${work}/unmapped_fq_v5
outDir=${work}/Aphid_mapping/BWA_Mpersicae_unmapped

## Make directories for individual job scripts, logs and output
mkdir -p ${outDir}
mkdir -p ${work}/work/bwa
mkdir -p ${work}/logs/bwa
mkdir -p ${work}/BWA_align_metrics	# Dir for alignment metrics
mkdir -p ${work}/tmp				# Temporary directory for MarkDuplicatesSpark
#mkdir -p ${work}/java_tmp			# Temporary directory for SortSam

## SUBMIT MULTIPLE ALIGNMENTS
for file in ${inDir}/unmapped_TA_DE*_1.fastq.gz;
do
	sample=$(basename $file _1.fastq.gz | cut -d"_" -f2-)
	jobName=${work}/work/bwa/bwa.${sample}_Mpersicae.sh
	(
	echo "#PBS -l nodes=1:ppn=4 #Nodes and cores"
	echo "#PBS -l walltime=10:00:00"
	echo "#PBS -l mem=8Gb"
	echo "#PBS -S /bin/bash"
	echo "#PBS -N bwa.${sample}"
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/bwa/bwa.${sample}_Mpersicae.out"
	echo "#PBS -e ${work}/logs/bwa/bwa.${sample}_Mpersicae.err"
	echo ""
	
	## Define input and output files
	echo "fin_1=${inDir}/$(basename $file)"
	echo "fin_2=${inDir}/$(basename $file 1.fastq.gz)2.fastq.gz"
	echo "samfile=${outDir}/${sample}.sam"
	echo "bamfile=${outDir}/${sample}.bam"
	echo "metrics=${work}/BWA_align_metrics/dups_metrics_${sample}_Mpersicae.txt"
	echo ""
	## Extract read group ID (this is necessary for GATK!!! But NB: It only works for my dataset where each sample was sequenced on a single lane!!!).
	echo "header=\$(zcat \${fin_1} | head -n 1)"
	echo "id=\$(echo \${header} | cut -f 1-4 -d: | sed 's/@//' | sed 's/:/_/g')"
	echo "sm=\$(echo \${header} | grep -Eo \"[ATGCN]+\$\")"
	echo "RG=\"@RG\tID:\${id}\tSM:${sample}\tLB:\${id}_\${sm}\tPL:ILLUMINA\""
	echo ""
	
	## Run bwa alignment (-R attaches read group, necessary for MarkDuplicatesSpark) removing unmapped reads (samtools -F 4)
	#echo "${bwa} mem -t 8 -R $(echo \${RG}) -c 1 ${genome} \${fin_1} \${fin_2} | $samtools view -h -F 4 -o \${samfile}"	# normal mapping
	## Accurate mapping. We increase seed length (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5643015/) and remove low MAPQ reads
	echo "${bwa} mem -t 8 -R $(echo \${RG}) -c 1 -k 25 ${genome} \${fin_1} \${fin_2} | $samtools view -h -F 4 -q 20 -o \${samfile}"
	
	echo "echo -e \"${sample} bwa alignment finished\n\""
	## MarkDuplicatesSpark (will output a coordinate-sorted bamfile)
	echo "${gatk} MarkDuplicatesSpark -I \${samfile} -O \${bamfile} -M \${metrics} --allow-multiple-sort-orders-in-input --conf 'spark.executor.cores=8' --conf 'spark.local.dir=${work}/tmp' --QUIET"
	echo "rm \${samfile} \${bamfile}.sbi"
	echo "[ -f \${bamfile} ] && echo sample ${sample} bamfile is mapped"
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}

done
exit

## Check results
#samtools=~/miniconda3/envs/bwa/bin/samtools
#for f in BWA_Mpersicae_unmapped/TA_*.bam;do n=$($samtools view $f | cut -f1 | sort | uniq | wc -l);echo -e $(basename $f .bam)"\t"${n} >> Mpersicae_reads.txt;done 

