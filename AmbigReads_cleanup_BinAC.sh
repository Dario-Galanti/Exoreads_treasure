#!/bin/bash

### Author: Dario Galanti August 2022
### Aim: Remove ambiguous reads from Thlaspi alignments. These are reads mapping to the target (Thlaspi), but also to either the aphid, buchnera or mildew genomes. If not removed, they create false positive SNPs strongly associated to the number of Aphid, Buchnera and Mildew reads
### Input 1): Dir containing Thlaspi alignment bam files
### Input 2): Dir containing Aphid alignment bam files - Contaminant 1
### Input 3): Dir containing Buchnera alignment bam files - Contaminant 2
### Input 4): Dir containing Mildew alignment bam files - Contaminant 3
### Run: bash AmbigReads_cleanup_BinAC.sh

## Define tools and other scripts
samtools=~/miniconda3/envs/bwa/bin/samtools

## Define files
work=/beegfs/work/bbmdg01
ThlaspiDir=${work}/GATK/BWA_alignments
AphidDir=${work}/ConservedRegs_TaGenome/Aphid_alignments
BuchDir=${work}/ConservedRegs_TaGenome/Buch_alignments
MildewDir=${work}/ConservedRegs_TaGenome/Mildew_alignments
outdir=${work}/GATK/BWA_align_noExoReads
#ref_index=/beegfs/work/bbmdg01/Tarvense_genome/v5_thlaspi_published/modified.fasta.fai	# V5 Necessary for opening sam files without headers in samtools
ref_index=/beegfs/work/bbmdg01/Tarvense_genome/v4_thlaspi_final/final.fasta.fai	# V4 Necessary for opening sam files without headers in samtools

## Make directories for individual job scripts, logs and output
mkdir -p $outdir
mkdir -p ${work}/work/Exo_cleanup
mkdir -p ${work}/logs/Exo_cleanup
mkdir -p ${outdir}/Cleanup_reports

for f in ${ThlaspiDir}/TA_*F1_HC0_M1_1.bam;
do
	spl=$(basename $f .bam)
	Aphid_bam=${AphidDir}/${spl}.bam					# Mapping to the Aphid ref genome
	Buch_bam=${BuchDir}/${spl}.bam						# Mapping to the Buchnera ref genome
	Mildew_bam=${MildewDir}/${spl}.bam					# Mapping to the Mildew ref genome
	excludedIDs_tmp=${outdir}/${spl}_excludedIDs_tmp.txt		# Will be created and deleted
	excludedIDs=${outdir}/${spl}_excludedIDs.txt		# Will be created and later compressed
	bamout=${outdir}/$(basename $f)
	report=${outdir}/Cleanup_reports/Cleanup_report.${spl}.txt

## Launch individual sample jobs for fast analysis
	jobName=${work}/work/Exo_cleanup/clean.${spl}_v4.sh
	(
	echo "#PBS -l nodes=1:ppn=4 #Nodes and cores"
	echo "#PBS -l walltime=20:00:00"
	echo "#PBS -l mem=24Gb"							#NB: Allocate more memory that the heaviest bam file!!!
	echo "#PBS -S /bin/bash"
	echo "#PBS -N clean.${spl}"
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/Exo_cleanup/clean.${spl}_v4.out"
	echo "#PBS -e ${work}/logs/Exo_cleanup/clean.${spl}_v4.err"
	echo ""
	
	## Combine reads mapping to Aphid, Buchnera and Mildew genomes
	echo "$samtools view ${Aphid_bam} -F 4 -o  | cut -f1 > $excludedIDs_tmp"
	echo "$samtools view ${Buch_bam} -F 4 -o  | cut -f1 >> $excludedIDs_tmp"
	echo "$samtools view ${Mildew_bam} -F 4 -o  | cut -f1 >> $excludedIDs_tmp"
	echo "cat $excludedIDs_tmp | sort | uniq > $excludedIDs"
	echo "rm $excludedIDs_tmp"
	## Cleanup thlaspi alignments: Remove unmapped, remove reads < 30bp and remove Exogenous reads
	echo "$samtools view -h -F 4 $f | awk '/^@/ || length(\$10) >= 35' | fgrep -v -w -f ${excludedIDs} | $samtools view -b -o ${bamout}"
	echo "$samtools index ${bamout}"
	echo ""
	## Fill report
	echo "mappings=\$($samtools view -F 4 $f | awk 'length(\$10) >= 35' | wc -l)"
	echo "mapped_fragments=\$($samtools view -F 4 $f | awk 'length(\$10) >= 35' | cut -f1 | sort | uniq | wc -l)"
	echo "excluded_fragments=\$(cat $excludedIDs | wc -l)"
	echo "target_mappings=\$($samtools view -F 4 $bamout | wc -l)"
	echo "echo -e ${spl}'\t'\${mappings}'\t'\${mapped_fragments}'\t'\${excluded_fragments}'\t'\${target_mappings} > $report"
	echo "gzip ${excludedIDs}"						# Compress list of excluded IDs
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}
done
exit


