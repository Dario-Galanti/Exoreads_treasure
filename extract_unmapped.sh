#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=10:00:00
#PBS -l mem=10Gb
#PBS -S /bin/bash
#PBS -N Extract_unmapped #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/Extract_unmapped_v5.out
#PBS -e /beegfs/work/bbmdg01/logs/Extract_unmapped_v5.err

## Author: Dario Galanti
## Aim: Extract unmapped reads from mapping bam files and recover them from the original fastq files
## Input: 1) Unmapped fastq files and 2) Alignment bam files. For mapping: https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/1_BWA_multi_align_BinAC.sh
## Run: qsub -q short extract_unmapped.sh

## NB: For fast processing parallelize the code for each sample. Eg. https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/1_BWA_multi_align_BinAC.sh

samtools=~/miniconda3/envs/bwa/bin/samtools
seqtk=~/miniconda3/envs/bwa/bin/seqtk

work=/beegfs/work/bbmdg01
fqDir=${work}/WGS_trimmed
bwaDir=${work}/GATK/BWA_alignments
unmappedDir=${work}/unmapped_fq_v5		# Output directory, will be created
report=${unmappedDir}/unmapped_report.txt

mkdir -p ${unmappedDir}


## 1) Extract unmapped paired end reads (both reads of the fragment unmapped).
# -b -> bam output
# -f n -> include reads flagged by n. SAM flag 12 --> Both mates unmapped
# -F n -> exclude reads flagged by n.
# SAM flag 256 --> NB: this usually flags secondary alignments (multimappings), but many tools (eg. Picard MarkDuplicates) use it instead of 2048 to flag supplementary alignments (chimeric/split reads) which are extra locations where pieces of reads map.
# My data have no multimapping reads (bwa mem -c 1) and, after running MarkDuplicatesSpark, splitreads flags are converted to 256. 
# sam flags info https://broadinstitute.github.io/picard/explain-flags.html
# See my comment on https://www.biostars.org/p/228787/#9519609

for f in ${bwaDir}/TA_DE_01_02*_M1_1.bam
do
 spl=$(basename $f .bam)
 
 ## 1)Extract unmapped paired end reads (both reads of the fragment unmapped) from alignment bam files
 unmapped=${unmappedDir}/unmapped_$(basename $f)
 $samtools view -b -f 12 -F 256 $f > $unmapped
 
 ## 2) Extract unmapped reads from fastq
 read_IDs=${unmappedDir}/unmappedIDs_${spl}.lst
 fin1=${fqDir}/trim.${spl}_1.fastq.gz
 fin2=${fqDir}/trim.${spl}_2.fastq.gz
 fout1=${unmappedDir}/unmapped_${spl}_1.fastq.gz
 fout2=${unmappedDir}/unmapped_${spl}_2.fastq.gz
 $samtools view $unmapped | cut -f1 | uniq > $read_IDs
 $seqtk subseq $fin1 $read_IDs | gzip > $fout1
 $seqtk subseq $fin2 $read_IDs | gzip > $fout2
done


## Make report
echo -e Fragments"\t"sample > $report
for f in ${unmappedDir}/unmappedIDs_*.lst
do
	wc -l $f >> $report
done

## Optional cleanup
rm ${unmappedDir}/unmapped_*.bam
#rm ${unmappedDir}/unmappedIDs_*.lst

