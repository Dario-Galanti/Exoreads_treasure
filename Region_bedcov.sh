#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=08:00:00
#PBS -l mem=20Gb 
#PBS -S /bin/bash
#PBS -N Region_bedcov #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/mosdepth/Bed_cov.out
#PBS -e /beegfs/work/bbmdg01/logs/mosdepth/Bed_cov.err

### Author: Dario Galanti Mar 2022
### Aim: Run samtools bedcov on a specific region of the genome for many samples
### Input: 1) region of interest and 2) mapping bam files of individual samples
### Output: unionbed file with positions/bins as rows and samples as columns
### Run: bash Region_mosdepth.sh Chr.start-end binsize
### Run: qsub -q short -F "Scaffold_6.21671000-21687000 20" Region_bedcov.sh


## Define input and output
samtools=~/miniconda3/envs/bwa/bin/samtools
work=/beegfs/work/bbmdg01
wDir=${work}/ConservedRegs_TaGenome/Coverage_Around_peak
bamDir=${work}/GATK/BWA_alignments		# Only files for which we also have methylation data

Chr=$(echo $1 | cut -d. -f1)
start=$(echo $1 | cut -d. -f2 | cut -d- -f1)
end=$(echo $1 | cut -d. -f2 | cut -d- -f2)
binsize=$2						# In theory it can also be 1 but for that mosdepth would have a more proper way to do it

regions=${wDir}/${1}_${2}bp_bins.bed		# Will be created and input to mosdepth
fout=${wDir}/Bedcov_union_${1}_${2}bp_bins.bed


## Create regions bed file
for bin_start in $(seq $start $binsize $end);
do
 bin_end=$(($bin_start + $binsize))
 echo -e $Chr"\t"$bin_start"\t"$bin_end
done > $regions


## Get coverage at every position within specified regions: samtools depth
#$samtools depth -b ${positions} -f bam_files.txt -H -o ${fout}

## Get average coverage of each region: bedtools bedcov for each region (bin) in the positions file (pre-made bins in my case)
##NB: bedtools bedcov reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file.
## Therefore this value has to be divided by the region length

bam_string=$(ls ${bamDir}/TA_*.bam | grep -v TA_AM | tr "\n" " ")
spls_string=$(ls ${bamDir}/TA_*.bam | grep -v TA_AM | rev | cut -d/ -f1 | rev | cut -d. -f1 | tr "\n" " ")
#mkdir -p samples_cov
echo -e chrom"\t"start"\t"end"\t"${spls_string} | tr " " "\t" > ${fout}
$samtools bedcov ${regions} ${bam_string} | awk 'OFS="\t"{len=$3-$2;for(i=4;i<=NF;i++){$i=$i/len};print}' >> ${fout}




