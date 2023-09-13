# Exoreads_treasure
This repo contains scripts to extract exogenous reads, not mapping to the target species, from large sequencing datasets.

We used them on a sequencing dataset of more than 200 Thlaspi arvense lines, published previously [here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010452)

[extract_unmapped.sh](https://github.com/Dario-Galanti/Exoreads_treasure/blob/main/extract_unmapped.sh)<br/>
Extract exogenous (unmapped) reads from mapping bam files and recovers them from the original fastq files. For mapping refer to my previous [script](https://github.com/Dario-Galanti/BinAC_varcalling/blob/main/1_BWA_multi_align_BinAC.sh). For large datasets parallelization can be implemented submitting each sample as a separate job.

[BWA_multi_accurate_align_BinAC.sh](https://github.com/Dario-Galanti/Exoreads_treasure/blob/main/BWA_multi_accurate_align_BinAC.sh)<br/>
Perform high confidence alignment of multiple samples in parallel.
We used it to map non-target reads (not mapping to the T. arvense genome) to the aphid genome Mizus persicae and its symbiont Buchnera aphidicola, to quantify aphid infestation of our T. arvense collection.

[AmbigReads_cleanup_BinAC.sh](https://github.com/Dario-Galanti/Exoreads_treasure/blob/main/AmbigReads_cleanup_BinAC.sh)<br/>
Remove ambiguous reads from the target species alignments (T. arvense in our case). These are reads mapping to the target (T. arvense), but also to contaminants (either the aphid, buchnera or mildew genomes). If not removed, these ambiguous reads can create false positive SNPs strongly associated to the number of Aphid, Buchnera and Mildew reads.

