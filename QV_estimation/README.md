# QV estimation
---

The **Q**uality **V**alue (**QV**) of a genome assembly is a [Phred-scaled](https://en.wikipedia.org/wiki/Phred_quality_score) estimate of the number of artifact SNP and INDEL errors present in the short-read mappable portion of the reference genome. For further details on the calculation of this parameter and to understand how it can be used to compare reference assembly fidelity, please see [Bickhart et al. 2017](https://www.nature.com/articles/ng.3802). Special thanks to [Sergey Koren](https://scholar.google.com/citations?user=QQBHMyAAAAAJ&hl=en) for developing this metric!

## Installation

Thankfully, you only need to have a handful of third party utilities installed on your path to calculate this metric:

* [Freebayes](https://github.com/ekg/freebayes)
* [Samtools](http://www.htslib.org/)
* [BWA](http://bio-bwa.sourceforge.net/)
* Perl

The script is designed for use on HPCs and has a header that specifies resources for the Slurm job scheduler. Please modify these settings for use on your cluster.

## Preparation to run

The script requires the following input:

1. Your Assembly fasta file
2. A Freebayes VCF file
3. A BWA MEM aligned BAM file

The following are instructions on how to generate these files:

### BWA MEM aligned BAM file

You will need to align whole-genome shotgun reads to your assembly first. We recommend the use of a set of short-reads from the same individual that were not used in the polishing or assembly of the fasta file in the first place. This follows a simplistic principle of data validation, where you are not using the same data used to assemble/polish the assembly to test its validity. In the absence of a separate pool of data, you could use reads from other datasets from the same or related species to test the error rate of the assembly in comparative alignment or resequencing. The choice is yours!

To get started, you need to index your assembly and align the reads like so:

```bash
# Indexing the reference
bwa index your_assembly.fa

# Aligning the reads
bwa mem your_assembly.fa wgs_reads.1.fq wgs_reads.2.fq | samtools sort -T temp -o your_bam.bam - 

# Indexing the BAM file (important!)
samtools index your_bam.bam
```

### Freebayes VCF

Next, you will need to estimate the number of SNP and INDELs present in your assembly as determined by analysis of the alignment of reads in the BAM file. To do this, we use Freebayes:

```bash
# This takes a while!
freebayes -C 2 -0 -O -q 20 -z 0.10 -E 0 -X -u -p 2 -F 0.75 -f your_assembly.fa -v your_output.vcf your_bam.bam
```

At the end, you should have a VCF file, and you'll be read to calculate the QV of your new assembly.

## Running the script

The script is meant to run on the bash shell and has three input parameters. They are the following:

1. A Freebayes VCF file
2. A BAM alignment file
3. The name of your QV output

You can invoke the script on the command line as follows:

```bash
# if not using a cluster
sh calculate_qv.sh your_output.vcf your_bam.bam your_final.qv

# if you are using a Slurm job scheduler
sbatch calculate_qv.sh your_output.vcf your_bam.bam your_final.qv
```

The output file is a plain text file containing only the QV of your assembly. This is a good comparative statistic for identifying the fidelity of your assembly compared to other assembly versions. Every 10 points of difference is an order of magnitude difference in base-pair-level assembly fidelity (e.g. 10 vs 20 is a difference of 90% to 99% base-level accuracy).