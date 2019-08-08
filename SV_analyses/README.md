# SV analyses
The count of SVs partitioned by their sizes and how these SVs overlapped between Brahman and Angus were analyzed. The basic workflow here is nucmer alignment of the two cattle assemblies broken at sequence gaps, followed by SV identification in Assemblytics. 

## nucmer alignment
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1            # number of nodes 
#SBATCH -n 4            # number of cores
#SBATCH --time=08:00:00 # time allocation (D-HH:MM:SS)
#SBATCH --mem=120GB       # memory pool for all cores

module load MUMmer/4.0.0beta-foss-2016b

#arsucd as the ref and align angus to it
nucmer --maxmatch -t 4 -l 100 -c 500 cattle_arsucd_chr_only_ungapped.fa bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil_chronly_ungapped.fa -p ARSUCD_VS_ANGUS
delta-filter -g ARSUCD_VS_ANGUS.delta > ARSUCD_VS_ANGUS.global.delta

#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1            # number of nodes 
#SBATCH -n 1            # number of cores
#SBATCH --time=05:00:00 # time allocation (D-HH:MM:SS)
#SBATCH --mem=40GB       # memory pool for all cores

module load MUMmer/4.0.0beta-foss-2016b

show-coords -r -c -l -T -H ARSUCD_VS_ANGUS.delta > ARSUCD_VS_ANGUS.delta.coords
show-coords -r -c -l -T -H ARSUCD_VS_ANGUS.global.delta > ARSUCD_VS_ANGUS.global.coords

```

Below are the relevant scripts for this analysis after getting Assemblytics output.
* brahman_angus_SV_Assemblytics.R
* brahman_angus_SV_Assemblytics_2.R
* brahman_angus_SV_Assemblytics_3.R
* brahman_angus_SV_Assemblytics_4.R

## CNV calls and statistics

Now we will download the SRA datasets for CNV calling using JaRMS and Lumpy. These will be queued up as sequential, simultaneous, tasks.

```bash
# The SRA read splitting doesn't work with BWA because of non-unique pair names! Fixing...
perl -lane 'print $F[2];' < beef_panel_pipeline_spreadsheet.tab | xargs -I {} sbatch --nodes=1 --mem=3000 --ntasks-per-node=1 -p short --wrap="python3 ~/python_toolchain/sequenceData/fixSRAFastqFiles.py -f {}_1.fastq.gz -r {}_2.fastq.gz -o {}.fmt -l {}.log"

perl -lane 'for($x = 0; $x < 2; $x++){$F[$x] =~ s/\_([12])/\.fmt\.$1/; $F[$x] =~ s/fastq/fq/;} print join("\t", @F);' < beef_panel_pipeline_spreadsheet.fullp.tab > beef_panel_pipeline_spreadsheet.fullp.rfmt.tab
```

Next, we need to align the reads to each assembly fasta to begin the SV calling.

```bash
# Angus
python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b angus_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -m -p msn

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b brahman_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -m -p msn

python3 ~/python_toolchain/sequenceData/slurmAlignScriptBWA.py -b arsucd_asm -t beef_panel_pipeline_spreadsheet.fullp.rfmt.tab -f /beegfs/project/bostauruscnv/assembly/ARS-UCD1.2_Btau5.0.1Y.fa -m -p msn

mkdir angus_jarms; module load java/1.8.0_121; for i in angus_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_angus_bionano_NCBI_full_corrected_gapfill_arrow_fil.fasta -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./angus_jarms/"; done
mkdir angus_lumpy; module load bwa samtools; for i in angus_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./angus_lumpy/;"; done

mkdir brahman_jarms; module load java/1.8.0_121; for i in brahman_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /project/cattle_genome_assemblies/angusxbrahman/asms/bostaurus_brahma_bionano_NCBI_full_corrected_gapfill_arrow_fil_withM.fasta -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./brahman_jarms/"; done
mkdir brahman_lumpy; module load bwa samtools; for i in brahman_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./brahman_lumpy/;"; done

mkdir arsucd_jarms; module load java/1.8.0_121; for i in arsucd_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=4 --mem=20000 -p short --wrap="java -Xmx19g -jar ~/rumen_longread_metagenome_assembly/binaries/JaRMS/store/JaRMS.jar call -i $i -f /beegfs/project/bostauruscnv/assembly/ARS-UCD1.2_Btau5.0.1Y.fa -o $name.JaRMs -t 4 -m 10000000; mv $name.JaRMs* ./arsucd_jarms/"; done
mkdir arsucd_lumpy; module load bwa samtools; for i in arsucd_asm/*/*.merged.bam; do name=`echo $i | cut -d'/' -f2`; echo $name; sbatch --nodes=1 --ntasks-per-node=1 --mem=15000 -p short --wrap="lumpyexpress -B $i -o $name.vcf; mv $name.vcf ./arsucd_lumpy/;"; done

# OK, now to convert to BEDPE
for i in *_lumpy/*.vcf; do echo $i; sbatch --nodes=1 --ntasks-per-node=1 --mem=5000 -p msn --wrap="/home/derek.bickharhth/lumpy-sv/scripts/vcfToBedpe -i $i -o $i.bedpe"; done

# First, let's try to generate counts for each event just to see if raw counts are less on the "native" assembly for each breed.
# I  need to generate a list for each file, but then I think that I can push the data through my tabfile column counter script
module load perl/5.24.1
for i in Angus Brahman Gelbvieh Hereford RedAngus Shorthorn Simmental; do for j in `seq 1 6`; do echo "$i$j"; done; done > combinations.list

mkdir raw_counts
for i in `cat combinations.list`; do echo $i; perl ~/rumen_longread_metagenome_assembly/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f angus_lumpy/$i.vcf.bedpe,brahman_lumpy/$i.vcf.bedpe,arsucd_lumpy/$i.vcf.bedpe -c 10 -e '#' -o raw_counts/$i.lumpy.counts.table; done

# I wrote a quick compilation script to generate raw counts from my table files
cd raw_counts/
perl calculate_averages.pl Angus
```

#### calculate_averages.pl

```perl
#!/usr/bin/perl
# This is a one-off script designed to process tabFileColumnCounter output from the CNV files

use strict;
my $usage = "perl $0 <prefix of file>\n";
chomp(@ARGV);
unless(scalar(@ARGV) == 1){
        print $usage;
        exit;
}

my @files = `ls $ARGV[0]*`;
chomp(@files);

my %data; # {SVtype}->{asm}->{file} = count
foreach my $f (@files){
        open(my $IN, "< $f");
        for(my $x = 0; $x < 5; $x++){
                <$IN>;
        }
        while(my $line = <$IN>){
                chomp $line;
                my @segs = split(/\t/, $line);
                for(my $x = 0; $x < scalar(@segs); $x += 2){
                        $data{$segs[$x]}->{"asm$x"}->{$f} = $segs[$x + 1] * 1;
                }
        }
        close $IN;
}

print "Asm\tSV\tsum\tavgnum\tstdev\n";
foreach my $svs (sort {$a cmp $b} keys(%data)){
        if($svs eq "" || length($svs) == 0){next;}
        foreach my $asm (sort {$a cmp $b} keys(%{$data{$svs}})){
                my $count = 0;
                my $sum = 0;
                foreach my $file (sort {$a cmp $b} keys(%{$data{$svs}->{$asm}})){
                        $count += 1;
                        $sum += $data{$svs}->{$asm}->{$file};
                }
                if($count == 0){
                        print "$asm\t$svs\t0\t0\t0\n";
                }else{
                        my $avg = $sum / $count;
                        if($count == 1){
                                $avg = sprintf("%.3f", $avg);
                                print "$asm\t$svs\t$sum\t$avg\t0\n";
                                next;
                        }
                        my $ss = 0;
                        foreach my $file (sort {$a cmp $b} keys(%{$data{$svs}->{$asm}})){
                                $ss += ($data{$svs}->{$asm}->{$file} - $avg)**2;
                        }
                        $avg = sprintf("%.3f", $avg);
                        my $stdev = sprintf("%.3f", sqrt($ss / ($count - 1)));
                        print "$asm\t$svs\t$sum\t$avg\t$stdev\n";
                }
        }
}
```


## VST analysis and gene coordinate conversion

Now we will apply the VST statistic to our dataset to identify copy number variable regions of the reference genome that are divergent between the taurine and indicine cattle. 

We will use the Brahman individuals as an outgroup for this comparison, and calculate the Vst on each respective assembly.

```bash
# Converting the gtf and NCBI files into bed files for coordinate intersections
module load bedtools/2.25.0

perl -e 'while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); @bsegs = split(/\"/, $s[8]); if($s[2] =~ /gene/){print "$s[0]\t$s[3]\t$s[4]\t$bsegs[1]\n";}}' < sire.UOA_angus_1.96.gtf | bedtools sort -i stdin > sire.UOA_angus_1.96.gene.bed
perl -e 'while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); @bsegs = split(/\"/, $s[8]); if($s[2] =~ /gene/){print "$s[0]\t$s[3]\t$s[4]\t$bsegs[1]\n";}}' < dam.UOA_brahman_1.96.gtf | bedtools sort -i stdin > dam.UOA_brahman_1.96.gene.bed

dos2unix ncbi_ars_ucd_1.2_refseq.txt
perl -e '<>; while(<>){chomp; @s = split(/\t/); if($s[12] eq "" || $s[13] eq ""){next;} print "$s[10]\t$s[12]\t$s[13]\t$s[5]\n";}' < ncbi_ars_ucd_1.2_refseq.txt | bedtools sort -i stdin > ncbi_ars_ucd_1.2_refseq.bed
perl -e 'while(<>){chomp; if($_ =~ /^\#/){next;} @s = split(/\t/); @bsegs = split(/\"/, $s[8]); if($s[2] =~ /gene/){print "$s[0]\t$s[3]\t$s[4]\t$bsegs[1]\n";}}' < bos_taurus.ars-ucd1.2.96.chr.gtf | bedtools sort -i stdin > bos_taurus.ars-ucd1.2.96.chr.gene.bed

echo -e 'dam.UOA_brahman_1.96.gene.bed\tensembl' > dam.UOA_brahman_1.96.gene.list
echo -e 'sire.UOA_angus_1.96.gene.bed\tensembl' > sire.UOA_angus_1.96.gene.list
echo -e 'ncbi_ars_ucd_1.2_refseq.bed\trefseq' > ncbi_ars_ucd_1.2_refseq.list
# I also added in the ensembl annotation for ARS-UCDv1.2

module load java/1.8.0_121
sbatch annotateScript.sh angus sire.UOA_angus_1.96.gene.list angus.cnvrs
sbatch annotateScript.sh brahman dam.UOA_brahman_1.96.gene.list brahman.cnvrs
sbatch annotateScript.sh arsucd ncbi_ars_ucd_1.2_refseq.list arsucd.cnvrs

head -n 1 brahman.cnvrs_windows_ensembl.tab | perl -lane 'for($x = 6; $x < scalar(@F); $x++){print "$F[$x]\t1";}' > pop_list_base.tab

cp pop_list_base.tab brahman_pop_list.tab
cp pop_list_base.tab angus_pop_list.tab
cp pop_list_base.tab arsucd_pop_list.tab

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c brahman.cnvrs_windows_ensembl.tab -o brahman.cnvrs_windows_vst_genes.bed
perl calculate_vst_differences_cn.pl -c brahman.cnvrs_windows_ensembl.tab -p brahman_pop_list.tab -o brahman.vst_test.bed

# I just updated the script to print out a melted file for plotting VST stats
# I also added a minimum filter for differences in CN count between averages
# the Brahman file will remain the same
python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c angus.cnvrs_windows_ensembl.tab -o angus.cnvrs_windows_vst_genes.tvi.bed -m angus.cnvrs_windows_vst_genes.tvi.melt -g 1.5

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c arsucd.cnvrs_windows_refseq.tab -o arsucd.cnvrs_windows_vst_genes.tvi.bed -m arsucd.cnvrs_windows_vst_genes.tvi.melt -g 1.5

python3 ~/python_toolchain/sequenceData/calculateVstDifferences.py -p brahman_pop_list.tab -c arsucd.cnvrs_windows_ensembl.tab -o arsucd.cnvrs_windows_vst_ensembl.bed -m arsucd.cnvrs_windows_vst_ensembl.melt -g 1.5
Melt    ENSBTAG00000023157      0.3742902428652937      1.8323636363636364
Melt    ENSBTAG00000053555      0.3865116874761947      1.7015151515151512
Melt    ENSBTAG00000053557      0.3045985843651203      1.5839151515151513
Melt    ENSBTAG00000033545      0.28134455365978633     1.5851090909090906
Melt    ENSBTAG00000053028      0.29194067204914553     1.621048484848485
Melt    ENSBTAG00000052940      0.7684926833817546      1.6270787878787878
Melt    ENSBTAG00000015575      0.9506884029779913      2.391454545454546
Melt    ENSBTAG00000037937      0.425415503351496       1.7660545454545453
Melt    ENSBTAG00000052878      0.5487232336595422      2.7513272727272735
Melt    ENSBTAG00000048030      0.226623521102344       1.5979272727272724
Melt    ENSBTAG00000013957      0.2068551932767507      1.647064516129033
Melt    ENSBTAG00000030838      0.25243450079686874     47.670748387096765
Melt    ENSBTAG00000001925      0.9079146432602809      1.7279090909090908
Melt    ENSBTAG00000050853      0.4019145213065941      7.678957575757577
Melt    ENSBTAG00000007649      0.26534987091180695     2.578490322580645
Dealt with 37 null fields
```

And here is the script I'm using to call the variants.

```bash
#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=10000
#SBATCH -p msn
# $1 = type
# $2 = annotation_db file
# $3 = output basename

lumpy=${1}"_lumpy"
jarms=${1}"_jarms"

cndata=${1}".cnlist"
cnvs=${1}".cnvs"

# Process bedpe file
for i in `ls $lumpy/*.bedpe`
do
        grep -v 'BND' $i | perl -e '<>; while(<>){chomp; @s = split(/\t/); print "$s[0]\t$s[1]\t$s[5]\n";' > $i.bed
        python3 ~/python_toolchain/utils/tabFileColumnCounter.py -f $i -c 10 -d '\t' > $i.stats
done

# Create cnv and cn lists
for i in `ls $lumpy/*.bed`
do
        sample=`basename $i | cut -d'.' -f1`
        echo -e "$i\t$sample"
done > $cnvs

for i in `ls $jarms/*.bed.levels`
do
        sample=`basename $i | cut -d'.' -f1`
        echo -e "$i\t$sample"
done > $cndata

java -Xmx10g -jar ~/rumen_longread_metagenome_assembly/binaries/AnnotateUsingGenomicInfo/store/AnnotateUsingGenomicInfo.jar -d $2 -i $cnvs -c $cndata -o $3 -t
```