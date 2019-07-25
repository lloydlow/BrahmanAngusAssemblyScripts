#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=1
#SBATCH -p short

# $1 = vcf file
# $2 = input bam file
# $3 = output QV file

module load samtools

NUM_BP=`samtools depth $2 | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); if(scalar(@s) >= 3){$c++;}} print "$c\n";'`
echo "num bp: "$NUM_BP

NUM_SNP=`cat $1 |grep -v "#" | awk -F "\t" '{if (!match($NF, "0/1")) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' | tr ';' ' ' | sed s/AB=//g | awk -v WEIGHT=0 '{if ($6 >= WEIGHT) print $0}' | awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`
echo "num snp: "$NUM_SNP

perl -e 'chomp(@ARGV); $ns = $ARGV[0]; $nb = $ARGV[1]; print (-10 * log($ns/$nb)/log(10)); print "\n";' $NUM_SNP $NUM_BP > $3
cat $3
