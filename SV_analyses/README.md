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

