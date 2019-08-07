# Scaffolding with optical map and HiC
---

## Optical map-based scaffolding
The methodology of how the scaffolds were built using optical map has been sufficiently described in the manuscript's Supplementary Information. So, this part will not be repeated here. The detail for the steps on how to process the scaffolds after getting the haplotype-resolved scaffolds is given here. As the optical map essentially only revealed locations of the 6 bp motif, whenever chimera contigs were identified, it was hard to determine the the precise break points. Where possible, we attemped to adjust the breakpoints. Briefly the procedures involved adjusting the break positions introduced to the contigs, which was done by two methods (i) adjust based on the precise drop in Illumina WGS short read coverage, and (ii) adjust based on breakpoint identified from aligning the haplotigs to the other breed haplotigs. 

* brahman_angus_bionano_align_orderSNP.R
* brahman_angus_bionano_check_agp.R
* brahman_angus_bionano_corrected_cut.R
* brahman_angus_bionano_cov_around_cutpt.R
* brahman_angus_bionano_cutpt.R
* brahman_angus_bionano_indi_alignProbes.R
* brahman_angus_bionano_mashmap.R
* brahman_angus_bionano_order_scaffold.R
* brahman_angus_bionano_order_scaffold2.R
* brahman_angus_bionano_order_scaffold2b.R
* brahman_angus_bionano_order_scaffold3.R

## Hi-C-based scaffolding
Three different scaffolders, 3D-DNA, Proximo (Phase Genomics) and SALSA2 were evaluated for building scaffolds and the input for scaffolding was PacBio-based contigs. Detail parameters used in each scaffolding programs were given in the Supplementary Information. The following scripts describe how the Hi-C scaffolds were ordered.

* brahman_angus_align_orderSNP.R
* brahman_angus_order_scaffold.R
* brahman_angus_order_scaffold2.R
* brahman_angus_order_scaffold2b.R
* brahman_angus_order_scaffold3.R

