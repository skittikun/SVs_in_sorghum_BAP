# SVs_in_sorghum_BAP

The SVs were called by svtools as described in https://github.com/hall-lab/svtools/blob/master/Tutorial.md. The R custom codes providing here were used to 

1. convert sv vcf format into 0 (absent) and 1 (present) SVs via "modify_lumpy.R" 
2. calculate 500-kb windowed SNP and SV diversity via "sv_snp_diversity.R"
3. plot a genome-wide circos via "8k_circos_5group.R"
4. plot global site frequency spectrum via "site_frequency_spectrum"
5. plot linkage disequilibrium decay via "linkage.R"
6. perform k-mean clustering via "k-clustering.R"
7. plot the k-mean clusters with geological location via "map_distribution.R"
8. plot ADMIXTURE result with k-mean clustering via "plot_admixture_347.R"
9. determine cluster-specific SVs via "unique_SV_group.R"
10. perform GO enrichment analysis via "GO_general.R", and "GO_cluster_specific.R"

