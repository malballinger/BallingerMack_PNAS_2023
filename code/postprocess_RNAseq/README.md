## The majority of the scripts housed here are used to analyze and plot various RNA-seq results.
The script [clean_parentalReadCounts.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/clean_parentalReadCounts.R) takes raw read count data (produced via HTSeq) and generates matrices used in [get_DE_parents_analysis.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/get_DE_parents_analysis.R) to model differential expression via DESEQ2.

- Figure 1D: [plot_sexPCA.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_sexPCA.R)
- Figure 1E: [plot_parental_DEdivergence.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_parental_DEdivergence.R)
- Figure 2A: [plot_parents_ternary_GxE.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_parents_ternary_GxE.R)
- Figure 2B:[plot_parental_GxE.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_parental_GxE.R)
- Figure 2C:[plot_direction_plasticity.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_direction_plasticity.R)
- Figure 3:[plot_cistrans.R](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_RNAseq/plot_cistrans.R)