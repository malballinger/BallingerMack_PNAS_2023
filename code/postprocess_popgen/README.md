## 1. Genetic PCA (Figure 4A) was generated via [SNPrelate](https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html)

## 2. Autosomal scans for selection were performed via [_PBSn1_](https://github.com/malballinger/BallingerMack_NYBZase_2022/blob/main/code/postprocess_popgen/PreparePBSN1_estimates.sh)
> Results of *PBSn1* (Figure 4B) were plotted using the package [qqman](https://r-graph-gallery.com/101_Manhattan_plot.html)
```R
manhattan(pbsn1_out,p = "zscore", logp = FALSE, main = " ", ylim = c(0, 1), cex = 0.6, cex.axis = 0.9, 
    col = c("#B0C5CA", "#A4A4A4"), suggestiveline = F, genomewideline = F, highlight = ASE_genes)
```
