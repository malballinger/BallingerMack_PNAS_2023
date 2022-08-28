##reading in tables, example
counts <- read.csv("counts_BAT_males.txt",sep="\t",header=TRUE,row.names=1)
condTable <- read.csv("condTable_BAT_males.txt",sep="\t",header=TRUE)
counts2<-as.matrix(counts)


##for trans, using LRT, contrasting FC between hybrid & parent pops
#here's an example -- we can specify w. or without sex + temp as additional variables if you include these individuals in design
#v1
design = ~F1_Parent + population + sex +population:F1_Parent
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
dds <- DESeq(dds, test="LRT", reduced= ~ population + F1_Parent + sex)
#v2
design = ~F1_Parent + population +population:F1_Parent
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
dds <- DESeq(dds, test="LRT", reduced= ~ population + F1_Parent)
res <- results(dds)


##looking at cis
#make sure there is a pseudo-sample column where samples are rep'd across treatments of interest
#pseudo-sample is sample trick from Michael Love, so that M,F or temp info can be incorporated
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
#m = # of individuals
design = ~temp + temp:pseudo_sampl + temp:allele
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
m <- 6
#we don't estimate size factors here b/c this is ASE
sizeFactors(dds) <- rep(1, 2*m)
dds <- DESeq(dds, fitType="parametric")
resultsNames(dds)


##once these are combined (e.g., need log fold & padj values between parents, F1s, and then LRT for trans component), sort out with AWK
##Note: we adjusted FDR for cases where genes were tested by DESeq2 for one test but not the other (i.e., FDR estimates for parents but not hybrid) -- this is done with p.adjust, method="BH"

#0.05 cutoff
awk -F' ' '{if($3<0.05 && $5<0.05 && $7>0.05){print $0,"CIS_only"}else{print}}' | awk -F' ' '{if($3>0.05 && $5<0.05 && $7<0.05){print $0,"Compensatory"}else{print  $0}}'| awk -F' ' '{if($3>0.05 && $5>0.05 && $7>0.05){print $0,"Conserved"}else{print $0}}' | awk -F' ' '{if($3<0.05 && $5>0.05 && $7<0.05){print $0,"Trans"}else{print  $0}}' |awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4>0) || ($2<0 && $4<0) ) && ( sqrt($2^2)>sqrt($4^2)) ){print $0,"cis+trans_same"}else{print}}' | awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4>0) || ($2<0 && $4<0) ) && ( sqrt($2^2)<sqrt($4^2)) ){print $0,"cis+trans_opp"}else{print}}' | awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4<0) || ($2<0 && $4>0) )){print $0,"cisxtrans"}else{print}}' | awk -F' ' '{if($3=="NA" || $5=="NA" || $7=="NA"){print $1,$2,$3,$4,$5,$6,$7,"NoPower"}else{print}}' | awk -F' ' '{if($8==""){print $0,"Ambiguous"}else{print}}' 

#add the following to condense cis+trans groups for plotting in R
 | awk -F' ' '{if($8=="Conserved" || $8=="Ambiguous"){print $1,$2,$3,$4,$5,$6,$7,"Ambiguous_orConserved"}else{print}}'  > categories.txt
