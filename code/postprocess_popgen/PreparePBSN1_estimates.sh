#!/bin/bash

#input is file with pairwise Fst for three populations (output of VCFtools --weir-fst-pop for each pair)
#Expected format is chr_location Fst_pop1v2 Fst_pop1v3 Fst_pop2v3
#Assumes space delimited format
#Pop 1 is our population of interest (here, NH)
#To prepare script for pbsn1 calculation: bash PreparePBSN1_estimates.sh > DoPreparePBSN1_estimates.sh
#To Run: bash DoPreparePBSN1_estimates.sh

##NOTE
#We calculate pbsn1 in 5 SNP blocks for each chromosome (can be changed below at NR%5)
#Assumes the number of sites on a chromsome is divisible by number of SNP blocks and that numbers are sequential (ordered files)

##For loop, execute separately per chromosome
for i in {1..19}
do
echo "
###We create summed SNP blocks for each comparison
awk -F' ' '{if(\$1 ~ /^${i}_/){print}}'  NHvGER_NHvIRA_GERvIRAN.weir.filt.sorted.txt  | awk '{s+=\$2}NR%5==0{print s;s=0}' > NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col1.chr${i}.txt
awk -F' ' '{if(\$1 ~ /^${i}_/){print}}'  NHvGER_NHvIRA_GERvIRAN.weir.filt.sorted.txt  | awk '{s+=\$3}NR%5==0{print s;s=0}' > NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col2.chr${i}.txt
awk -F' ' '{if(\$1 ~ /^${i}_/){print}}'  NHvGER_NHvIRA_GERvIRAN.weir.filt.sorted.txt | awk '{s+=\$4}NR%5==0{print s;s=0}' > NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col3.chr${i}.txt
awk -F' ' '{if(\$1 ~ /^${i}_/){print \$1}}' NHvGER_NHvIRA_GERvIRAN.weir.filt.sorted.txt | awk 'NR%5==1' >  NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_colnames.chr${i}.txt


paste -d\" \" NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_colnames.chr${i}.txt NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col1.chr${i}.txt NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col2.chr${i}.txt NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col3.chr${i}.txt > NHvGER_NHvIRA_GERvIRAN.chr${i}.pasted.txt

#rm NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_colnames.chr${i}.txt
#rm NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col1.chr${i}.txt
#rm NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col2.chr${i}.txt
#rm NHvGER_NHvIRA_GERvIRAN.weir.filt.sum5_col3.chr${i}.txt

###we average FST values per block, exclude cases where FST is 1 in any comparison, normalize FST, and calculate PBS and PBSn1
cat  NHvGER_NHvIRA_GERvIRAN.chr${i}.pasted.txt|  awk -F' ' '{print \$1,\$2/5,\$3/5,\$4/5}' |   awk -F' ' '{if(\$2==\"1\" || \$3==\"1\" || \$4==\"1\"){}else{print}}'   | awk -F' ' '{print \$1,-log(1-\$2), -log(1-\$3), -log(1-\$4)}' | awk -F' ' '{print \$1,(\$2+\$3-\$4)/2, (\$2+\$4-\$3)/2, (\$3+\$4-\$2)/2}'  | awk -F' ' '{print \$1,\$2,\$2/(1+\$2+\$3+\$4)}' >NHvGER_NHvIRA_GERvIRAN.chr${i}.pbsn1.txt
"
done

