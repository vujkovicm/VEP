#!/bin/bash

# dichotomous phenotype
/software/PLINK/PLINK_1.9/plink --bim ../output/region.bim --bed ../output/region.bed --fam ../input/mi.tab --hardy --out ../output/counts_MI
/software/PLINK/PLINK_1.9/plink --bim ../output/region.bim --bed ../output/region.bed --fam ../input/dm.tab --hardy --out ../output/counts_DM

# continuous phenotype
/software/PLINK/PLINK_1.9/plink --bfile ../output/region --pheno ../input/pheno.tab --allow-no-sex --pheno-name ldlc --assoc qt-means --out ../output/means_LDL
/software/PLINK/PLINK_1.9/plink --bfile ../output/region --pheno ../input/pheno.tab --allow-no-sex --pheno-name hdlmgdl --assoc qt-means --out ../output/means_HDL
/software/PLINK/PLINK_1.9/plink --bfile ../output/region --pheno ../input/pheno.tab --allow-no-sex --pheno-name tgmgdl --assoc qt-means --out ../output/means_TG
/software/PLINK/PLINK_1.9/plink --bfile ../output/region --pheno ../input/pheno.tab --allow-no-sex --pheno-name cholesterolmgdl --assoc qt-means --out ../output/means_TC
