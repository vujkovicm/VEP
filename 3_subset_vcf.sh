#!/bin/bash

BED="../output/gene.bed"
VCF="/data/exome_dataset_22k_vcf.bgz"
FA="../meta/human_g1k_v37.fasta"
TMP="../output/region-tmp.vcf"
OUT="../output/region.vcf"

# norm: Left-align and normalize indels
# -m both: Merge, both SNP and indel records can be multiallelic
# -R regions file

# left-normalize and multi-allelic split
bcftools norm -Ou -m -both -R $BED -o $TMP $VCF 

# align to to reference
bcftools norm -f $FA -o $OUT $TMP

# rename SNP/rsids to chr:bp:ref:alt
bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -o $OUT.anno $OUT

# now all the variants are already bi-allelic so plink has no problem
/project/saleheenlab/software/PLINK/PLINK_1.9/plink --vcf $OUT.anno --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr --split-x b37 no-fail --make-bed --out ../output/region
/project/saleheenlab/software/PLINK/PLINK_1.9/plink --bfile ../output/region --freq --missing --hardy --out ../output/region
/project/saleheenlab/software/PLINK/PLINK_1.9/plink --bfile ../output/region --recode AD --out ../output/region
