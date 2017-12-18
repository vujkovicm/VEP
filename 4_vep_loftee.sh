#!/bin/bash

###==============================================
### VEP specific paths and plugins

### the VEP perl script
inPerlScriptPath="../software/VEP/ensembl-vep/vep"

# CACHE FOR LOFTEE IS VEP-83
vepCachePluginPath="/software/VEP/vep_83/vep_cache_and_plugins"

### the path to the human_ancester files
HApath="/software/LOFTEE/human_ancestor.fa.gz"

### the conservation file
CONSpath="/software/LOFTEE/phylocsf.sql"

### LOFTEE path
lofteePath="/software/LOFTEE/loftee"

###==============================================
### Polyphen, HUMVAR

module load samtools-1.1

perl -I $vepCachePluginPath $inPerlScriptPath \
	--dir $vepCachePluginPath \
	--cache \
	--offline \
	--input_file ../output/region.vcf.anno \
	--fasta /project/saleheenlab/snakemake/annotation/human_g1k_v37.fasta \
	--format vcf \
	--everything \
	--force_overwrite \
	--stats_file $outStatsPath \
	--pick \
	--output_file ../output/region.loftee \
	--plugin LoF,human_ancestor_fa:$HApath,conservation_file:$CONSpath,loftee_path:$lofteePath
