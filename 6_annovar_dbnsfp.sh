#!/bin/sh

########### setup input file and output and annovar directory ########### 
DIR=${PWD}
input=$DIR"/../output/region.vcf.anno"

annovar_dir="/software/ANNOVAR/annovar_07162017/annovar"
 
########### run annovar ########### 
cd $annovar_dir

perl convert2annovar.pl -format vcf4 $input -allsample -withfreq > $input.avinput

perl table_annovar.pl   \
 $input.avinput    \
 humandb/   \
 -vcfinput \
 -buildver hg19  \
 -protocol refGene,dbnsfp33a,clinvar_20170905  \
 -operation g,f,f   \
 -outfile $input.avout \
 -otherinfo

cd -
