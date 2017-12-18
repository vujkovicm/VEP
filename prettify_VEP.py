#!~/miniconda3/bin/python

import os
import sys
import gzip
import re
import string
import pandas as pd
import json


###=================================
### test values for parameters

if False:
	in_file_path = '../output/region.loftee'
	out_file_path = '../output/test_main.tab'

###=================================
### read parameters from command line

if True:
        in_file_path = sys.argv[1]
        out_file_path = sys.argv[2]

###=================================

#if os.path.exists(out_file_path) and out_file_path != in_file_path:
#        os.remove(out_file_path)

#in_file_path = "/project/saleheenlab/snakemake/VEP/1_vep/output/test.txt"
#out_file_path = '../output/test_main.tab'

# import table
df = pd.read_table(in_file_path, sep = "\t", comment = "#", header = None)
df.columns = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra']

# split dataframe into fixed and tmp
tmp = df[['Extra']].copy()
df.drop('Extra', axis = 1, inplace = True)

# convert to dictionary-looking string
tmp = pd.DataFrame(tmp.Extra.str.replace('=', '" : "'))
tmp = pd.DataFrame(tmp.Extra.str.replace(';', '", "'))
tmp['Extra'] = '{"' + tmp['Extra'].astype(str) + '"}'

# transform dictionary-looking string into to an indexed dataframe
tmp = pd.DataFrame([json.loads(line.strip()) for line in tmp.Extra])

# merge fixed with variable columns
df = pd.concat([df, tmp], axis = 1)

df.to_csv(out_file_path, sep='\t', na_rep=".", index = False, header = True)

