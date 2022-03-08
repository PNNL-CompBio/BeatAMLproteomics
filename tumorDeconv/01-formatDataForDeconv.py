'''
takes the python beataml data and writes to filese
then runs the cwl framework
'''


###get data write to filese
import sys

sys.path.append("..")
from pybeataml import load_data
import pandas as pd

data = load_data.AMLData()
mat = pd.pivot(data.proteomics, index='gene_symbol', \
               columns='sample_id', values='exp_value')
mat.to_csv('AML-tumor-prot-raw.tsv', sep='\t')

##write out proteomics files

mat = pd.pivot(data.rna, index='gene_symbol', \
               columns='sample_id', values='exp_value')
mat.to_csv('AML-tumor-mrna-raw.tsv', sep='\t')
#bdb_file = 'http://'
#sig_mat_file = ''
