import synapseclient as syn
import pandas as pd


import os

allfiles = [a for a in os.listdir(".") if 'protData-' in a]

tablist = []

for a in allfiles:
    #read in file
    tab = pd.read_csv(a, sep='\t', index_col=0)
    tab['Cell type'] = tab.index
    tab = pd.melt(tab, id_vars='Cell type')
    #append algorithm matrix
    [dat, alg, mat] = os.path.splitext(a)[0].split('-')
    #add to tablist
    tab = tab.rename({'variable':'sample', 'value':'cellPop'}, axis=1)
    tab['algorithm'] = alg
    tab['matrix'] = mat
    tablist.append(tab)


fulltab = pd.concat(tablist)
print(fulltab.head())
fulltab.to_csv('full_tab.csv')
sync = syn.login()

projid = "syn22128879"
table = syn.table.build_table('Tumor deconvolution results',projid,'full_tab.csv')
sync.store(table)
