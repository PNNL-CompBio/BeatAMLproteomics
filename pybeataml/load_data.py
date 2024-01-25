"""
First pass at data centered class.
Ideally a single class instance can handle loading of the entire dataset.
TODO
    - Add cache mechanism for quick storage/loading rather than going to synap
    - Think of a way to get current feature cols and drug cols from
    a subset data, add functionality
"""
import os

import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split

from pybeataml.data import ExperimentalData
from pybeataml.load_data_from_synpase import load_table, load_file, load_excel

# colors for syncing across plots/R/python
cluster_colors = list(sns.color_palette("Dark2", 8)[1:3]) + \
                 list(sns.color_palette("Dark2", 8)[4:6])
data_colors = ["#5495CF", "#F5AF4D", "#DB4743", "#7C873E", "#FEF4D5"]

# current synapse ids, check with Camilo to see if these are the final (
# I know there are some other corrected/v2/uncorrected in the R code)
global_id = 'syn25808020'
phospho_id = 'syn26477193'  # syn25808662
rnaseq_id = 'syn26545877'
drug_response_id = 'syn25830473'
meta_file_id = 'syn26534982'
wes_id = 'syn26428827'
clusters_id = 'syn26642544'
clinical_summary_id = 'syn25796769'
metabolomics_id = 'syn52224584'
lipidomics_id = 'syn52121001'

def prep_rnaseq():
    f_name = 'data/rna.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name)
    cols = ['display_label', 'labId',
            'RNA counts']

    mapper = {
        'display_label': 'gene_symbol',
        'RNA counts': 'exp_value',
        'labId': 'sample_id',
    }
    data = load_table(rnaseq_id)
    subset = data.loc[:, cols]
    subset.rename(mapper, axis=1, inplace=True)
    subset['source'] = 'rna_seq'
    subset['label'] = subset.gene_symbol + '_rna'
    if not os.path.exists(f_name):
        subset.to_csv(f_name)
    return subset

def prep_metabolomics():
    cols = ['Name', 'labId',
            'area']

    mapper = {
        'display_label': 'Name',
        'area': 'area',
        'labId': 'labId',
    }

    # import HILIC pos & neg and drop extra rows & columns
    data_pos = load_excel(metabolomics_id, 0)
    data_pos = data_pos.iloc[:-135] # drop unknowns
    data_pos = data_pos.drop(columns=['Blank_BEAT_AML_01_HILIC_POS',
                                      'Blank_BEAT_AML_02_HILIC_POS',
                                      'Blank_BEAT_AML_03_HILIC_POS',
                                      'Blank_BEAT_AML_04_HILIC_POS',
                                      'Blank_BEAT_AML_05_HILIC_POS']) # drop blanks
    data_pos = data_pos.drop(columns=['m/z', 'RT [min]', 'Tags',
                              'Standardized name', 'Super class',
                              'Main class', ' Sub class', 'Formula',
                              'Annot. DeltaMass [ppm]',
                              'Annotation MW', 'Reference Ion'])

    data_neg = load_excel(metabolomics_id, 1)
    data_neg = data_neg.iloc[:-98] # drop unknowns
    data_neg = data_neg.drop(columns=['Blank_BEAT_AML_01_HILIC_NEG_2uL_18Apr23_Olympic_WBEH-8588_r1.raw (F169)',
                                      'Blank_BEAT_AML_01_HILIC_NEG_2uL_18Apr23_Olympic_WBEH-8588_r1_20230418144054.raw (F170)',
                                      'Blank_BEAT_AML_02_HILIC_NEG',
                                      'Blank_BEAT_AML_02_HILIC_NEG2',
                                      'Blank_BEAT_AML_02_HILIC_NEG3',
                                      'Blank_BEAT_AML_03_HILIC_NEG',
                                      'Blank_BEAT_AML_04_HILIC_NEG',
                                      'Blank_BEAT_AML_05_HILIC_NEG']) # drop blanks
    data_neg = data_neg.drop(columns=['m/z', 'RT [min]', 'Tags',
                              'Standardized name', 'Super class',
                              'Main class', ' Sub class', 'Formula',
                              'Annot. DeltaMass [ppm]',
                              'Annotation MW', 'Reference Ion'])

    # reformat to long format, normalize, and combine pos & neg data
    data_pos = pd.melt(data_pos, id_vars=['Name'], 
                        var_name = 'labId', value_name='area')
    data_pos['area'] = data_pos['area'] / data_pos['area'].abs().max()

    data_neg = pd.melt(data_neg, id_vars=['Name'], 
                        var_name = 'labId', value_name='area')
    data_neg['area'] = data_neg['area'] / data_neg['area'].abs().max()

    data = pd.concat([data_pos, data_neg])

    # extract sample IDs from labID column
    data['labId_og'] = data['labId']
    data['labId'] = data['labId_og'].apply(lambda st: st[st.find("BEAT_AML_PNL_") + 1:st.find("_M")])
    data['labId'] = pd.to_numeric(data['labId'], errors = 'coerce').astype(pd.Int16Dtype())
    
    # reformat data
    subset = data.loc[:, cols]
    subset.rename(mapper, axis=1, inplace=True)
    subset['source'] = 'metabolomics'
    subset['label'] = subset.Name + '_met'
    return subset

def prep_lipidomics():
    cols = ['Metabolite name', 'labId',
            'area']

    mapper = {
        'display_label': 'Metabolite name',
        'area': 'area',
        'labId': 'labId',
    }

    # import pos & neg and drop extra columns
    data_pos = load_excel(lipidomics_id, 1)
    data_pos = data_pos.drop(columns=['CPTAC4_AML_BM_L_QC_01_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_02_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_03_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_04_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_05_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_06_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_07_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_01_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_02_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_03_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_04_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_05_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_06_Lumos_Pos_18Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_07_Lumos_Pos_18Feb23_Crater-WCSH315305']) # drop CPTAC4
    data_pos = data_pos.drop(columns=['Alignment ID', 'Average Rt(min)',
                                      'Average Mz', 'Adduct type',
                                      'Reference m/z', 'Formula', 
                                      'Ontology', 'MS/MS spectrum'])

    data_neg = load_excel(lipidomics_id, 0)
    data_neg = data_neg.drop(columns=['CPTAC4_AML_BM_L_QC_01_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_02_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_03_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_04_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_05_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_06_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_BM_L_QC_07_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_01_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_02_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_03_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_04_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_05_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_06_Lumos_Neg_22Feb23_Crater-WCSH315305',
                                      'CPTAC4_AML_WB_L_QC_07_Lumos_Neg_22Feb23_Crater-WCSH315305']) # drop CPTAC4
    data_neg = data_neg.drop(columns=['Alignment ID', 'Average Rt(min)',
                                      'Average Mz', 'Adduct type',
                                      'Reference m/z', 'Formula', 
                                      'Ontology', 'MS/MS spectrum'])
    
    # average across duplicate compound names
    data_pos = data_pos.groupby(['Metabolite name'], as_index = False).mean()
    data_neg = data_neg.groupby(['Metabolite name'], as_index = False).mean()

    # reformat to long format, normalize, and combine pos & neg data
    data_pos = pd.melt(data_pos, id_vars=['Metabolite name'], 
                        var_name = 'labId', value_name='area')
    data_pos['area'] = data_pos['area'] / data_pos['area'].abs().max()

    data_neg = pd.melt(data_neg, id_vars=['Metabolite name'], 
                        var_name = 'labId', value_name='area')
    data_neg['area'] = data_neg['area'] / data_neg['area'].abs().max()

    data = pd.concat([data_pos, data_neg])

    # extract sample IDs from labID column
    data['labId_og'] = data['labId']
    data['labId'] = data['labId_og'].apply(lambda st: st[st.find("BEAT_AML_PNL_") + 1:st.find("_L")])
    data['labId'] = pd.to_numeric(data['labId'], errors = 'coerce').astype(pd.Int16Dtype())
    
    # reformat data
    subset = data.loc[:, cols]
    subset.rename(mapper, axis=1, inplace=True)
    subset['source'] = 'lipidomics'
    subset['label'] = subset['Metabolite name'] + '_lip'
    return subset

def prep_phosph():
    f_name = 'data/phospho.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name)
    pho_cols = ['Gene', 'SiteID', 'LogRatio',
                'SampleID.full', 'Barcode.ID']
    phosp_mapper = {
        'Gene': 'gene_symbol',
        'SiteID': 'label',
        'LogRatio': 'exp_value',
        'SampleID.full': 'sample_id_full',
        'Barcode.ID': 'sample_id',
    }
    phospho_data = load_table(phospho_id)
    phosph_subset = phospho_data.loc[:, pho_cols]
    phosph_subset.rename(phosp_mapper, axis=1, inplace=True)
    phosph_subset['source'] = 'phospho'
    if not os.path.exists(f_name):
        phosph_subset.to_csv(f_name)
    return phosph_subset


def prep_proteomics():
    f_name = 'data/global.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name)
    proteomics_mapper = {
        'Gene': 'gene_symbol',
        'SiteID': 'label',
        'LogRatio': 'exp_value',
        'SampleID.full': 'sample_id_full',
        'Barcode.ID': 'sample_id',
    }
    global_data = load_table(global_id)
    # remove empty gene columns? Is this safe
    proteomics = global_data.loc[~global_data.Gene.isna(), :].copy()
    proteomics.rename(proteomics_mapper, axis=1, inplace=True)
    pho_cols = ['gene_symbol', 'exp_value',
                'sample_id_full', 'sample_id']
    proteomics = proteomics.loc[:, pho_cols]

    # add source and label column for MAGINE
    proteomics['label'] = proteomics.gene_symbol + '_prot'
    proteomics['source'] = 'proteomics'
    if not os.path.exists(f_name):
        proteomics.to_csv(f_name)
    return proteomics


def load_drug_response():
    f_name = 'data/drug_response.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name)
    response_data = load_table(drug_response_id)
    response_data.auc = response_data.auc.astype(float)
    response_data.aic = response_data.aic.astype(float)
    response_data.deviance = response_data.deviance.astype(float)

    # create this column as filter column
    response_data[
        'new_col'] = response_data.proteomic_lab_id + '_' + response_data.inhibitor
    to_remove = []
    # If multiple measurements, remove if std(auc) >50 (seems to be skewed
    # to remove larger auc)
    for i, j in response_data.groupby(['proteomic_lab_id', 'inhibitor']):
        if j.shape[0] > 1:
            if j['auc'].std() > 50:
                to_remove.append('_'.join(i))

    response_data = response_data.loc[
        ~response_data.new_col.isin(to_remove)].copy()
    response_data = response_data.groupby(
        ['proteomic_lab_id', 'inhibitor']).mean().reset_index()
    response_data = response_data.loc[
        ~(
                (response_data.aic > 12) &
                (response_data.deviance > 2)
        )
    ].copy()
    new_auc = response_data[['proteomic_lab_id', 'inhibitor', 'auc']].copy()
    new_auc.rename({'proteomic_lab_id': 'sample_id'}, axis=1, inplace=True)
    if not os.path.exists(f_name):
        new_auc.to_csv(f_name)
    return new_auc


def load_mutations():
    """
    Loads WES data.

    Processes mutational status into two levels. First one is at the gene level,
    second one gene with amino acid level.


    Returns
    -------

    """
    f_name = 'data/wes.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name)
    df = load_table(wes_id)
    mapper = {
        'symbol': 'gene_symbol',
        'labId': 'sample_id',
    }
    df.rename(mapper, axis=1, inplace=True)
    df['exp_value'] = 1
    wes_gene_level = pd.pivot_table(
        df,
        index='sample_id',
        values='exp_value',
        columns='gene_symbol',
        fill_value=0
    )
    wes_gene_level = wes_gene_level.melt(ignore_index=False).reset_index()
    wes_gene_level['label'] = wes_gene_level.gene_symbol + '_mut'
    wes_gene_level['exp_value'] = wes_gene_level['value']
    wes_gene_level['source'] = 'wes'

    wes_aa_level = df.copy()
    wes_aa_level['label'] = wes_aa_level['hgvsp'].str.split(':p.').str.get(1)
    wes_aa_level['label'] = wes_aa_level['gene_symbol'] + '_' + wes_aa_level['label']

    wes_aa_level = pd.pivot_table(
        wes_aa_level,
        index='sample_id',
        values='exp_value',
        columns='label',
        fill_value=0
    )
    wes_aa_level = wes_aa_level.melt(ignore_index=False).reset_index()
    wes_aa_level['gene_symbol'] = wes_aa_level.label.str.split('_').str.get(0)
    wes_aa_level['source'] = 'wes_protein_level'
    wes_aa_level['exp_value'] = wes_aa_level['value']
    merged = pd.concat([wes_gene_level, wes_aa_level])
    if not os.path.exists(f_name):
        merged.to_csv(f_name)
    return merged


class AMLData(object):
    def __init__(self):
        self.drug_names = None
        self._auc_table = None
        self.proteomics = prep_proteomics()
        self.phospho = prep_phosph()
        self.rna = prep_rnaseq()
        self.functional = load_drug_response()
        self.wes = load_mutations()
        self.flat_data = pd.concat(
            [self.phospho, self.proteomics, self.rna, self.wes]
        )

        self.meta = add_cluster_plus_meta()
        self.meta = self.meta.join(load_cluster_pred())
        self.flt3 = self.meta.loc[self.meta['FLT3-ITDcalls']].index.values
        self.non_flt3 = self.meta.loc[~self.meta['FLT3-ITDcalls']].index.values
        self.all_data = self.convert_to_matrix(self.flat_data)
        self.feature_names = list(self.all_data.columns.values)
        self.feature_names.remove('sample_id')

        # format for magine.ExperimentalData class
        d = self.flat_data.rename({'gene_symbol': 'identifier'}, axis=1)
        d['species_type'] = 'gene'
        self.exp_data = ExperimentalData(d)

    def add_meta(self, pivoted_table):
        return self.meta.join(pivoted_table)

    @property
    def auc_table(self):
        if self._auc_table is None:
            self._auc_table = pd.pivot_table(
                self.functional,
                index='sample_id', columns='inhibitor', values='auc'
            )
            self.drug_names = list(self._auc_table.columns.unique())
            self._auc_table = self._auc_table.join(self.meta, on='sample_id')
        return self._auc_table

    @auc_table.setter
    def auc_table(self, new_val):
        self._auc_table = new_val

    def subset(self, source, with_meta=False):
        if isinstance(source, str):
            source = [source]
        subset = self.flat_data.loc[self.flat_data.source.isin(source)]
        if with_meta:
            return self.add_meta(
                self.convert_to_matrix(subset).set_index('sample_id'))
        else:
            return self.convert_to_matrix(subset)

    def subset_flat(self, source):
        if isinstance(source, str):
            source = [source]
        return self.flat_data.loc[self.flat_data.source.isin(source)]

    def convert_to_matrix(self, flat_dataframe):
        df = pd.pivot_table(
            flat_dataframe,
            index='sample_id',
            values='exp_value',
            columns='label'
        ).reset_index()
        return df
        # until meta info is added to tables, commenting this out
        # return df.join(self.meta, on='sample_id').reset_index()

    def get_trainable_data(self, source, drug_name, new_format=False,
                           flt3_only=False, non_flt3_only=False,
                           cluster=None):
        # Filter experimental platform
        mol_data = self.subset(source)
        feature_names = list(mol_data.columns.values)
        if 'sample_id' in feature_names:
            feature_names.remove('sample_id')
        # merge with auc table to get row=patient, col=genes + drug_auc
        joined = mol_data.join(
            self.auc_table[drug_name],
            on='sample_id'
        ).set_index('sample_id')

        # remove rows without an AUC measurement
        df_subset = joined[~joined[drug_name].isna()].copy()

        # require 50% of the data be present for any given column
        df_subset.dropna(
            axis=0,
            how='any',
            thresh=df_subset.shape[1] * .50,
            inplace=True
        )
        # filter down if missing any measured cols
        # TODO Ask Camilo about the data filling
        n_features_remaining = df_subset.shape[0]

        df_subset.dropna(
            axis=1,
            how='any',
            thresh=n_features_remaining,
            inplace=True
        )
        if new_format:
            return SampleByClusterDataSet(
                self,
                df_subset,
                drug_name,
                flt3_only=flt3_only,
                non_flt3_only=non_flt3_only,
                cluster=cluster
            )
        return df_subset


class SampleByClusterDataSet(object):
    def __init__(self,
                 data,
                 df,
                 target_name,
                 flt3_only=False,
                 non_flt3_only=False,
                 cluster=None
                 ):
        self._df = df

        if flt3_only:
            index_vals = set(self._df.index.values).intersection(data.flt3)
            self.df = self._df.loc[index_vals].copy()
        elif non_flt3_only:
            index_vals = set(self._df.index.values).intersection(data.non_flt3)
            self.df = self._df.loc[index_vals].copy()

        if cluster:
            k5_cluster = data.meta['k=4'].copy()
            sub = k5_cluster[k5_cluster == cluster].index.values
            index_vals = set(self._df.index.values).intersection(sub)
            self._df = self._df.loc[index_vals]

        feat_names = list(set(self._df.columns.values))
        if target_name in feat_names:
            feat_names.remove(target_name)
        self.features = self._df[feat_names].copy()
        self.target = self._df[target_name].values * 1

    def train_test_split(self):

        x_train, x_test, y_train, y_test = train_test_split(
            self.features,
            self.target,
            test_size=0.2,
            shuffle=True,
            random_state=101,
        )

        return x_train, x_test, y_train, y_test

    def remove_features(self, feature_names):
        current_features = set(self.features.columns.values)
        self.features = self._df[
            current_features.difference(set(feature_names))]

    def require_features(self, feature_names):
        current_features = set(self.features.columns.values)
        fn = [i + '_rna' for i in feature_names]
        fn += [i + '_prot' for i in feature_names]
        self.features = self._df[current_features.intersection(set(fn))]

    def require_features_by_label(self, feature_names):
        current_features = set(self.features.columns.values)
        self.features = self._df[
            current_features.intersection(set(feature_names))]


def add_cluster_plus_meta():
    f_name = 'data/meta_labels.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name, index_col='sample_id')
    clusters = load_file(clusters_id)
    del clusters['Barcode.ID']
    clusters.reset_index(inplace=True)
    clusters.rename({'index': 'sample_id'}, axis=1, inplace=True)

    summary = load_excel(clinical_summary_id)
    summary['sample_id'] = summary['labId']
    summary.set_index('sample_id', inplace=True)
    summ_cols = [
        'FLT3-ITDcalls',
        'NPM1calls',
        'cumulativeChemo',
        'overallSurvival'
    ]
    summary = summary[summ_cols].copy()
    rename = {
        'positive': True,
        'negative': False,
        'y': True,
        'n': False,
        np.nan: False
    }
    summary.replace(rename, inplace=True)
    merged = summary.join(clusters.set_index('sample_id'))
    if not os.path.exists(f_name):
        merged.to_csv(f_name)
    return merged


def load_cluster_pred():
    f_name = 'data/cluster_pred.csv'
    f_name = os.path.join(os.path.dirname(__file__), f_name)
    if os.path.exists(f_name):
        return pd.read_csv(f_name, index_col=0)
    cluster_pred = load_file('syn30030154')
    cluster_pred['sample_id'] = cluster_pred['Barcode.ID']
    del cluster_pred['Barcode.ID']
    cluster_pred['Cluster'] = cluster_pred['Cluster'].str.split(' ').str.get(
        1).astype(int)
    cluster_pred.set_index('sample_id', inplace=True)
    if not os.path.exists(f_name):
        cluster_pred.to_csv(f_name)
    return cluster_pred


if __name__ == '__main__':
    d = AMLData()


