"""
First pass at data centered class.
Ideally a single class instance can handle loading of the entire dataset.
TODO
    - Add cache mechanism for quick storage/loading rather than going to synap
    - Think of a way to get current feature cols and drug cols from
    a subset data, add functionality
"""
import pandas as pd

from pybeataml.data import ExperimentalData
from pybeataml.load_data_from_synpase import load_table, load_file

# current synapse ids, check with Camilo to see if these are the final (
# I know there are some other corrected/v2/uncorrected in the R code)
global_id = 'syn25808020'
phospho_id = 'syn26477193'  # syn25808662
rnaseq_id = 'syn26545877'
drug_response_id = 'syn25830473'
meta_file_id = 'syn26534982'
wes_id = 'syn26428827'


def prep_rnaseq():
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
    return subset


def prep_phosph():
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

    return phosph_subset


def prep_proteomics():
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
    return proteomics


def load_drug_response():
    response_data = load_table(drug_response_id)
    new_auc = response_data[['lab_id', 'inhibitor', 'auc']].copy()
    new_auc.rename({'lab_id': 'sample_id'}, axis=1, inplace=True)
    new_auc.auc = new_auc.auc.astype(float)
    return new_auc


def load_meta_data():
    meta = load_file(meta_file_id)


def load_mutations():
    """
    Loads WES data.

    Processes mutational status into two levels. First one is at the gene level,
    second one gene with amino acid level.


    Returns
    -------

    """
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
    return pd.concat([wes_gene_level, wes_aa_level])


class AMLData(object):
    def __init__(self):
        self._drug_names = None
        self._auc_table = None
        self.proteomics = prep_proteomics()
        self.phospho = prep_phosph()
        self.rna = prep_rnaseq()
        self.functional = load_drug_response()
        self.wes = load_mutations()
        self.flat_data = pd.concat(
            [self.phospho, self.proteomics, self.rna, self.wes]
        )

        # Until i find the tables, commenting this
        meta_info_cols = ['sample_id', 'InitialAMLDiagnosis',
                          'PostChemotherapy', 'FLT3.ITD']
        # meta = self.flat_data[meta_info_cols].drop_duplicates()
        # meta.set_index('sample_id', inplace=True)
        # self.meta = meta * 1
        self.all_data = self.convert_to_matrix(self.flat_data)
        self.feature_names = list(self.all_data.columns.values)
        self.feature_names.remove('sample_id')

        # format for magine.ExperimentalData class
        d = self.flat_data.rename({'gene_symbol': 'identifier'}, axis=1)
        d['species_type'] = 'gene'
        self.exp_data = ExperimentalData(d)

    @property
    def drug_names(self):
        if self._drug_names is None:
            self._drug_names = list(set(self.functional['inhibitor'].unique()))
        return self._drug_names

    @property
    def auc_table(self):
        if self._auc_table is None:
            self._auc_table = pd.pivot_table(
                self.functional,
                index='sample_id', columns='inhibitor', values='auc'
            )
        return self._auc_table

    def subset(self, source):
        if isinstance(source, str):
            source = [source]
        subset = self.flat_data.loc[self.flat_data.source.isin(source)]
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

    def get_trainable_data(self, source, drug_name):
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
        # df_subset = joined.loc[:, feature_names + [drug_name]]
        # remove rows without a AUC measurement
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
        return df_subset


if __name__ == '__main__':
    data = AMLData()
