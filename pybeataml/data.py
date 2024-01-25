import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import table

flag = 'significant'

# pandas.set_option('display.max_colwidth', -1)
# column definitions

exp_value = 'exp_value'
flag = 'significant'
exp_method = 'source'
p_val = 'p_value'
rna = 'rna_seq'
gene = 'gene'
protein = 'protein'
metabolites = 'metabolites'
species_type = 'species_type'
sample_id = 'sample_id'
identifier = 'identifier'
label = 'label'
valid_cols = [exp_value, flag, p_val, species_type, sample_id]


class BaseData(pd.DataFrame):
    """
    This class derived from pd.DataFrame
    """
    _index = None

    def __init__(self, *args, **kwargs):
        super(BaseData, self).__init__(*args, **kwargs)

    @property
    def _constructor(self):
        return BaseData

    @property
    def sig(self):
        """ terms with significant flag """
        return self.loc[self[flag]].copy()

    def pivoter(self, convert_to_log=False, columns='sample_id',
                values='fold_change', index=None, fill_value=None, min_sig=0):
        """ Pivot data on provided axis.

        Parameters
        ----------
        convert_to_log : bool
            Convert values column to log2
        index : str
            Index for pivot table
        columns : str
            Columns to pivot
        values : str
            Values of pivot table
        fill_value : float, optional
            Fill pivot table nans with
        min_sig : int
            Required number of significant terms to keep in a row, default 0

        Returns
        -------

        """
        d_copy = self.copy()
        if index is None:
            index = self._index

        if convert_to_log:
            d_copy.log2_normalize_df(values, inplace=True)

        if min_sig:
            if not isinstance(min_sig, int):
                raise AssertionError()
            if 'significant' not in d_copy.columns:
                print('In order to filter based on minimum sig figs, '
                      'please add a "significant" column')

            d_copy.require_n_sig(index=index, columns=columns,
                                 n_sig=min_sig,
                                 inplace=True)
            if not d_copy.shape[0]:
                return pd.DataFrame()

        array = pd.pivot_table(d_copy, index=index, fill_value=fill_value,
                               columns=columns, values=values)
        if isinstance(values, list):
            return array
        if isinstance(columns, list):
            array.sort_values(
                by=sorted(tuple(map(tuple, d_copy[columns].values))),
                ascending=False, inplace=True
            )
        elif isinstance(columns, str):
            array.sort_values(by=sorted(d_copy[columns].unique()),
                              ascending=False, inplace=True)
        return array

    def require_n_sig(self, columns='sample_id', index=None,
                      n_sig=3, inplace=False,
                      verbose=False):
        """ Filter index to have at least "min_terms" significant species.

        Parameters
        ----------
        columns : str
            Columns to consider
        index : str, list
            The column with which to filter by counts
        n_sig : int
            Number of terms required to not be filtered
        inplace : bool
            Filter in place or return a copy of the filtered data
        verbose : bool

        Returns
        -------
        new_data : BaseData
        """
        if index is None:
            index = self._index
        # create safe copy of array
        new_data = self.copy()

        # get list of columns
        cols_to_check = list(new_data[columns].unique())
        # convert boolean to numeric (didn't used to need to do this, but for some reason
        # pandas changed?
        new_data[flag] = pd.to_numeric(new_data[flag])

        if flag not in new_data.columns:
            raise AssertionError('Requires significant column')
        # pivot
        sig = pd.pivot_table(new_data,
                             index=index,
                             fill_value=0,
                             values=flag,
                             columns=columns
                             )[cols_to_check]

        # convert everything that's not 0 to 1
        sig[sig > 0] = 1
        sig = sig[sig.T.sum() >= n_sig]
        if isinstance(index, list):
            keepers = {i[0] for i in sig.index.values}
            new_data = new_data[new_data[index[0]].isin(keepers)]
        elif isinstance(index, str):
            n_before = len(new_data[index].unique())
            keepers = {i for i in sig.index.values}
            new_data = new_data.loc[new_data[index].isin(keepers)]
            n_after = len(new_data[index].unique())
            if verbose:
                print("Number in index went from {} to {}"
                      "".format(n_before, n_after))
        else:
            print("Index is not a str or a list. What is it?")

        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def present_in_all_columns(self, columns='sample_id',
                               index=None, inplace=False):
        """ Require index to be present in all columns

        Parameters
        ----------
        columns : str
            Columns to consider
        index : str, list
            The column with which to filter by counts
        inplace : bool
            Filter in place or return a copy of the filtered data

        Returns
        -------
        new_data : BaseData
        """
        if index is None:
            index = self._index
        # create safe copy of array
        new_data = self.copy()
        n_before = len(new_data[index].unique())
        # get list of columns
        cols_to_check = list(new_data[columns].unique())

        if flag not in new_data.columns:
            raise AssertionError("Missing {} column in data".format(flag))

        # pivot
        pivoted_df = pd.pivot_table(new_data, index=index, fill_value=np.nan,
                                    values=flag, columns=columns
                                    )[cols_to_check]

        # sig = pivoted_df.loc[~np.any(np.isnan(pivoted_df.values), axis=1)]
        sig = pivoted_df.loc[~pivoted_df.isnull().T.any()]
        if isinstance(index, list):
            keepers = {i[0] for i in sig.index.values}
            new_data = new_data[new_data[index[0]].isin(keepers)]
        elif isinstance(index, str):
            keepers = {i for i in sig.index.values}
            new_data = new_data.loc[new_data[index].isin(keepers)]
        else:
            print("Index is not a str or a list. What is it?")
        n_after = len(new_data[index].unique())

        print("Number in index went from {} to {}".format(n_before, n_after))

        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data

    def log2_normalize_df(self, column='fold_change', inplace=False):
        """ Convert "fold_change" column to log2.

        Does so by taking log2 of all positive values and -log2 of all negative
        values.

        Parameters
        ----------
        column : str
            Column to convert
        inplace : bool
            Where to apply log2 in place or return new dataframe

        Returns
        -------

        """
        new_data = self.copy()
        greater = new_data[column] > 0
        less = new_data[column] < 0
        new_data.loc[greater, column] = np.log2(new_data[greater][column])
        new_data.loc[less, column] = -np.log2(-new_data[less][column])
        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data


class Sample(BaseData):
    """ Provides tools for subsets of data types

    """

    def __init__(self, *args, **kwargs):
        super(Sample, self).__init__(*args, **kwargs)
        # self.drop_duplicates(inplace=True)
        self._index = identifier
        self._identifier = identifier
        self._value_name = exp_value
        self._sample_id_name = sample_id
        self._label = label
        self._up = None
        self._down = None
        self._sig = None

    @property
    def _constructor(self):
        return Sample

    @property
    def exp_methods(self):
        """ List of sample_ids in data"""
        return sorted(set(self[exp_method].values))

    @property
    def sample_ids(self):
        """ List of sample_ids in data"""
        return sorted(set(self[sample_id].values))

    @property
    def up(self):
        """return up regulated species"""
        return self.loc[self[flag] & (self[self._value_name] > 0)]

    @property
    def down(self):
        """return down regulated species"""
        return self.loc[self[flag] & (self[self._value_name] < 0)]

    @property
    def id_list(self):
        """ Set of species identifiers """
        return set(self[self._identifier].values)

    @property
    def label_list(self):
        """ Set of species labels """
        return set(self[self._label].values)

    @property
    def up_by_sample(self):
        """List of up regulated species by sample"""
        return [self.loc[self[sample_id] == i].up.id_list
                for i in self.sample_ids]

    @property
    def down_by_sample(self):
        """List of down regulated species by sample"""
        return [self.loc[self[sample_id] == i].down.id_list
                for i in self.sample_ids]

    @property
    def by_sample(self):
        """List of significantly flagged species by sample"""
        return [self.loc[self[sample_id] == i].id_list
                for i in self.sample_ids]

    def subset(self, species=None, index='identifier', sample_ids=None,
               exp_methods=None):
        """

        Parameters
        ----------
        species : list, str
            List of species to create subset dataframe from
        index : str
            Index to filter based on provided 'species' list
        sample_ids : str, list
            List or string to filter sample
        exp_methods : str, list
            List or string to filter sample

        Returns
        -------
        magine.data.experimental_data.Species
        """
        df = self.copy()
        if isinstance(species, str):
            df = df.loc[df[index].str.contains(species)]
        elif isinstance(species, (list, tuple, set)):
            df = df.loc[df[index].isin(species)]
        if sample_ids is not None:
            if isinstance(species, str):
                df = df.loc[df[sample_id].str.contains(sample_ids)]
            else:
                df = df.loc[df[sample_id].isin(sample_ids)]
        if exp_methods is not None:
            if isinstance(species, str):
                df = df.loc[df[exp_method].str.contains(exp_methods)]
            else:
                df = df.loc[df[exp_method].isin(exp_methods)]
        return df




class ExperimentalData(object):
    """
    Manages all experimental data

    """

    def __init__(self, data_file):
        """

        Parameters
        ----------
        data_file : str, pandas.DataFrame
            Name of file, generally csv.
            If provided a str, the file will be read in as a pandas.DataFrame


        """
        if isinstance(data_file, pd.DataFrame):
            df = data_file.copy()
        else:
            df = pd.read_csv(data_file, parse_dates=False, low_memory=False)
        df.reset_index(drop=True, inplace=True)
        df.drop_duplicates(inplace=True)
        for i in valid_cols:
            if i not in df.dtypes:
                print("{} not in columns.".format(i))

        self.data = BaseData(df)
        self._index = 'identifier'
        self.__proteins = None
        self.__genes = None
        self.__species = None
        self.__rna = None
        self.__compounds = None
        for i in self.exp_methods:
            self.__setattr__(i, Sample(
                self.data.loc[self.data[exp_method] == i]))

    def __setattr__(self, name, value):
        super(ExperimentalData, self).__setattr__(name, value)

    def __getitem__(self, name):
        return super(ExperimentalData, self).__getattribute__(name)

    @property
    def genes(self):
        """ All data tagged with gene

        Includes protein and RNA.

        Returns
        -------

        """
        if self.__genes is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[species_type].isin([protein, gene])]
            self.__genes = Sample(tmp)
        return self.__genes

    @property
    def proteins(self):
        """ Protein level data

        Tagged with "gene" identifier that is not RNA

        Returns
        -------

        """
        if self.__proteins is None:
            tmp = self.data.copy()
            tmp = tmp.loc[(self.data[species_type].isin([protein, gene])) &
                          ~(tmp[exp_method] == rna)]
            self.__proteins = Sample(tmp)
        return self.__proteins

    @property
    def rna(self):
        """ RNA level data

        Tagged with "RNA"

        Returns
        -------

        """
        if self.__rna is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[exp_method] == rna]
            self.__rna = Sample(tmp)
        return self.__rna

    @property
    def compounds(self):
        """ Only compounds in data

        Returns
        -------
        Sample

        """
        if self.__compounds is None:
            tmp = self.data.copy()
            tmp = tmp.loc[tmp[species_type] == metabolites]
            self.__compounds = Sample(tmp)
        return self.__compounds

    @property
    def species(self):
        """ Returns data in Sample format

        Returns
        -------
        Sample

        """
        if self.__species is None:
            self.__species = Sample(self.data.copy())
        return self.__species

    @property
    def exp_methods(self):
        """ List of source columns """
        return list(self.data[exp_method].unique())

    @property
    def sample_ids(self):
        """ List of sample_ids """
        return sorted(list(self.data[sample_id].unique()))

    def subset(self, species, index='identifier'):
        """

        Parameters
        ----------
        species : list, str
            List of species to create subset dataframe from
        index : str
            Index to filter based on provided 'species' list

        Returns
        -------
        magine.data.experimental_data.Species
        """
        df = self.species.copy()
        if isinstance(species, str):
            df = df.loc[df[index].str.contains(species)]
        else:
            df = df.loc[df[index].isin(species)]
        return df

    def get_measured_by_datatype(self):
        """
        Returns dict of species per data type

        Returns
        -------
        dict

        """
        return get_measured_by_datatype(self)

    def create_summary_table(self, sig=False, index=identifier, save_name=None,
                             plot=False):
        """
        Creates a summary table of data.


        Parameters
        ----------
        sig: bool
            Flag to summarize significant species only
        save_name: str
            Name to save csv and .tex file
        index: str
           Index for counts
        plot: bool
            If you want to create a plot of the table
        write_latex: bool
            Create latex file of table


        Returns
        -------
        pandas.DataFrame

        """
        return create_table_of_data(self, sig=sig, index=index,
                                    save_name=save_name, plot=plot)

    def volcano_analysis(self, out_dir, use_sig_flag=True,
                         p_value=0.1, fold_change_cutoff=1.5):
        """
        Creates a volcano plot for each experimental method

        Parameters
        ----------
        out_dir: str, path
            Path to where the output figures will be saved
        use_sig_flag: bool
            Use significant flag of data
        p_value: float, optional
            p value criteria for significant
            Will not be used if use_sig_flag
        fold_change_cutoff: float, optional
            fold change criteria for significant
            Will not be used if use_sig_flag

        Returns
        -------

        """
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for i in self.exp_methods:
            self[i].volcano_plot(
                i, out_dir=out_dir, sig_column=use_sig_flag,
                p_value=p_value, fold_change_cutoff=fold_change_cutoff
            )


def get_measured_by_datatype(data):
    """ Get unique list of species for each 'source' label in data.

    Parameters
    ----------
    data : ExperimentalData

    Returns
    -------
    measured, sig_measured : dict, dict
        Dictionaries where keys are 'source' and values are sets of ids.

    """

    measured = dict()
    sig_measured = dict()
    for i in data.exp_methods:
        sig_measured[i] = set(data[i].sig.id_list)
        measured[i] = set(data[i].id_list)
    return measured, sig_measured


def create_table_of_data(data, sig=False, index='identifier', save_name=None,
                         plot=False):
    """
    Creates a summary table of data.


    Parameters
    ----------
    data : ExperimentalData
    sig: bool
        Flag to summarize significant species only
    save_name: None, str
        Name to save csv and .tex file
    index: str
        Index to create counts
    plot: bool
        If you want to create a plot of the table


    Returns
    -------
    pandas.DataFrame

    """

    if sig:
        data_copy = data.species.sig.copy()
    else:
        data_copy = data.species.copy()

    count_table = data_copy.pivot_table(values=index, index=exp_method,
                                        columns=sample_id, fill_value=np.nan,
                                        aggfunc=lambda x: x.dropna().nunique())

    # This just makes sure things are printed as ints, not floats
    for i in count_table.columns:
        count_table[i] = count_table[i].fillna(-1).astype(int).replace(-1, '-')
    unique_col = {}
    for i in data.exp_methods:
        if sig:
            unique_col[i] = len(set(data[i].sig[index].values))
        else:
            unique_col[i] = len(set(data[i][index].values))
    count_table['Total Unique Across'] = pd.Series(unique_col,
                                                   index=count_table.index)
    if plot:
        ax = plt.subplot(111, frame_on=False)

        table(ax, count_table, loc='center')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.tight_layout()
        if save_name is not None:
            plt.savefig('{}.png'.format(save_name), dpi=300,
                        bbox_inches='tight')

    if save_name is not None:
        count_table.to_csv('{}.csv'.format(save_name))
    return count_table
