import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import table

from magine.data.base import BaseData
from magine.data.tools import log2_normalize_df
from magine.plotting import volcano_plots as v_plot
from magine.plotting.species_plotting import plot_dataframe, plot_species

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

    def plot_pie_sig_ratio(self, save_name=None, ax=None, fig=None,
                           figsize=None):
        """

        Parameters
        ----------
        save_name : str
        ax : matplotlib.axes, optional
        fig : matplotlib.figure
        figsize : tuple
            Size of figure

        Returns
        -------

        """
        x = len(self.id_list)
        y = len(self.sig.id_list)
        total = x + y
        if fig is None and ax is None:
            if figsize is None:
                figsize = (3, 3)
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        wedges, texts, autotexts = ax.pie(
            [x, y],
            explode=(0.05, 0.05),
            textprops={'fontsize': 16},
            autopct=lambda p: '{:.0f}'.format(p * total / 100),
            shadow=True,
            startangle=140
        )

        plt.setp(autotexts, size=20)
        plt.axis('equal')
        if save_name is not None:
            plt.savefig('{}.png'.format(save_name), dpi=300,
                        bbox_inches='tight')
        return fig

    def volcano_plot(self, save_name=None, out_dir=None, sig_column=False,
                     p_value=0.1, fold_change_cutoff=1.5, x_range=None,
                     y_range=None):
        """ Create a volcano plot of data


        Parameters
        ----------
        save_name: str
            name to save figure
        out_dir: str, directory
            Location to save figure
        sig_column: bool, optional
            If to use significant flags of data
        p_value: float, optional
            Criteria for significant
        fold_change_cutoff: float, optional
            Criteria for significant
        y_range: array_like
            upper and lower bounds of plot in y direction
        x_range: array_like
            upper and lower bounds of plot in x direction

        Returns
        -------
        matplotlib.Figure

        """
        fig = v_plot.volcano_plot(self, save_name=save_name, out_dir=out_dir,
                                  sig_column=sig_column, p_value=p_value,
                                  fold_change_cutoff=fold_change_cutoff,
                                  x_range=x_range, y_range=y_range)
        return fig

    def volcano_by_sample(self, save_name=None, p_value=0.1,
                          out_dir=None, fold_change_cutoff=1.5, y_range=None,
                          x_range=None, sig_column=False):
        """
        Creates a figure of subplots of provided experimental method

        Parameters
        ----------
        save_name: str
            name to save figure
        out_dir: str, directory
            Location to save figure
        sig_column: bool, optional
            If to use significant flags of data
        p_value: float, optional
            Criteria for significant
        fold_change_cutoff: float, optional
            Criteria for significant
        y_range: array_like
            upper and lower bounds of plot in y direction
        x_range: array_like
            upper and lower bounds of plot in x direction

        Returns
        -------

        """

        data = self.copy()
        n_sample = np.sort(data[sample_id].unique())

        if len(n_sample) > 8:
            n_cols = 3
        else:
            n_cols = 2
        n_rows = int(np.rint(np.rint(len(n_sample) / float(n_cols))))
        if n_cols * n_rows < len(n_sample):
            if n_cols >= n_rows:
                n_rows += 1
            else:
                n_cols += 1

        fig = plt.figure(figsize=(4 * n_rows, 3 * n_cols))
        for n, i in enumerate(n_sample):
            sample = data[data[sample_id] == i].copy()

            sample = sample.dropna(subset=[p_val])
            sample = sample[np.isfinite(sample[exp_value])]
            sample = sample.dropna(subset=[exp_value])
            sec_0, sec_1, sec_2 = v_plot.create_mask(sample, sig_column,
                                                     p_value,
                                                     fold_change_cutoff)
            ax = fig.add_subplot(n_rows, n_cols, n + 1)
            ax.set_title(i)
            v_plot.add_volcano_plot(ax, sec_0, sec_1, sec_2)
            if not sig_column:
                fc = np.log2(fold_change_cutoff)
                log_p_val = -1 * np.log10(p_value)
                ax.axvline(x=fc, linestyle='--')
                ax.axvline(x=-1 * fc, linestyle='--')
                ax.axhline(y=log_p_val, linestyle='--')
            if y_range is not None:
                ax.set_ylim(y_range[0], y_range[1])
            if x_range is not None:
                ax.set_xlim(x_range[0], x_range[1])
        fig.tight_layout()
        if save_name is not None:
            v_plot.save_plot(fig, save_name=save_name, out_dir=out_dir)
        return fig

    def plot_species(self, species_list=None, subset_index=None,
                     save_name=None, out_dir=None, title=None,
                     plot_type='plotly', image_format='png'):
        """
        Create scatter plot of species list

        Parameters
        ----------
        species_list : list
            list of compounds
        subset_index : list
            Column to filter based on species_list
        save_name : str
            Name of html output file
        out_dir : str
            Location to place plots
        title : str
            Title for HTML page
        plot_type : str
            Type of plot outputs, can be "plotly" or "matplotlib"
        image_format : str
            pdf or png, only used if plot_type="matplotlib"

        Returns
        -------
        matplotlib.Figure or plotly.Figure
        """
        df = self.copy()
        if species_list is not None:
            if subset_index is None:
                subset_index = self._index
            df = df.subset(species_list, index=subset_index)
        return plot_species(
            df, save_name=save_name, out_dir=out_dir, title=title,
            plot_type=plot_type, image_format=image_format
        )

    def plot_all(self, html_file_name, out_dir='out',
                 plot_type='plotly', run_parallel=False):
        """
        Creates a plot of all metabolites

        Parameters
        ----------
        html_file_name : str
            filename to save html of all plots
        out_dir: str, path
            Directory that will contain all proteins
        plot_type : str
            plotly or matplotlib output
        run_parallel : bool
            Create the plots in parallel
        Returns
        -------

        """

        plot_dataframe(self, html_filename=html_file_name,
                       out_dir=out_dir, plot_type=plot_type,
                       run_parallel=run_parallel)

    def plot_histogram(self, save_name=None, y_range=None, out_dir=None):
        """
        Plots a histogram of data

        Parameters
        ----------
        save_name: str
            Name of figure
        out_dir: str, path
            Path to location to save figure
        y_range: array_like
            range of data


        Returns
        -------

        """
        data = self.copy()
        data = data.dropna(subset=[p_val])
        data = data[np.isfinite(data[exp_value])]
        data = data.dropna(subset=[exp_value])

        tmp = np.array(log2_normalize_df(data, exp_value)[exp_value])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(tmp, 50, color='gray')
        if y_range is not None:
            plt.xlim(y_range[0], y_range[1])

        ax.set_yscale('log', basey=10)
        ax.set_xlabel('log$_2$ Fold Change', fontsize=16)
        ax.set_ylabel('Count', fontsize=16)
        fig.tight_layout()
        if save_name is not None:
            v_plot.save_plot(fig, save_name, out_dir)
        return fig


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
