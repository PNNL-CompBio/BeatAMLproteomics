"""
Script to run models for each drug using cross validation.
Resulting output file can be viewed with plot_model_summaries.ipynb

"""
import pandas as pd
import sklearn.linear_model as LM
from sklearn import preprocessing, pipeline
from sklearn.model_selection import cross_validate

from pybeataml.load_data import AMLData

data = AMLData()


def run_k_fold(sources, drug_name, model):
    """

    Parameters
    ----------
    sources : list
        List of experimental sources to use
    drug_name : str
        Name of drug
    model : LinearModel
        Model for k-fold

    Returns
    -------

    """
    if not isinstance(sources, str):
        out_name = '_'.join(sorted(sources))
    else:
        out_name = sources

    df_subset = data.get_trainable_data(sources, drug_name)

    cols = list(set(df_subset.columns.values))
    cols.remove(drug_name)

    features = df_subset[cols].copy()
    target = df_subset[drug_name].values

    n_features_before = features.shape[1]

    features = features.loc[:, features.mean() > 0]
    features = features.loc[:, features.std() > 0]
    n_features = features.shape[1]
    print(f"Using {n_features} out of {n_features_before}"
          f" ({n_features_before - n_features} removed)")

    if features.shape[0] < 100:
        return []

    pipe = pipeline.Pipeline(
        [
            ('scaler', preprocessing.StandardScaler()),
            ('model', model)
        ]
    )
    scores = cross_validate(pipe, features, target, cv=5, n_jobs=5,
                            scoring=['neg_mean_squared_error', 'r2'])

    return [out_name, drug_name, scores['test_r2'],
            scores['test_neg_mean_squared_error']]


if __name__ == '__main__':
    save = True

    model = LM.ElasticNet(
        max_iter=10000,
        fit_intercept=True,
        l1_ratio=.5,
        alpha=.9
    )
    # model = svm.SVR(C=1.0, epsilon=0.2, kernel='linear')
    # model = LGBMRegressor()
    # Require 10% (20 samples) have an AUC lower than 100
    counts = data.auc_table[data.auc_table < 100].count().sort_values()
    drugs = counts[counts > 20].index.values
    # only run single drugs for now
    drug_solo = [i for i in drugs if ' - ' not in i]
    data_sources = [
        'rna_seq', 'proteomics', 'phospho', 'wes',  # each by themselves
        ['proteomics', 'phospho'],
        ['rna_seq', 'proteomics', 'phospho', 'wes'],
    ]

    models = []
    for i in drug_solo:
        print(f"Working on {i}")
        for j in data_sources:
            models.append(run_k_fold(j, i, model))

    df = pd.DataFrame(models,
                      columns=['data_type', 'drug', 'r2', 'mse'])
    df.dropna(inplace=True)
    df['r2'] = df['r2'].apply(lambda x: '|'.join([str(i) for i in x]))
    df['mse'] = df['mse'].apply(lambda x: '|'.join([str(i) for i in x]))
    if save:
        df.to_csv("all_models_performance_en.csv")
