import lightgbm as lgb
import numpy as np
import pandas as pd
import sklearn.linear_model as LM
import itertools as it
from scipy.stats import pearsonr, spearmanr
from sklearn import metrics
from sklearn.model_selection import RepeatedKFold

from pybeataml.load_data import AMLData
from pybeataml.load_data_from_synpase import load_file
from datetime import date

# point of access for all data
data = AMLData()

en_model = LM.ElasticNet(
    random_state=0,
    max_iter=100000,
    fit_intercept=True,
    l1_ratio=.7,
    alpha=.9,
)

def run_sklearn(x_train, y_train, x_test, y_test, model, model_name):
    model.fit(x_train, y_train)
    preds = model.predict(x_test)
    coef = np.array(model.coef_).flatten()
    selected_feat = model.feature_names_in_[coef > 0]
    error, r2, pearson, spearman, pr, sr = score_all(y_test, preds)

    return {
        'model': model_name,
        'feature_names': sorted(selected_feat),
        'n_feat': len(selected_feat),
        'rmse': error,
        'r2': r2,
        'pearson': pearson,
        'spearman': spearman,
        'pr': pr,
        'sr': sr,
    }


good_param = dict(
    device_type='cpu',
    boosting_type='gbdt',
    num_threads=8,
    n_jobs=None,
    objective='regression',
    metric='rmse',
    lambda_l1=100,
    lambda_l2=10,
    reg_alpha=None,
    reg_lambda=None,
    learning_rate=.1,

    max_bin=128,
    num_leaves=21,  # 21
    max_depth=-1,  # -1

    zero_as_missing=True,

    feature_fraction=.8,  # .8

    bagging_freq=1,
    bagging_fraction=.75,
    subsample_freq=None,
    subsample=None,

    min_data_in_leaf=1,
    min_child_samples=None,
    colsample_bytree=None,
    min_split_gain=None,
    n_estimators=10000,
    verbose=-1,
    deterministic=True,
    random_state=10,

)
good_param = dict(
    device_type='cpu',
    boosting_type='gbdt',
    num_threads=8,
    n_jobs=None,
    objective='regression',
    metric='rmse',
    lambda_l1=1000,
    lambda_l2=1,
    reg_alpha=None,
    reg_lambda=None,
    learning_rate=.1,
    tree_learner='serial',
    max_bin=128,
    num_leaves=5,
    max_depth=-1,

    feature_fraction=1,  # .8

    bagging_freq=1,
    bagging_fraction=.8,
    subsample=None,
    subsample_freq=None,

    min_child_weight=0.2,
    min_data_in_leaf=2,
    min_child_samples=None,
    min_gain_to_split=None,
    colsample_bytree=None,
    min_split_gain=None,
    n_estimators=10000,
    verbose=-1,
    deterministic=True,
    random_state=10,

)


def run_gbt(x_train, y_train, x_test, y_test, feature_names):
    eval_metric = 'rmse'
    model_name = 'gbt'

    lgb_model = lgb.LGBMRegressor(**good_param)
    lgb_model.fit(
        x_train, y_train,
        callbacks=[lgb.early_stopping(
            500, first_metric_only=True, verbose=False)
        ],
        eval_metric=eval_metric,
        eval_set=[(x_train, y_train), (x_test, y_test)],
    )

    feats = pd.Series(lgb_model.feature_importances_, index=feature_names)
    selected_feat = feats[feats > 0].index.values

    preds = lgb_model.predict(x_test, num_iteration=lgb_model.best_iteration_)
    error, r2, pearson, spearman, pr, sr = score_all(y_test, preds)

    return {
        'model': model_name,
        'feature_names': sorted(selected_feat),
        'n_feat': len(selected_feat),
        'rmse': error,
        'r2': r2,
        'pearson': pearson,
        'spearman': spearman,
        'pr': pr,
        'sr': sr,
    }


def score_all(y_test, preds):
    error = np.sqrt(metrics.mean_squared_error(y_test, preds))
    r2 = metrics.r2_score(y_test, preds)
    pearson, pr = pearsonr(y_test, preds)
    spearman, sr = spearmanr(y_test, preds)
    return error, r2, pearson, spearman, pr, sr


def run_model(d_sets, drug_name):
    df_subset = data.get_trainable_data(d_sets, drug_name, new_format=True)
    features = df_subset.features
    target = df_subset.target
    feature_names = list(set(features.columns.values))

    all_results = []
    # kf = KFold(n_splits=10, shuffle=True, random_state=101)
    kf = RepeatedKFold(n_splits=5, n_repeats=5, random_state=101)

    for n, (train_index, test_index) in enumerate(kf.split(features)):
        x_train, x_test = features.iloc[train_index], features.iloc[test_index]
        y_train, y_test = target[train_index], target[test_index]

        args = dict(
            x_train=x_train,
            y_train=y_train,
            x_test=x_test,
            y_test=y_test
        )

        gbt_results = run_gbt(feature_names=feature_names, **args)

        enet_results = run_sklearn(
            model=en_model, model_name='EN', **args
        )

        # svm_results = run_sklearn(
        #     model=svm_model, model_name='SVM', **args
        # )

        results = pd.DataFrame(
            [enet_results, gbt_results]
        )

        # results = pd.DataFrame([gbt_results])
        results['k'] = n
        all_results.append(results)
        # print(results[['model', 'rmse', 'r2', 'pearson', 'auc']])

    if not isinstance(d_sets, str):
        out_name = '_'.join(sorted(d_sets))
    else:
        out_name = d_sets
    all_results = pd.concat(all_results)
    all_results['drug_name'] = drug_name
    all_results['data_type'] = out_name
    cols = ['model', 'n_feat', 'rmse', 'r2',
            'pearson', 'spearman', ]

    print('\t', d_sets, drug_name)
    with pd.option_context("display.precision", 2):
        print(all_results[cols])
    all_results.feature_names = all_results.feature_names.str.join('|')
    return all_results


if __name__ == '__main__':
    table = data.auc_table[data.drug_names].copy()

    # at least 100 samples
    drug_counts = table.describe().T['count']
    high_occ_drugs = drug_counts[drug_counts > 100].index.values

    counts = table[table < 100].count()
    responsive_drugs = counts[counts > 10].index.values
    # only run single drugs for now
    drug_solo = [i for i in responsive_drugs if ' - ' not in i]
    drugs_to_focus = [
        'Gilteritinib',
        'Quizartinib (AC220)',
        'Trametinib (GSK1120212)',
        'Sorafenib',
        'Panobinostat',
        'Venetoclax',
    ]
    good_drugs = set(drug_solo).intersection(high_occ_drugs)
    print(len(good_drugs))

    # adding FLT3 inhibitors back in.
    for i in drugs_to_focus:
        if i not in good_drugs:
            good_drugs.add(i)

    # generate all possible combinations of input data
    all_sources = ['wes', 'rna_seq', 'proteomics', 'phospho', 'metabolomics', 'lipidomics']
    data_sources = []
    for l in range(len(all_sources) + 1):
        for subset in it.combinations(all_sources, l):
            data_sources = [data_sources, subset]
    old_data_sources = [
        'rna_seq',
        'proteomics',
        'phospho',
        'wes',
        ['proteomics', 'phospho'],
        ['proteomics', 'rna_seq', ],
        ['proteomics', 'wes'],
        ['phospho', 'rna_seq', ],
        ['phospho', 'wes'],
        ['rna_seq', 'wes'],
        ['phospho', 'rna_seq', 'wes'],
        ['proteomics', 'rna_seq', 'wes'],
        ['proteomics', 'phospho', 'wes'],
        ['proteomics', 'phospho', 'rna_seq', 'wes'],
    ]
    data_sources = list(set(data_sources) - set(old_data_sources))

    models = []
    for i in reversed(sorted(good_drugs)):
        print(f"Working on {i}")
        for j in data_sources:
            models.append(run_model(j, i))

    df = pd.concat(
        models,
    )
    old_models = load_file('syn52299998')
    df_final = pd.concat(df, old_models)

    f_name = "regression_all_models_all_data_combos_cv_5v5" + date.today() + ".csv"
    df_final.to_csv(f_name)