import lightgbm as lgb
import numpy as np
import pandas as pd
import sklearn.linear_model as LM
from scipy.stats import pearsonr
from sklearn import metrics
from sklearn import svm
from sklearn.model_selection import KFold

from pybeataml.load_data import AMLData

# point of access for all data
data = AMLData()

en_model = LM.ElasticNet(
    random_state=0,
    max_iter=100000,
    fit_intercept=True,
    l1_ratio=.7,
    alpha=.9,
)

lasso_model = LM.Lasso(
    alpha=0.1,
    max_iter=100000,
    fit_intercept=True,
)

svm_model = svm.SVR(C=1.0, epsilon=0.2, kernel='linear')
svc_model = svm.SVC(C=1.0, kernel='linear')


def run_sklearn(x_train, y_train, x_test, y_test, model, model_name, binarize=False):
    # don't think we need scaling? Can add in pretty quickly if we do
    #     pipe = pipeline.Pipeline([
    # #         ('scaler', preprocessing.StandardScaler()),
    #         ('model', model)
    #     ])
    if binarize:
        y_train, y_test = convert_to_binary(y_train, y_test)

    model.fit(x_train, y_train)
    preds = model.predict(x_test)
    error, r2, pearson = score_all(y_test, preds)
    coef = np.array(model.coef_).flatten()
    feature_names = model.feature_names_in_[coef > 0]
    auc = np.nan
    if binarize:
        auc = metrics.average_precision_score(y_test, preds)
    return {
        'pearsonr': pearsonr(y_test, preds)[0],
        'mse': error,
        'r2': r2,
        'model': model_name,
        'feature_names': feature_names,
        'auc': auc
    }


def convert_to_binary(y_train, y_test):
    """ Binarize AUC """
    y_train_c = np.copy(y_train)
    y_test_c = np.copy(y_test)
    y_train_c[y_train_c < 100] = 1
    y_train_c[y_train_c > 100] = 0
    y_test_c[y_test_c < 100] = 1
    y_test_c[y_test_c > 100] = 0
    return y_train_c, y_test_c


def run_gbt(x_train, y_train, x_test, y_test, feature_names, binarize=False):
    param = dict(
        device_type='cpu',
        boosting='gbdt',
        nthread=8,
        objective='regression',
        metric='rmse',
        lambda_l1=1,
        lambda_l2=1,
        learning_rate=.01,
        tree_learner='serial',
        max_bin=63,
        num_leaves=10,
        max_depth=10,
        feature_fraction=.5,
        min_data_in_leaf=1,
        min_gain_to_split=1,
        verbose=-1

    )
    model_name = 'gbt'

    if binarize:
        param['objective'] = 'binary'
        param['metric'] = 'auc'
        y_train, y_test = convert_to_binary(y_train, y_test)
        model_name = 'gbt_binary'

    train_data = lgb.Dataset(x_train, label=y_train, feature_name=feature_names)
    validation_data = lgb.Dataset(x_test, label=y_test, feature_name=feature_names)
    num_round = 1000
    bst = lgb.train(
        param,
        train_data,
        num_round,
        valid_sets=validation_data,
        callbacks=[lgb.early_stopping(stopping_rounds=100, verbose=0)]
    )
    table = bst.feature_importance()
    feats = pd.Series(table, index=feature_names)
    selected_feat = feats[feats > 0].index.values

    preds = bst.predict(x_test, num_iteration=bst.best_iteration)
    error, r2, pearson = score_all(y_test, preds)
    auc = np.nan
    if binarize:
        auc = metrics.average_precision_score(y_test, preds)
    return {
        'mse': error,
        'r2': r2,
        'pearsonr': pearsonr(y_test, preds)[0],
        'model': model_name,
        'feature_names': selected_feat,
        'auc': auc
    }


def score_all(y_test, preds):
    error = np.sqrt(metrics.mean_squared_error(y_test, preds))
    r2 = metrics.r2_score(y_test, preds)
    pearson = pearsonr(y_test, preds)[0]
    #     print(f"RMSE: {error:0.3f} | R^2 {r2:0.3f} | R {pearson:0.3f}")
    return error, r2, pearson


def run_model(d_sets, drug_name):
    df_subset = data.get_trainable_data(d_sets, drug_name)
    cols = list(set(df_subset.columns.values))
    cols.remove(drug_name)

    features = df_subset[cols].copy()
    target = df_subset[drug_name].values

    n_features_before = features.shape[1]
    # features = features.loc[:, features.mean() > 0]
    # features = features.loc[:, features.std() > 0]
    n_features = features.shape[1]
    print(f"Using {n_features} out of {n_features_before}"
          f" ({n_features_before - n_features} removed)")
    if features.shape[0] < 100:
        return pd.DataFrame()
    #     features = pd.DataFrame(features, columns=feature_names)
    feature_names = list(set(features.columns.values))

    all_results = []
    kf = KFold(n_splits=5, shuffle=True, random_state=101)
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
        gbt_binary_results = run_gbt(
            feature_names=feature_names, binarize=True, **args
        )

        enet_results = run_sklearn(
            model=en_model, model_name='EN', **args
        )

        lasso_results = run_sklearn(
            model=lasso_model, model_name='LASSO', **args
        )
        svm_results = run_sklearn(
            model=svm_model, model_name='SVM', **args
        )
        svc_results = run_sklearn(
            model=svm_model, model_name='SVC', binarize=True, **args
        )

        results = pd.DataFrame(
            [gbt_results, gbt_binary_results, enet_results, lasso_results, svm_results, svc_results]
        )
        results['k'] = n
        all_results.append(results)
        print(results[['model', 'mse', 'r2', 'pearsonr', 'auc']])

    if not isinstance(d_sets, str):
        out_name = '_'.join(sorted(d_sets))
    else:
        out_name = d_sets
    all_results = pd.concat(all_results)
    all_results['drug'] = drug_name
    all_results['data_type'] = out_name
    all_results.feature_names = all_results.feature_names.str.join('|')
    return all_results[['model', 'k', 'mse', 'r2', 'pearsonr', 'auc', 'drug', 'data_type', 'feature_names']]


if __name__ == '__main__':

    counts = data.auc_table[data.auc_table < 100].count().sort_values()
    drugs = counts[counts > 20].index.values
    # only run single drugs for now
    drug_solo = [i for i in drugs if ' - ' not in i]
    data_sources = [
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

    models = []
    for i in drug_solo:
        print(f"Working on {i}")
        for j in data_sources:
            models.append(run_model(j, i))

    df = pd.concat(
        models,
    )

    df.to_csv("ultimate_output_all_data.csv")
