"""
Script to run models for each drug using cross validation with lightGBM.
Still in progress.

# TODO need to refactor xgboost to lightgbm.


Resulting output file can be viewed with plot_model_summaries.ipynb

"""
import itertools

import lightgbm as lgb
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn import preprocessing
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split


def run_lgb(sources, drug_name):
    if not isinstance(sources, str):
        out_name = '_'.join(sorted(sources))
    else:
        out_name = sources
    df_subset = data.get_trainable_data(sources, drug_name)

    cols = list(set(df_subset.columns.values))
    cols.remove(drug_name)

    features = df_subset[cols].copy()
    target = df_subset[drug_name].values.reshape(-1, 1).ravel()

    n_features_before = features.shape[1]

    features = features.loc[:, features.mean() > 0]
    features = features.loc[:, features.std() > 0]
    feature_names = list(set(features.columns.values))
    n_features = features.shape[1]
    print(f"Using {n_features} out of {n_features_before}"
          f" ({n_features_before - n_features} removed)")
    if features.shape[0] < 100:
        return {}

    X_train, X_test, y_train, y_test = train_test_split(
        features,
        target,
        test_size=0.2,
        shuffle=True,
        random_state=101,
    )

    scaler = preprocessing.StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.fit_transform(X_test)

    train_data = lgb.Dataset(X_train, label=y_train,
                             feature_name=feature_names)
    validation_data = lgb.Dataset(X_test, label=y_test,
                                  feature_name=feature_names, )

    param = dict(
        device_type='gpu',
        boosting='gbdt',
        nthread=1,
        objective='regression',
        metric='rmse',
        #         lambda_l1=.5,
        #         lambda_l2=.5,
        learning_rate=.01,
        tree_learner='serial',
        max_bin=63,
        num_leaves=6,
        max_depth=6,
        feature_fraction=.5,
        min_data_in_leaf=1,
        min_gain_to_split=1,
        verbose=-1

    )

    num_round = 1000
    bst = lgb.train(
        param,
        train_data,
        num_round,
        valid_sets=validation_data,
        callbacks=[lgb.early_stopping(stopping_rounds=100)]
    )
    # bst.save_model('model.txt', num_iteration=bst.best_iteration)
    #     lgb.plot_importance(bst, figsize =(4, 8))
    #     plt.show()

    # t_preds = bst.predict(X_train, num_iteration=bst.best_iteration)
    preds = bst.predict(X_test, num_iteration=bst.best_iteration)
    error = np.sqrt(mean_squared_error(y_test, preds))
    r2 = r2_score(y_test, preds)
    print(f"MSE: {error:0.3f} | $R^2$ {r2}")
    return {
        'data_sets': out_name,
        'drug_name': drug_name,
        'mse': error,
        'r2': r2
    }


def run_kfold(sources, drug_name):
    df_subset = data.get_trainable_data(sources, drug_name)

    cols = list(set(df_subset.columns.values))
    cols.remove(drug_name)

    features = df_subset[cols].copy()
    target = df_subset[drug_name].values.reshape(-1, 1).ravel()

    n_features_before = features.shape[1]
    features = features.loc[:, features.mean() > 0]
    features = features.loc[:, features.std() > 0]
    feature_names = list(set(features.columns.values))
    n_features = features.shape[1]
    print(f"Using {n_features} out of {n_features_before}"
          f" ({n_features_before - n_features} removed)")

    features = features.values

    # way 1
    features_output = []
    output_tracker = []
    scaler = preprocessing.StandardScaler()
    kf = KFold(n_splits=5, shuffle=True, random_state=101)
    for n, (train_index, test_index) in enumerate(kf.split(features)):
        X_train, X_test = features[train_index], features[test_index]
        y_train, y_test = target[train_index], target[test_index]

        X_train = scaler.fit_transform(X_train)
        X_test = scaler.fit_transform(X_test)

        # organize data into xgb data matrix
        train = xgb.DMatrix(data=X_train, label=y_train)
        test = xgb.DMatrix(data=X_test, label=y_test)

        # add gene names as feature labels
        train.feature_names = feature_names
        test.feature_names = feature_names

        num_round = 100
        results = dict()
        model = xgb.train(
            params, train, num_round,
            verbose_eval=100,
            early_stopping_rounds=100,
            evals=[(train, 'train'), (test, 'valid')],
            evals_result=results,

        )

        feature_scores = model.get_fscore()
        s = pd.Series(list(feature_scores.values()), index=feature_scores)

        s = s.to_frame().reset_index()
        s.rename(columns={0: 'count'}, inplace=True)
        s.index.name = 'feature'
        s['nfold'] = n
        features_output.append(s)

        # trained
        t_preds = model.predict(train,
                                iteration_range=(0, model.best_iteration + 1))
        # predictions
        preds = model.predict(test)

        error = np.sqrt(mean_squared_error(y_test, preds))
        r2 = r2_score(y_test, preds)
        print(f"MSE: {error:0.3}")
        print(f"$R^2$ {r2}")
        output_tracker.append({'nfold': n, 'r2': r2, 'error': error})
    return pd.DataFrame(output_tracker), pd.concat(features_output,
                                                   ignore_index=True)


def run_input_sets(drug_name):
    output = []
    input_sets = [
        'rna_seq',
        'proteomics',
        'phospho',
        # ['rna_seq', 'proteomics'],
        # ['rna_seq', 'phospho'],
        ['proteomics', 'phospho'],
        ['rna_seq', 'proteomics', 'phospho'],
    ]
    for i in input_sets:
        print(f"\t {i}")
        result = run_lgb(i, drug_name)
        if 'model' in result:
            del result['model']
        output.append(result)
    return output


def run_for_all():
    all_output = []
    drugs = ['Venetoclax', 'Gilteritinib', 'Quizartinib (AC220)',
             'Trametinib (GSK1120212)', 'Sorafenib']
    counts = data.auc_table[data.auc_table < 100].count().sort_values()
    drugs = counts[counts > 20].index.values
    drug_solo = [i for i in drugs if ' - ' not in i]
    # for i in data.drug_names:
    for i in sorted(drug_solo):
        print(f"Working on {i}")
        x = run_input_sets(i)
        all_output.append(x)
    all_output = list(itertools.chain.from_iterable(all_output))

    return pd.DataFrame(all_output)


if __name__ == '__main__':
    from pybeataml.load_data import AMLData

    data = AMLData()
    df = run_for_all()
    df.to_csv('gbt_all_drugs_all_dtypes4.csv')
