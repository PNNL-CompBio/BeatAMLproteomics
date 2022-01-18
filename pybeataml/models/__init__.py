import lightgbm as lgb
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sklearn.linear_model as LM
from scipy.stats import pearsonr
from sklearn import pipeline
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

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


def run_sklearn(x_train, y_train, x_test, y_test, model, model_name):
    pipe = pipeline.Pipeline([
        #         ('scaler', preprocessing.StandardScaler()),
        ('model', model)
    ])

    pipe.fit(x_train, y_train)
    train_pred = pipe.predict(x_train)
    preds = pipe.predict(x_test)
    error, r2, pearson = score_all(y_test, preds)
    return {
        'test_prediction': preds,
        'train_prediction': train_pred,
        'pearsonr': pearsonr(y_test, preds)[0],
        'mse': error,
        'r2': r2,
        'model': model_name
    }


def run_gbt(x_train, y_train, x_test, y_test, feature_names):
    #     scaler = preprocessing.StandardScaler()
    #     x_train = scaler.fit_transform(x_train)
    #     x_test = scaler.fit_transform(x_test)
    y_train[y_train < 100] = 1
    y_train[y_train > 100] = 0
    y_test[y_test < 100] = 1
    y_test[y_test > 100] = 0
    train_data = lgb.Dataset(x_train, label=y_train,
                             feature_name=feature_names)
    validation_data = lgb.Dataset(x_test, label=y_test,
                                  feature_name=feature_names, )

    param = dict(
        device_type='gpu',
        boosting='gbdt',
        nthread=1,
        objective='binary',
        metric='auc',
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

    num_round = 1000
    bst = lgb.train(
        param,
        train_data,
        num_round,
        valid_sets=validation_data,
        callbacks=[lgb.early_stopping(stopping_rounds=100, verbose=1)]
    )

    train_pred = bst.predict(x_train, num_iteration=bst.best_iteration)
    preds = bst.predict(x_test, num_iteration=bst.best_iteration)
    error, r2, pearson = score_all(y_test, preds)

    return {
        'test_prediction': preds,
        'train_prediction': train_pred,
        'mse': error,
        'r2': r2,
        'pearsonr': pearsonr(y_test, preds)[0],
        'model': 'gbt'
    }


def score_all(y_test, preds):
    error = np.sqrt(mean_squared_error(y_test, preds))
    r2 = r2_score(y_test, preds)
    pearson = pearsonr(y_test, preds)[0]
    print(f"RMSE: {error:0.3f} | R^2 {r2:0.3f} | R {pearson:0.3f}")
    return error, r2, pearson


def run_model(data, drug_name):
    #     df_subset = data.get_trainable_data(['rna_seq'], drug_name)
    df_subset = data.get_trainable_data(['phospho', 'proteomics', 'rna_seq'],
                                        drug_name)

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

    feature_names = list(set(features.columns.values))

    x_train, x_test, y_train, y_test = train_test_split(
        features,
        target,
        test_size=0.2,
        shuffle=True,
        random_state=101,
    )
    gbt_results = run_gbt(x_train, y_train, x_test, y_test, feature_names)
    enet_results = run_sklearn(x_train, y_train, x_test, y_test, en_model,
                               'EN')
    lasso_results = run_sklearn(x_train, y_train, x_test, y_test, lasso_model,
                                'LASSO')

    avg = np.zeros(len(enet_results['test_prediction']))
    for i in [gbt_results, enet_results, lasso_results]:
        sns.regplot(y_test, i['test_prediction'], label=i['model'])
        avg += i['test_prediction']
    avg = avg / 3
    score_all(y_test, avg)
    sns.regplot(y_test, avg, label='Average')
    plt.xlabel("Actual")
    plt.ylabel("Predicted")
    plt.legend()
    plt.suptitle(f"Drug: {drug_name}")
