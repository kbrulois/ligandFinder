





import pandas as pd
import numpy as np

import xgboost as xgb
from xgboost import XGBClassifier

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score

from scipy.io import mmread
from scipy import sparse
import scipy.sparse

import json
import shap


from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.regularizers import l2



dat = mmread("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_new.mtx")
dat = dat.tocsr()
dat = dat.todense()

#dat = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/con_rep.csv", engine="pyarrow")

feat_names = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/feat_names.csv")
feat_names = feat_names.columns.tolist()

res = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno.csv")
res["known_pair"] = res["known_pair"].map({"known": 1, "unknown": 0})

with open("/oak/stanford/groups/ebutcher/kevin/interaction_cons.json", "r") as f:
    int_cons = json.load(f)



splits = [f"data_split_{i}" for i in range(1, 21)]

df2 = pd.DataFrame(columns=['data_split', 'auc_train', 'auc_test', 'auc_not_train', 'auc_all', 'prec_train', 'prec_test', 'prec_not_train', 'prec_all'])

for x in splits:
  x_train = dat[res[res[x] == "train"].index,:]
  y_train = np.array(res.loc[res[x] == "train", "known_pair"].to_numpy())
  x_test = dat[res[res[x] == "test"].index,:]
  y_test = np.array(res.loc[res[x] == "test", "known_pair"].to_numpy())
  x_not_train = dat[res[res[x] != "train"].index,:]
  y_not_train = np.array(res.loc[res[x] != "train", "known_pair"].to_numpy())
  y_all = np.array(res["known_pair"].to_numpy())
  dtrain = xgb.DMatrix(x_train, label=y_train, feature_names = feat_names)
  dtest = xgb.DMatrix(x_test, label=y_test, feature_names = feat_names)
  dnot_train = xgb.DMatrix(x_not_train, label=y_not_train, feature_names = feat_names)
  dall = xgb.DMatrix(dat, label = y_all, feature_names = feat_names)
  evallist = [(dtrain, 'test'), (dtest, 'train')]
  evallist = [(dtrain, 'train'), (dtest, 'test')]
  param = {
  'max_depth': 6, 
  'eta': 0.08, 
  'tree_method': 'exact', #default 'auto'
  'objective': 'binary:logistic', 
  'subsample': 0.5, #default 1
  'colsample_bytree': 1, #default 1
  'min_child_weight': 0,
  'gamma': 0,#default 0
  'lambda': 0,#default 1
  'alpha': 1,#default 1
  'nthread': 16,
  'eval_metric': 'auc'}
  param['interaction_constraints'] = int_cons
  print(f"using {x}")
  num_round = 200
  bst = xgb.train(param, dtrain, num_round, evallist, num_boost_round=400, verbose_eval=1)
  pred_train = bst.predict(dtrain)
  pred_not_train = bst.predict(dnot_train)
  pred_test = bst.predict(dtest)
  pred_all = bst.predict(dall)
  roc_auc_score(y_train, pred_train)
  roc_auc_score(y_not_train, pred_not_train)
  roc_auc_score(y_test, pred_test)
  roc_auc_score(y_all, pred_all)
  new_row = pd.DataFrame([[x,
  roc_auc_score(y_train, pred_train),
  roc_auc_score(y_test, pred_test),
  roc_auc_score(y_not_train, pred_not_train),
  roc_auc_score(y_all, pred_all),
  average_precision_score(y_train, pred_train),
  average_precision_score(y_test, pred_test),
  average_precision_score(y_not_train, pred_not_train),
  average_precision_score(y_all, pred_all)]], 
  columns=df.columns)
  df2 = pd.concat([df, new_row], ignore_index=True)




