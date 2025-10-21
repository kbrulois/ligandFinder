




import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
import pandas as pd
import polars as pl
import numpy as np
from scipy.io import mmread
from scipy import sparse
import scipy.sparse
import json
import shap


res = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno.csv")
res["known_pair"] = res["known_pair"].map({"known": 1, "unknown": 0})


dat = pl.read_csv("/scratch/groups/ebutcher/deorphan/analysis/con_rep.csv", infer_schema_length=res.shape[0])
dat = dat.to_pandas()
# feat_names = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/feat_names.csv")
# dat = mmread("/scratch/groups/ebutcher/deorphan/analysis/xg_dat.mtx")
# dat = dat.tocsr()
# dat = dat.todense()
# dat = pd.DataFrame(dat)
# dat.columns = feat_names.columns


splits = [f"data_split_{i}" for i in range(1, 31)]

with open("/oak/stanford/groups/ebutcher/kevin/interaction_cons.json", "r") as f:
    int_cons = json.load(f)
    

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


df = pd.DataFrame({'data_set': ['train', 'not_train', 'test']})

for x in splits:
  x_train = dat.iloc[res[res[x] == "train"].index,:]
  y_train = np.array(res.loc[res[x] == "train", "known_pair"].to_numpy())
  x_test = dat.iloc[res[res[x] == "test"].index,:]
  y_test = np.array(res.loc[res[x] == "test", "known_pair"].to_numpy())
  x_not_train = dat.iloc[res[res[x] != "train"].index,:]
  y_not_train = np.array(res.loc[res[x] != "train", "known_pair"].to_numpy())
  dtrain = xgb.DMatrix(x_train, label=y_train)
  dnot_train = xgb.DMatrix(x_not_train, label=y_not_train)
  dtest = xgb.DMatrix(x_test, label=y_test)
  evallist = [(dtrain, 'train'), (dtest, 'test')]
  print(f"{x}")
  num_round = 370
  bst = xgb.train(param, dtrain, num_round, evallist, num_boost_round=400, verbose_eval=2)
  pred_train = bst.predict(dtrain)
  pred_not_train = bst.predict(dnot_train)
  pred_test = bst.predict(dtest)
  a = roc_auc_score(y_train, pred_train)
  b = roc_auc_score(y_not_train, pred_not_train)
  c = roc_auc_score(y_test, pred_test)
  df[x] = [a,b,c]

df.to_csv('/oak/stanford/groups/ebutcher/kevin/train_eval30.csv')



x = 'data_split_18'

x_train = dat.iloc[res[res[x] == "train"].index,:]
y_train = np.array(res.loc[res[x] == "train", "known_pair"].to_numpy())
x_test = dat.iloc[res[res[x] == "test"].index,:]
y_test = np.array(res.loc[res[x] == "test", "known_pair"].to_numpy())
x_not_train = dat.iloc[res[res[x] != "train"].index,:]
y_not_train = np.array(res.loc[res[x] != "train", "known_pair"].to_numpy())
dtrain = xgb.DMatrix(x_train, label=y_train)
dnot_train = xgb.DMatrix(x_not_train, label=y_not_train)
dtest = xgb.DMatrix(x_test, label=y_test)
evallist = [(dtrain, 'train'), (dtest, 'test')]
print(f"{x}")
num_round = 370
bst = xgb.train(param, dtrain, num_round, evallist, num_boost_round=400, verbose_eval=2)
pred_train = bst.predict(dtrain)
pred_not_train = bst.predict(dnot_train)
pred_test = bst.predict(dtest)
a = roc_auc_score(y_train, pred_train)
b = roc_auc_score(y_not_train, pred_not_train)
c = roc_auc_score(y_test, pred_test)


y_all = np.array(res["known_pair"].to_numpy())
dall = xgb.DMatrix(dat, label = y_all)

x_traintest = dat.iloc[res[res[x].isin(["train", "test"])].index,:]


explainer = shap.TreeExplainer(bst, x_train)
shap_values = explainer(x_traintest)

np.savetxt("/scratch/groups/ebutcher/deorphan/analysis/shap_tree_traintest.csv", np.squeeze(shap_values.values).T, delimiter=",")





col_index = res.columns.get_loc('code') + 1
res.insert(col_index, 'contact_score_all', np.nan)

res['contact_score_all'] = pred_all

res.loc[res["data_split"] != "train", 'contact_score_sep'] = pred_not_train

pred_train = bst_test.predict(dtrain)

res.loc[res["data_split"] == "train", 'contact_score_sep'] = pred_train


res.to_csv("/oak/stanford/groups/ebutcher/kevin/bm_contact_score_latest3.csv")











