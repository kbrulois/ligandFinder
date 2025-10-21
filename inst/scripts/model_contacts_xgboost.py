




col_index = res.columns.get_loc('code') + 1
res.insert(col_index, 'contact_score_all', np.nan)

res['contact_score_all'] = pred_all

res.loc[res["data_split"] != "train", 'contact_score_sep'] = pred_not_train

pred_train = bst_test.predict(dtrain)

res.loc[res["data_split"] == "train", 'contact_score_sep'] = pred_train


res.to_csv("/oak/stanford/groups/ebutcher/kevin/bm_contact_score_latest4.csv")



dat2 = mmread("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_new.mtx")
dat2 = dat2.tocsr()
dat2 = dat2.todense()

#dat = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/con_rep.csv", engine="pyarrow")

feat_names = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/feat_names.csv")
feat_names = feat_names.columns.tolist()

res2 = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno_new.csv")
res2["known_pair"] = res2["known_pair"].map({"known": 1, "unknown": 0})


dnew = xgb.DMatrix(dat2, label=np.array(res2.loc[:, "known_pair"].to_numpy()), feature_names = feat_names)

pred_new = bst.predict(dnew)

res2['contact_score_all'] = pred_new


res2.to_csv("/oak/stanford/groups/ebutcher/kevin/cxcl14_jh_w_contact_score.csv")







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

int_cons2 = [dtrain.feature_names for y in int_cons]

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






bst_test = xgb.train(param, dtest, num_round, evallist, num_boost_round=400)




pred_all = bst.predict(dall)



roc_auc_score(y_all, pred_all)


explainer = shap.TreeExplainer(bst, x_train)
shap_values = explainer(x_test)

np.savetxt("/oak/stanford/groups/ebutcher/kevin/shap_tree_no_int.csv", np.squeeze(shap_values.values).T, delimiter=",")



col_index = res.columns.get_loc('code') + 1
res.insert(col_index, 'contact_score_all', np.nan)

res['contact_score_all'] = pred_all

res.loc[res["data_split"] != "train", 'contact_score_sep'] = pred_not_train

pred_train = bst_test.predict(dtrain)

res.loc[res["data_split"] == "train", 'contact_score_sep'] = pred_train


res.to_csv("/oak/stanford/groups/ebutcher/kevin/bm_contact_score_latest3.csv")




x_train = dat[res[res["data_split"] == "train"].index,:]
y_train = np.array(res.loc[res["data_split"] == "train", "known_pair"].to_numpy())

x_not_train = dat[res[res["data_split"] != "train"].index,:]
y_not_train = np.array(res.loc[res["data_split"] != "train", "known_pair"].to_numpy())


x_val = dat[res[res["data_split"] == "validate"].index,:]
y_val = np.array(res.loc[res["data_split"] == "validate", "known_pair"].to_numpy())













import os

import xgboost as xgb

dtrain = xgb.DMatrix()
dtest = xgb.DMatrix(
    os.path.join(CURRENT_DIR, "../data/agaricus.txt.test?format=libsvm")
)

param = [
    ("max_depth", 2),
    ("objective", "binary:logistic"),
    ("eval_metric", "logloss"),
    ("eval_metric", "error"),
]

num_round = 2
watchlist = [(dtest, "eval"), (dtrain, "train")]

evals_result = {}
bst = xgb.train(param, dtrain, num_round, watchlist, evals_result=evals_result)

print("Access logloss metric directly from evals_result:")
print(evals_result["eval"]["logloss"])

print("")
print("Access metrics through a loop:")
for e_name, e_mtrs in evals_result.items():
    print("- {}".format(e_name))
    for e_mtr_name, e_mtr_vals in e_mtrs.items():
        print("   - {}".format(e_mtr_name))
        print("      - {}".format(e_mtr_vals))

print("")
print("Access complete dictionary:")
print(evals_result)






















dat = pd.read_csv("/scratch/groups/ebutcher/deorphan/analysis/8M_bm_cons.csv", engine="pyarrow")
dat["known_pair"] = dat["known_pair"].map({"known": 1, "unknown": 0})

y_all = np.array(dat["known_pair"].to_numpy())

to_change=['gpcr_family', 'lig1_end', 'protein_segment', 'BW', 'AA_pair_type']

for item in to_change:
    dat[item]=dat[item].astype('category')
    
dat.drop('known_pair', axis=1, inplace=True)
  
x='data_split_1'

dtrain = xgb.DMatrix(dat[dat[x] == "train"], label=y_all[dat[x] == "train"], enable_categorical=True)
dtest = xgb.DMatrix(x_test, label=y_test, feature_names = feat_names)
dnot_train = xgb.DMatrix(x_not_train, label=y_not_train, feature_names = feat_names)
dall = xgb.DMatrix(dat, label = y_all, feature_names = feat_names)


evallist = [(dtrain, 'test')]
evallist = [(dtrain, 'train')]

param = {
'max_depth': 6, 
'eta': 0.08, 
'tree_method': 'auto', #default 'auto'
'objective': 'binary:logistic', 
'subsample': 0.5, #default 1
'colsample_bytree': 1, #default 1
'min_child_weight': 0,
'gamma': 0,#default 0
'lambda': 0,#default 1
'alpha': 1,#default 1
'nthread': 16,
'eval_metric': 'auc'}





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




