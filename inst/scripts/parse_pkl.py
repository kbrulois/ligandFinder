
import pickle as pk
import numpy as np
import os


file_path = "/Users/kbrulois/peptide_alg/ACKR2_and_CXCL14_68W_czmm/result_model_1_multimer_v3_pred_0.pkl"

file_path = "/oak/stanford/groups/ebutcher/catherine/alphapulldown/test/test_pulldown_output/Q14240_and_P60842/result_model_1_multimer_v3_pred_0.pkl"

file_path = "/oak/stanford/groups/ebutcher/catherine/alphapulldown/pulldown_output/040824/BDKRB2_human_full_and_CXCL14-61-F/result_model_1_multimer_v3_pred_0.pkl"




with open(file_path, 'rb') as file:  
    data = pk.load(file)
    
dir = os.path.splitext(file_path)[0] + "_parsed_pkl"

os.makedirs(dir)

for key in data.keys():
	dat_type = type(data[key])
	if(dat_type == dict):
		for key2 in data[key].keys():
			dat_type2 = type(data[key][key2])
			print(key, key2, dat_type2)
			np.save(dir + "/" + key + "." + key2 + ".npy", np.asarray(data[key][key2]))
	else:
		print(key, dat_type)
		if(key == 'seqs'):
			with open(dir + "/" + key + ".txt", "w") as file:
				for item in data[key]:
					file.write(f"{item}\n")
		else:
			np.save(dir + "/" + key + ".npy", data[key])
