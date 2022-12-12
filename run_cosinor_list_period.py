import pandas as pd
import pickle as pkl
import numpy as np
import os
import warnings
from multiprocessing import Pool
import statsmodels as sm
from CosinorPy import cosinor_float_period

warnings.simplefilter(action='ignore', category=FutureWarning) # cosinorPy uses soon-to-be deprecated method pd.append...
warnings.simplefilter(action='ignore', category=RuntimeWarning)

from CosinorPy import cosinor
# Open the data dictionary 
with open('data/data_dict.dat', 'rb') as f:
    data = pkl.load(f)
    
data['name'] = "complete_dictionary"

# conds : control or nlrp3
# cts : circadian times
# nreps : number of replicates
# nprots : number of proteins
conds, cts, nreps, nprots = ["ctrl", "nlrp3"], data["CTs"], 4, len(data["ctrl"])
nprots = 10570
# Code name of the proteins
names = list(data["ctrl"].keys())

# Formatting the underlying dictionaries :
ctrl = data["ctrl"]
nlrp3 = data["nlrp3"]

pas = 7
lenght = int(nprots/pas)

dico1_ctrl = dict(list(ctrl.items())[0:lenght])
dico2_ctrl = dict(list(ctrl.items())[lenght:2*lenght])
dico3_ctrl = dict(list(ctrl.items())[2*lenght:3*lenght])
dico4_ctrl = dict(list(ctrl.items())[3*lenght:4*lenght])
dico5_ctrl = dict(list(ctrl.items())[4*lenght:5*lenght])
dico6_ctrl = dict(list(ctrl.items())[5*lenght:6*lenght])
dico7_ctrl = dict(list(ctrl.items())[6*lenght:7*lenght])

dico1_nlrp3 = dict(list(nlrp3.items())[0:lenght])
dico2_nlrp3 = dict(list(nlrp3.items())[lenght:2*lenght])
dico3_nlrp3 = dict(list(nlrp3.items())[2*lenght:3*lenght])
dico4_nlrp3 = dict(list(nlrp3.items())[3*lenght:4*lenght])
dico5_nlrp3 = dict(list(nlrp3.items())[4*lenght:5*lenght])
dico6_nlrp3 = dict(list(nlrp3.items())[5*lenght:6*lenght])
dico7_nlrp3 = dict(list(nlrp3.items())[6*lenght:7*lenght])


dico1 = {"ctrl" : dico1_ctrl, "nlrp3" : dico1_nlrp3, 'name' : 'serie1'}
dico2 = {"ctrl" : dico2_ctrl, "nlrp3" : dico2_nlrp3, 'name' : 'serie2'}
dico3 = {"ctrl" : dico3_ctrl, "nlrp3" : dico3_nlrp3, 'name' : 'serie3'}
dico4 = {"ctrl" : dico4_ctrl, "nlrp3" : dico4_nlrp3, 'name' : 'serie4'}
dico5 = {"ctrl" : dico5_ctrl, "nlrp3" : dico5_nlrp3, 'name' : 'serie5'}
dico6 = {"ctrl" : dico6_ctrl, "nlrp3" : dico6_nlrp3, 'name' : 'serie6'}
dico7 = {"ctrl" : dico7_ctrl, "nlrp3" : dico7_nlrp3, 'name' : 'serie7'}
dico = {'dico1' : dico1, 'dico2' : dico2, 'dico3' : dico3, 'dico4' : dico4, 'dico5' : dico5,
        'dico6' : dico6, 'dico7' : dico7}


# Runames_protsing Cosinor for a range of period and storing the results for each 
def cosinor_results(dictionary) :
    # formatting data for CosinorPy, df is a data frame with 3 columns (name, circadian times with repetition, value)
    df = pd.DataFrame()
    names_prots = list(dictionary["ctrl"].keys())
    # Sum of the circadian times considerinf the replicates
    x = sum([[float(4 * i)] * nreps for i in range(1, len(cts) + 1)], [])
    # Iteration on the number of proteins
    for k in range(len(dictionary["ctrl"])):
        # Iteration on [1,"ctrl"] and [2,"nlrp3"] 
        for idx_cond, cond in enumerate(conds):
            # y store the protein's data for each circadian time and each replicates
            y = np.array(dictionary[cond][names_prots[k]]["pts"]).reshape(-1,1).squeeze()
            
            notnan = np.argwhere(~np.isnan(y)) # cosinorPy does not handle NaNs: they are removed.
            y, xnew = y[notnan].squeeze(), np.array(x)[notnan].squeeze()
            
            df = pd.concat((df, pd.DataFrame({"test":[f'{names_prots[k]}/{cond}'] * len(notnan),
                                              "x": xnew,
                                              "y": y})))
    # this may take a while
    df_results = cosinor.fit_group(df, n_components = 1, 
                                   period=[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],  # period can also be a list, e.g. [22, 23, 24, 25, 26]
                                   plot=False)

    # Solve the pvalue problem for mulitple testing
    pval_threshold = .05
    names = np.unique(df_results.test)   
    for test in names :
        df_test = df_results[df_results.test == test]
        df_qval = df_test.p
        df_results.loc[df_results.test == test, "q"] = sm.stats.multitest.multipletests(df_qval,alpha = pval_threshold, method ='fdr_bh',
        is_sorted=False, returnsorted=False)[1]

    
    # Pick the best period 
    df_best_fits = cosinor_float_period.get_best_fits(df,df_results)    

    
    # Export to CSV
    name_dico = dictionary['name']
    output_path_1 = 'data/Cosinor/cosinor_data'
    output_path_2 = 'data/Cosinor/cosinor_results'
    output_path_3 = 'data/Cosinor/cosinor_best_models'
    filepath = os.path.join(output_path_1, 'cosinor_data_'+str(name_dico)+'.csv')
    df.to_csv(filepath,index=False)
    filepath = os.path.join(output_path_2, 'cosinor_results_'+str(name_dico)+'.csv')
    df_results.to_csv(filepath, index=False)
    filepath = os.path.join(output_path_3, 'cosinor_best_models_'+str(name_dico)+'.csv')
    df_best_fits.to_csv(filepath, index=False)
 
if __name__ == "__main__":
    pool = Pool(7)
    results = pool.map(cosinor_results,dico.values())
    print(results)

df_bm_serie1 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie1.csv")
df_bm_serie2 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie2.csv")
df_bm_serie3 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie3.csv")
df_bm_serie4 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie4.csv")
df_bm_serie5 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie5.csv")
df_bm_serie6 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie6.csv")
df_bm_serie7 = pd.read_csv("data/Cosinor/cosinor_best_models/cosinor_best_models_serie7.csv")

df_bm = pd.concat([df_bm_serie1, df_bm_serie2, df_bm_serie3, df_bm_serie4, df_bm_serie5, df_bm_serie6, df_bm_serie7], axis=0)
output_path = 'data/Cosinor/cosinor_best_models'
filepath = os.path.join(output_path, 'cosinor_best_models_concat_serie.csv')
df_bm.to_csv(filepath,index=False)
