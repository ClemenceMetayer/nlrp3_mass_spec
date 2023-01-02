import pandas as pd
import pickle as pkl
import numpy as np
import os
import warnings
from multiprocessing import Pool
import statsmodels as sm
import cosinor_simple

warnings.simplefilter(action='ignore', category=FutureWarning) # cosinorPy uses soon-to-be deprecated method pd.append...
warnings.simplefilter(action='ignore', category=RuntimeWarning)

# Open the data dictionary 
with open('data/data_dict.dat', 'rb') as f:
    data = pkl.load(f)
    
data['name'] = "complete_dictionary"

# conds : control or nlrp3
# cts : circadian times
# nreps : number of replicates
# nrna : number of proteins
conds, cts, nreps, nrna = ["ctrl", "nlrp3"], data["CTs"], 4, len(data["ctrl"])
nrna = 34700
# Code name of the proteins
names = list(data["ctrl"].keys())

# Formatting the underlying dictionaries :
ctrl = data["ctrl"]
nlrp3 = data["nlrp3"]

lenght = 4957

dico1_ctrl = dict(list(ctrl.items())[0:lenght])
dico2_ctrl = dict(list(ctrl.items())[lenght:2*lenght])
dico3_ctrl = dict(list(ctrl.items())[2*lenght:3*lenght])
dico4_ctrl = dict(list(ctrl.items())[3*lenght:4*lenght])
dico5_ctrl = dict(list(ctrl.items())[4*lenght:5*lenght])
dico6_ctrl = dict(list(ctrl.items())[5*lenght:6*lenght])
dico7_ctrl = dict(list(ctrl.items())[6*lenght:7*lenght+1])

dico1_nlrp3 = dict(list(nlrp3.items())[0:lenght])
dico2_nlrp3 = dict(list(nlrp3.items())[lenght:2*lenght])
dico3_nlrp3 = dict(list(nlrp3.items())[2*lenght:3*lenght])
dico4_nlrp3 = dict(list(nlrp3.items())[3*lenght:4*lenght])
dico5_nlrp3 = dict(list(nlrp3.items())[4*lenght:5*lenght])
dico6_nlrp3 = dict(list(nlrp3.items())[5*lenght:6*lenght])
dico7_nlrp3 = dict(list(nlrp3.items())[6*lenght:7*lenght+1])


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
    # Cosinor results
    df_results = pd.DataFrame(columns = ['test', 'period','n_components','p','q','RSS','amplitude','acrophase', 'mesor','n_points'], dtype=float)
    range_period = [20.5, 20.6, 20.7, 20.8, 20.9, 21, 21.1, 21.2, 21.3, 21.4, 21.5, 23.5, 23.6, 23.7, 23.8, 23.9, 24, 24.1, 24.2, 24.3, 24.4, 24.5]
    for period in range_period :        
        df_results_single_per = cosinor_simple.fit_group(df, period= period)
        df_results = df_results.append(df_results_single_per)
    
    # Calculate the q_value for multiple testing
    pval_threshold = .05
    names = np.unique(df_results.test)   
    for test in names :
        df_test = df_results[df_results.test == test]
        df_qval = df_test.p
        df_results.loc[df_results.test == test, "q"] = sm.stats.multitest.multipletests(df_qval,alpha = pval_threshold, method ='fdr_bh',
        is_sorted=False, returnsorted=False)[1]

    # df_periodic_pval results
    df_periodic_pval = cosinor_simple.get_periodic_pval(df_results)
   
    # df_periodic_qval results
    df_periodic_qval = cosinor_simple.get_periodic_qval(df_results)
     
   
    # Export to CSV
    name_dico = dictionary['name']
    output_path_1 = 'data/Cosinor/cosinor_data'
    output_path_2 = 'data/Cosinor/cosinor_results'
    output_path_4 = 'data/Cosinor/cosinor_periodic_pval'
    output_path_5 = 'data/Cosinor/cosinor_periodic_qval'
    filepath = os.path.join(output_path_1, 'cosinor_data_'+str(name_dico)+'.csv')
    df.to_csv(filepath,index=False)
    filepath = os.path.join(output_path_2, 'cosinor_results_'+str(name_dico)+'.csv')
    df_results.to_csv(filepath, index=False)
    filepath = os.path.join(output_path_4, 'cosinor_periodic_pval_'+str(name_dico)+'.csv')
    df_periodic_pval.to_csv(filepath, index=False)
    filepath = os.path.join(output_path_5, 'cosinor_periodic_qval_'+str(name_dico)+'.csv')
    df_periodic_qval.to_csv(filepath, index=False)
 
if __name__ == "__main__":
    pool = Pool(7)
    results = pool.map(cosinor_results,dico.values())
    print(results)

# df_results concatenation
df_results_serie1 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie1.csv")
df_results_serie2 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie2.csv")
df_results_serie3 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie3.csv")
df_results_serie4 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie4.csv")
df_results_serie5 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie5.csv")
df_results_serie6 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie6.csv")
df_results_serie7 = pd.read_csv("data/Cosinor/cosinor_results/cosinor_results_serie7.csv")

df_results = pd.concat([df_results_serie1, df_results_serie2, df_results_serie3, df_results_serie4, df_results_serie5, 
                        df_results_serie6, df_results_serie7], axis=0)
output_path = 'data/Cosinor/cosinor_results'
filepath = os.path.join(output_path, 'cosinor_results_concat_serie.csv')
df_results.to_csv(filepath,index=False)


# df_periodic_pval concatenation
df_periodic_pval_serie1 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie1.csv")
df_periodic_pval_serie2 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie2.csv")
df_periodic_pval_serie3 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie3.csv")
df_periodic_pval_serie4 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie4.csv")
df_periodic_pval_serie5 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie5.csv")
df_periodic_pval_serie6 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie6.csv")
df_periodic_pval_serie7 = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_serie7.csv")

df_periodic_pval = pd.concat([df_periodic_pval_serie1, df_periodic_pval_serie2, df_periodic_pval_serie3, df_periodic_pval_serie4, 
                              df_periodic_pval_serie5, df_periodic_pval_serie6, df_periodic_pval_serie7], axis=0)
output_path = 'data/Cosinor/cosinor_periodic_pval'
filepath = os.path.join(output_path, 'cosinor_periodic_pval_concat_serie.csv')
df_periodic_pval.to_csv(filepath,index=False)


# df_periodic_qval concatenation
df_periodic_qval_serie1 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie1.csv")
df_periodic_qval_serie2 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie2.csv")
df_periodic_qval_serie3 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie3.csv")
df_periodic_qval_serie4 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie4.csv")
df_periodic_qval_serie5 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie5.csv")
df_periodic_qval_serie6 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie6.csv")
df_periodic_qval_serie7 = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_serie7.csv")

df_periodic_qval = pd.concat([df_periodic_qval_serie1, df_periodic_qval_serie2, df_periodic_qval_serie3, df_periodic_qval_serie4, 
                              df_periodic_qval_serie5, df_periodic_qval_serie6, df_periodic_qval_serie7], axis=0)
output_path = 'data/Cosinor/cosinor_periodic_qval'
filepath = os.path.join(output_path, 'cosinor_periodic_qval_concat_serie.csv')
df_periodic_qval.to_csv(filepath,index=False)

