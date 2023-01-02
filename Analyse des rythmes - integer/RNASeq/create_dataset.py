"""Create data dict from the raw MS data and format it for RAIN."""

import pandas as pd
import numpy as np
import pickle as pkl
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning) # cosinorPy uses soon-to-be deprecated method pd.append...
warnings.simplefilter(action='ignore', category=RuntimeWarning)


### Create data dictionnary from raw MS data ##################################

# For every protein and for every time point, we have 3 control replicates and 3 sinlrp3 replicates 
ms = pd.read_csv('data/F067_HUBER_raw_counts.txt', sep='\t')

# Keep only the non null rows
ms_full = pd.DataFrame()
for row in ms.iterrows():
    row = row[1]
    count = 0
    for idx in range(len(row)) :
        if row[idx] != 0 :
            count = count +1
    if count !=0 :
        ms_full = ms_full.append(row)
        

ms_index = ms_full.index # Line names
#ms.drop('SUM_INT.txt', axis=1, inplace=True) # Name of the first column

        
# Way of filling up the list (fill the list with the name of the columns wih some conditions)
NLRP3_cols = [col for col in ms_full.columns if 'NLRP3' in col] # Column names siNRLP3
CTRL_cols = [col for col in ms_full.columns if 'NLRP3' not in col] # Column names siCTRL
nrna = len(ms_full)
    

conds = ["ctrl", "nlrp3"]
# Dictionaries are non-ordered collections of objects. 
# Starting with 2 keys (ctrl and nlrp3) but no associated values
data = {c: {} for c in conds}
cts = [4, 8, 12, 16, 20, 24] #Circadian time
# Creation of a list with 6 list in it corresponding to the name of colummns by replicates at a certain circadian time.
ct_ctrl_columns = [[col for col in CTRL_cols if f'CT{c}' in col] for c in cts]
ct_nlrp3_columns = [[col for col in NLRP3_cols if f'CT{c}' in col] for c in cts]
# Both list combine into one 
ct_cols = [ct_ctrl_columns, ct_nlrp3_columns]


# The function enumerate() return the component of conds with a counter [(0, 'ctrl'), (1, 'nlrp3')]
# Iteration on (0, 'ctrl') or (1, 'nlrp3')
for idx_cond, cond in enumerate(conds):
    # Iteration on the number of proteins 
     for n in range(nrna):
          print(cond,n)
          # Extraction of the values for every circadian time on a list of the columns corresponding to ctrl/nlrp3 and one protein 
          li = [ms_full.loc[:, ct_cols[idx_cond][j]].iloc[n].values for j in range(len(cts))]  
          li = np.array(li,dtype=np.float)
          # Creation of a dictionnary prot with several keys (pts, mean_vector, std_vector) and the corresponding values
          # pts : list with the values for the n protein and the condition cond classify by circadian time
          # mean_vector : list with mean of values per circadian time
          # std_vector : list with std of values per circadian time
          prot = {"pts": li, 
                  "mean_vector": np.array([np.nanmean(list(l)) if not pd.isna(list(l)).all() else -1 for l in li]),  # all-nans-CTs will not be shown at plot time
                  "std_vector": np.array([np.nanstd(list(l)) if not pd.isna(list(l)).all() else 0 for l in li])} 
          # The data dictionnary is filled with pts, mean_vector, std_vector for all protein and for ctrl/nrlp3
          data[cond][ms_index[n]] = prot 
# New key in the data dictionnary which store the circadian times
data["CTs"] = cts
# New key in the data dictionnary which store the names of ms_full columns
data["cols"] = ms_full.columns


# Save in binary format a python object ('wb') 
with open('data/data_dict.dat', 'wb') as f:
    pkl.dump(data, f)

### Formatting data to run with rain ##########################################

nreps = 3
# Creation of ncols (6 circadian times * 3 replicates) and names which store all the names of protein
ncols, names = len(data["CTs"]) * nreps, list(data["ctrl"].keys())
# Creation of the cols list which store the names ZT_i_j+1 for every circadian time and every replicate
cols = ["Name"] + sum([[f"ZT_{i}_{j+1}" for j in range(nreps)] for i in data['CTs']], [])
# A panda's serie is a one-dimensional ndarray with axis labels (including time series).
series = []
series_ctrl = []
series_nlrp3 = []
# Iteration on the number of proteins
for k in range(nrna):
    print(k)
    se_ctrl = pd.Series(np.zeros(ncols), dtype=object)
    se_nlrp3 = pd.Series(np.zeros(ncols), dtype=object)
    idx = 0
    # Iteration on the circadian times
    for i in range(len(data["CTs"])):
        # Iteration on the number of replicates
        for j in range(nreps): 
            # se_ctrl/se_nlrp3 store the values per circadian times and replicates for one protein at a time 
            se_ctrl.iloc[idx] = float(data["ctrl"][names[k]]["pts"][i][j])
            se_nlrp3.iloc[idx] = float(data["nlrp3"][names[k]]["pts"][i][j])
            idx += 1
    # For each protein, "series" store succesively series with values for each circadian times and replicates first for the control group and second for the nlrp3 group
    series.append(se_ctrl)
    series.append(se_nlrp3)
    series_ctrl.append(se_ctrl)
    series_nlrp3.append(se_nlrp3)

# Creation of the data frame with circadian times and replicates as columns and proteins as rows.
df = pd.DataFrame(series)
df_ctrl = pd.DataFrame(series_ctrl)
df_nlrp3 = pd.DataFrame(series_nlrp3)
df.columns = cols[1:]
df_ctrl.columns = cols[1:]
df_nlrp3.columns = cols[1:]
df.to_csv("data/ms_data_rain.csv", sep="\t", index=False)
df_ctrl.to_csv("data/ms_data_rain_ctrl.csv",index=False)
df_nlrp3.to_csv("data/ms_data_rain_nlrp3.csv", index=False)

# Creation of the data frame with the row's names for the RAIN analysis
names_wcond = sum([[f"{n}_ctrl", f"{n}_nlrp3"] for n in names[:nrna]], [])
names_wcond = pd.DataFrame(names_wcond)
names_wcond.to_csv("data/names_wcond.csv", sep="\t", index=False)

names_wcond_ctrl = sum([[f"{n}_ctrl"] for n in names[:nrna]], [])
names_wcond_ctrl = pd.DataFrame(names_wcond_ctrl)
names_wcond_ctrl.to_csv("data/names_wcond_ctrl.csv", index=False)

names_wcond_nlrp3 = sum([[f"{n}_nlrp3"] for n in names[:nrna]], [])
names_wcond_nlrp3 = pd.DataFrame(names_wcond_nlrp3)
names_wcond_nlrp3.to_csv("data/names_wcond_nlrp3.csv", index=False)

