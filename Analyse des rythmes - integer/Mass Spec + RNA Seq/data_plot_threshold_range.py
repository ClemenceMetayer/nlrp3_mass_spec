#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 11:22:39 2022

@author: clemencemetayer
"""

# For each period tested, we store the number of rhythmic proteins according to the pvalues
from multiprocessing import Pool
import numpy as np 
import pandas as pd
import os

df_periodic_pval = pd.read_csv("data/Cosinor/cosinor_periodic_pval/cosinor_periodic_pval_concat_serie.csv")
df_periodic_qval = pd.read_csv("data/Cosinor/cosinor_periodic_qval/cosinor_periodic_qval_concat_serie.csv")

Range_period = range(10,29)


def rhythm_prots_pval(period) :
    rhythm_prots_period = pd.DataFrame(columns=['pval','nb_prots_ctrl','nb_prots_nlrp3','nb_prots_both'])
    for pval in np.arange(0,0.06,0.001) :
        cosinor_rhythms = [[], []]
        new_row = {"pval" : [pval]}
        for i in range(len(df_periodic_pval)):
            [name, cond] = df_periodic_pval['test'].iloc[i].split('/')
            if df_periodic_pval['period'].iloc[i] == period and cond == "ctrl" and df_periodic_pval['p'].iloc[i] < pval :
                cosinor_rhythms[0].append(name)
            if df_periodic_pval['period'].iloc[i] == period and cond == "nlrp3" and df_periodic_pval['p'].iloc[i] < pval :
                cosinor_rhythms[1].append(name) 
        new_row.update({"nb_prots_ctrl" : [len(cosinor_rhythms[0])]})
        new_row.update({"nb_prots_nlrp3" : [len(cosinor_rhythms[1])]})
        
        nb_prots_both = len([n for n in cosinor_rhythms[0] if n in cosinor_rhythms[1]])
        new_row.update({"nb_prots_both" : [nb_prots_both]})
        new_row = pd.DataFrame(new_row)
        rhythm_prots_period = pd.concat([rhythm_prots_period,new_row], ignore_index=True)

    filepath = os.path.join("data/Cosinor", 'rhythm_prots_period_'+str(period)+'_pval.csv')
    rhythm_prots_period.to_csv(filepath,index=False) 
    
if __name__ == "__main__":
    pool = Pool(7)
    results = pool.map(rhythm_prots_pval,Range_period)
    print(results)
    

def rhythm_prots_qval(period) :
    rhythm_prots_period = pd.DataFrame(columns=['qval','nb_prots_ctrl','nb_prots_nlrp3','nb_prots_both'])
    for qval in np.arange(0,0.06,0.001) :
        cosinor_rhythms = [[], []]
        new_row = {"qval" : [qval]}
        for i in range(len(df_periodic_qval)):
            [name, cond] = df_periodic_qval['test'].iloc[i].split('/')
            if df_periodic_qval['period'].iloc[i] == period and cond == "ctrl" and df_periodic_qval['q'].iloc[i] < qval :
                cosinor_rhythms[0].append(name)
            if df_periodic_qval['period'].iloc[i] == period and cond == "nlrp3" and df_periodic_qval['q'].iloc[i] < qval :
                cosinor_rhythms[1].append(name) 
        new_row.update({"nb_prots_ctrl" : [len(cosinor_rhythms[0])]})
        new_row.update({"nb_prots_nlrp3" : [len(cosinor_rhythms[1])]})
        
        nb_prots_both = len([n for n in cosinor_rhythms[0] if n in cosinor_rhythms[1]])
        new_row.update({"nb_prots_both" : [nb_prots_both]})
        new_row = pd.DataFrame(new_row)
        rhythm_prots_period = pd.concat([rhythm_prots_period,new_row], ignore_index=True)

    filepath = os.path.join("data/Cosinor", 'rhythm_prots_period_'+str(period)+'_qval.csv')
    rhythm_prots_period.to_csv(filepath,index=False) 
    
if __name__ == "__main__":
    pool = Pool(7)
    results = pool.map(rhythm_prots_qval,Range_period)
    print(results)