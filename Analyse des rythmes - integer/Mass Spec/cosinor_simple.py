import pandas as pd
import numpy as np
import scipy.stats as stats
import math
import statistics
import statsmodels.regression.linear_model as sm

# fit_group consist in providing the data for each test and then storing the final results in the dataframe
def fit_group(df, period = 24):
    df_results_single_per = pd.DataFrame(columns = ['test', 'period','n_components','p','q','RSS','amplitude','acrophase', 'mesor','n_points'], dtype=float)
    # names store all the test in df, ex PER2/ctrl
    names = np.unique(df.test) 

    for test in names:
        # X store the circadian times 
        # Y store the values 
        X, Y = np.array(df[df.test == test].x), np.array(df[df.test == test].y)  
        
        # Do a function to calculate the results, the parameters and the statistics
        parameters, statistics = fit_me(X, Y, period = period, name = test)
        
        # Store the results in the df_results_single_per dataframe 
        df_results_single_per = df_results_single_per.append({'test': test, 
                                            'period': period,
                                            'n_components': 1,
                                            'p': statistics['p'], 
                                            'q': np.nan,
                                            'RSS': statistics['RSS'], 
                                            'amplitude': parameters['amplitude'],
                                            'acrophase': parameters['acrophase'],
                                            'mesor': parameters['mesor'],
                                            'n_points' : len(X)
                                            }, ignore_index=True)

    return df_results_single_per


# fit_me consist in providing all the cosinor results, statistics and parameters
def fit_me(X, Y, period = 24, name = ''):
    
   
    # FIRST STEP : generate X_fit = cos(2pit/T) and Z_fit = sin(2pit/T) 
    # Both data are combined into the data_fit dataframe
    X_fit = np.cos(2*np.pi*X/period)
    Z_fit = np.sin(2*np.pi*X/period) 
    
    data_fit = pd.DataFrame({"constant" : [1]*len(X_fit),"X_fit" : X_fit, "Z_fit" : Z_fit})  
    

    # SECOND STEP : calculate the parameters estimation
    # û = S^(-1)*d
    # d =np.array([[np.sum(Y)],
    #               [np.sum(Y*X_fit)],
    #               [np.sum(Y*Z_fit)]])
    
    # S = np.array([[len(X), np.sum(X_fit), np.sum(Z_fit)],
    #               [np.sum(X_fit), np.sum(X_fit*X_fit), np.sum(X_fit*Z_fit)],
    #               [np.sum(Z_fit), np.sum(X_fit*Z_fit), np.sum(Z_fit*Z_fit)]])
    
    # u_hat = np.dot(inv(S),d)
    
    # This method can lead to some issues with singular matrix,
    # we use a multiple linear regression to estimate M_hat, beta_hat gamma_hat.
    
    model = sm.OLS(Y, data_fit)
    results = model.fit()
    Y_fit = results.fittedvalues
    
    MESOR = results.params["constant"]
    AMPLITUDE = (np.square(results.params["X_fit"]) + np.square(results.params["Z_fit"]))**(1/2)
    
    if results.params["X_fit"] == 0 :
        ACROPHASE = np.nan
    
    elif results.params["Z_fit"] == 0 :
        ACROPHASE = 0
        
    elif results.params["Z_fit"] > 0 and results.params["X_fit"] > 0 : 
        ACROPHASE = - np.arctan(abs(results.params["Z_fit"]/results.params["X_fit"]))
    
    elif results.params["Z_fit"] > 0 and results.params["X_fit"] < 0 : 
        ACROPHASE = - np.pi + np.arctan(abs(results.params["Z_fit"]/results.params["X_fit"]))
    
    elif results.params["Z_fit"] < 0 and results.params["X_fit"] > 0 : 
        ACROPHASE = - 2*np.pi + np.arctan(abs(results.params["Z_fit"]/results.params["X_fit"]))
    
    elif results.params["Z_fit"] < 0 and results.params["X_fit"] < 0 : 
        ACROPHASE = - np.pi - np.arctan(abs(results.params["Z_fit"]/results.params["X_fit"]))

    

    parameters = {'period': period, 'amplitude': AMPLITUDE, 'acrophase': ACROPHASE, 'mesor':MESOR}
    
    
    
    # THIRD STEP :  calculate the statistics of the model
    MSS = sum((Y_fit - Y.mean())**2) 
    RSS = sum((Y - Y_fit)**2) 
    
    if len(X) != 3 and RSS != 0 : 
        F = (MSS/2) / (RSS/(len(X) - 3)) 
        p = 1 - stats.f.cdf(F, 2, len(X) - 3)
        statistics = {'p':p, 'RSS': RSS}
    
    else :
        statistics = {'p': np.nan, 'RSS': np.nan}


    return parameters, statistics


# get_best_fits consists in providing the model with the best period for each test 
def get_best_fits(df_results):
    df_best = pd.DataFrame(columns = df_results.columns, dtype=float)
    names = np.unique(df_results.test)
    
    for name in names:
        df_test_models = df_results[df_results.test == name]
        
        # If there is less than 3 n_points we consider that the periodic model is not enough significative
        if df_test_models.n_points.iloc[0] > 3 : 
            # Best model is the one with the minimum RSS
            M = df_test_models.RSS.min()
            best_row_period =  df_test_models[df_test_models.RSS == M]
            
            # If there is more than one case with the same RSS, we peak the middle model
            if len(best_row_period) > 1 :
                if len(best_row_period)%2 !=0 :
                    median = math.ceil(statistics.median(best_row_period.period))
                    best_row_period = best_row_period[best_row_period.period == median]
                else : 
                    median = best_row_period.period.iloc[int(len(best_row_period)/2)-1]
                    best_row_period = best_row_period[best_row_period.period == median]
            
            # Then, we check if the periodic model is better than the MESOR one
            if float(best_row_period["p"]) > 0.05 :
                best_row = {"test" : name, "n_components" : 0, "period" : 0, "RSS" : best_row_period.RSS,
                "mesor" : best_row_period.mesor, "amplitude" : np.nan, "acrophase" : np.nan, "p" : np.nan,
                "q" : np.nan, "n_points" : best_row_period.n_points}
                best_row = pd.DataFrame(best_row)
            else : 
                best_row = best_row_period
    
            df_best = df_best.append(best_row, ignore_index = True)
            
        else : 
            best_row = {"test" : name, "n_components" : [0], "period" : [0], "RSS" : np.nan,
            "mesor" : np.nan, "amplitude" : np.nan, "acrophase" : np.nan, "p" : np.nan,
            "q" : np.nan, "n_points" : [df_test_models.n_points.iloc[0]]}
            best_row = pd.DataFrame(best_row)
            
            df_best = df_best.append(best_row, ignore_index = True)
            
    return df_best


# get_periodic_pval consists in providing the significative models for each protein/cond with a pval threshold
def get_periodic_pval(df_results):
    df_periodic = pd.DataFrame(columns = df_results.columns, dtype=float)
    names = np.unique(df_results.test)
    
    for name in names :
        df_test_models = df_results[df_results.test == name]
        for i in range(len(df_test_models)) :
            new_row = df_test_models.iloc[i,:]
            if float(new_row["p"]) < 0.05 :
                df_periodic = df_periodic.append(new_row, ignore_index = True)
                
    return df_periodic

# get_periodic_qval consists in providing the significative models for each protein/cond with a qval threshold
def get_periodic_qval(df_results):
    df_periodic = pd.DataFrame(columns = df_results.columns, dtype=float)
    names = np.unique(df_results.test)
    
    for name in names :
        df_test_models = df_results[df_results.test == name]
        for i in range(len(df_test_models)) :
            new_row = df_test_models.iloc[i,:]
            if float(new_row["q"]) < 0.05 :
                df_periodic = df_periodic.append(new_row, ignore_index = True)
                
    return df_periodic
        
