import pandas as pd
from pandas.core.indexes.numeric import IntegerIndex
import numpy as np
from scipy import optimize
from numpy import sqrt
#*****************************************************************************
# LINEAR_REG
#*****************************************************************************

def chisq( X, Y, SY, a, b):
    chisq = 0
    
    for x, y, sy in zip( X, Y, SY ):
        chisq += ((y - a - x*b)/sy)**2
    return chisq

def fitfunc(x,a,b):
    return a + b*x

def linear_reg( X, Y): 
    
    sigma_X = [0.4]*len(X) # 400 um std
    if ((X[0]-X[1]) == 0):
        mGuess = 100
    else: 
        mGuess = (Y[0]-Y[1])/(X[0]-X[1])
    qGuess = Y[0]-mGuess*X[0]
    p_init = [qGuess, mGuess] # valori iniziali 
    p_best, pcov = optimize.curve_fit( fitfunc, Y, X, sigma=sigma_X, p0=p_init)

    chi2 = chisq(Y, X, sigma_X, p_best[0], p_best[1])
    dof   = len(X) -2#- 1
    chisq_comp = abs( chi2 - dof )/sqrt(2*dof)

    m = p_best[1]
    q = p_best[0]
    return { 'm':m, 'q':q, 'chisq_comp':chisq_comp }
    

#*****************************************************************************
# COMBINATE LOCAL
#*****************************************************************************

from itertools import combinations

def compute(df, useTrigger):
    
    comb = []
    if ( len(df.LAYER.unique())==3 ):
        comb.append(df)
        tot_Hits = 3
    else:
        for index in list(combinations(df.index,4)):
            tmp_df = df.loc[index,:]
            if ( len(tmp_df.LAYER.unique())==4 ): comb.append( tmp_df ) # comb[] contains combinations of data
        tot_Hits = 4

    #saving ORBIT_CNT
    orbit = np.array(df['ORBIT_CNT'])[0]
    
    #saving SL
    sl = np.array(df['SL'])[0]
    
    if ( useTrigger):
        #saving T_ANGLE
        t_angle = np.array(df['T_ANGLE'])[0]
        
        #saving T_POSITION
        t_position = np.array(df['T_POSITION'])[0]

    flag = True

    for data in comb: 
        X = np.array(pd.concat([data['X_RIGHT_GLOB'], data['X_LEFT_GLOB']]))
        Y = np.array(pd.concat([data['WIRE_Z_GLOB'], data['WIRE_Z_GLOB']]))
        for indexes_comb in list(combinations(range(len(X)), tot_Hits)):
                new_X = []
                new_Y = []
                for i in indexes_comb:
                    new_X.append(X[i])
                    new_Y.append(Y[i])
                    
                if(len(np.unique(new_Y))==tot_Hits) :
                    regr_tuple = linear_reg(new_X, new_Y)
                    if ( flag ):  
                        min_lambda = abs(regr_tuple['chisq_comp'])
                        xdata = new_X
                        ydata = new_Y
                        res_dict = regr_tuple
                        flag = False
                        best_comb = indexes_comb
                        best_data = data
                    elif ( abs(regr_tuple['chisq_comp']) < min_lambda ):
                        best_comb = indexes_comb
                        min_lambda = abs(regr_tuple['chisq_comp'])
                        xdata = new_X
                        ydata = new_Y
                        res_dict = regr_tuple
                        best_data = data

    big_df = pd.concat([best_data, best_data], axis=0, ignore_index=True)
    reco_df = big_df.loc[best_comb, :]
    reco_df['m'] = np.full(len(reco_df), res_dict['m']) 
    res_dict['ORBIT_CNT'] = orbit
    res_dict['SL'] = sl
    if (useTrigger):
        res_dict['T_ANGLE'] = t_angle
        res_dict['T_POSITION'] = t_position
    return res_dict, xdata, ydata, reco_df

#*****************************************************************************
# COMBINATE GLOBAL and PARTIAL
#*****************************************************************************

def compute_tot(df_E, data, useTrigger):
    
    res_df = pd.DataFrame()
    
    #saving ORBIT_CNT
    orbit = np.array(df_E['ORBIT_CNT'])[0]
    
    #saving SL
    sl = 4 
    if ( useTrigger):
        #saving T_ANGLE
        t_angle = np.array(df['T_ANGLE'])[0]
        
        #saving T_POSITION
        t_position = np.array(df['T_POSITION'])[0]
        
    X1 = np.array(data.X)
    Y1 = np.array(data.Z)

    if (len(X1)==0 or len(Y1)==0): return
    res_dict = linear_reg(X1, Y1)

    res_dict['ORBIT_CNT'] = orbit
    res_dict['SL'] = sl
    if ( useTrigger):
        res_dict['T_ANGLE'] = t_angle
        res_dict['T_POSITION'] = t_position
    res_df=res_df.append(res_dict, ignore_index=True)

    return res_df


#*****************************************************************************
# COMPUTE EVENT
#*****************************************************************************

class Hitsdf():
    def __init__(self, dataframe = None):
        self.df    = dataframe
    def load(self, dataframe):
        self.df    = dataframe
        
        
def computeEvent(df_E, df_Hits, useTrigger = True):
    
    res_df = pd.DataFrame()
    hits_data = []
    
    chamber=[]
    chamber.append( df_E[ df_E['SL'] == 0] )
    #chamber.append( df_E[ df_E['SL'] == 1] )
    chamber.append( df_E[ df_E['SL'] == 2] )
    chamber.append( df_E[ df_E['SL'] == 3] )

    event_reco_df = pd.DataFrame()

    for df in chamber:
        if( len(pd.unique(df.LAYER)) < 3 ): continue

        regrtuple, xvec, yvec, chamber_reco_df = compute(df, useTrigger)
        event_reco_df = pd.concat([event_reco_df, chamber_reco_df], axis = 0, ignore_index = True)
        hit_data = pd.DataFrame({'X':xvec, 'Z':yvec})
        hits_data.append(hit_data)
        res_df=res_df.append(regrtuple, ignore_index=True)

    if (len(hits_data) == 0): return

    hits_tot  = pd.concat(hits_data, ignore_index = True)
    regrtuple = compute_tot(df_E, hits_tot, useTrigger)
    res_df    = res_df.append(regrtuple, ignore_index=True)
    df_Hits.load(hits_tot)
    res_df = res_df.astype({'SL' : 'int8'})

    if (res_df is None) or (hits_tot is None) or (event_reco_df is None) : return 

    return [res_df, event_reco_df]


#*****************************************************************************
# GET RESULTS
#*****************************************************************************

def getRecoResults(events, useTrigger = True):
    resultsList = []
    resultsData = []
    resultsHits = []
    resultsDf   = []

    for df_E in events:
        if (len(df_E)>32): 
            continue
        resultsHits.append(Hitsdf())
        reco = computeEvent(df_E, resultsHits[-1], useTrigger)
        if (reco is None): continue
        res_df, event_reco_df = reco[0], reco[1]
        resultsList.append(res_df)
        resultsDf.append(event_reco_df)
        resultsData.append(df_E)
    
    return resultsList, resultsData, resultsHits, resultsDf
