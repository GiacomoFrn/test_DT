import pandas as pd
import numpy as np
from numpy import arctan
from modules.constants import XCELL
#*****************************************************************************
# TIME RESOLUTION
#*****************************************************************************

def getTimeRes(events, timeRes):
    if (timeRes is None):
        timeRes = []
    for x in events:
        timeRes.append(x.loc[0, 'T0']*25-x.loc[0, 'T_T0'])
    for d in timeRes:
        if (d<-500 or d>1000):
            timeRes.remove(d)
    return timeRes


#*****************************************************************************
# EFFICIENCY
#*****************************************************************************

def efficiency ( S_runLists, runLists, cfg):
    z = cfg['mc_z']
    tmpdf = pd.DataFrame()
    effD = 0
    effN = 0
    for resultsList in runLists:
        for event in resultsList:
    
            if not ( event is  None): 
                tmpdf = event.set_index('SL')
                if( 4 in tmpdf.index ): 
                    x = tmpdf.loc[4,'m']*z+tmpdf.loc[4,'q']
                    if ( (x > (cfg['mc_x_shift'] - XCELL*1.5 )) and (x < (cfg['mc_x_shift'] + XCELL*1.5 )) ):
                        effN += 1
    for resultsList in S_runLists:
        for event in resultsList:
    
            if not ( event is  None): 
                tmpdf = event.set_index('SL')
                if( 4 in tmpdf.index ): 
                    x = tmpdf.loc[4,'m']*z+tmpdf.loc[4,'q']
                    if ( (x > (cfg['mc_x_shift'] - XCELL*1.5 )) and (x < (cfg['mc_x_shift'] + XCELL*1.5 ))  ):
                        effD += 1
    return effN/effD, effN, effD
    
#*****************************************************************************
# POSITION RESOLUTION
#*****************************************************************************

def getPosRes( runLists, cfg ):
    diffPosGlob = []
    diffPosLoc  = []

    z = cfg['mc_z']
    x_shift = cfg['mc_x_shift']
    tmpdf = pd.DataFrame()

    for resultsList in runLists:
        for event in resultsList:
            if( event is None): 
                continue
            else:
                tmpdf = event.set_index('SL')

            if( (2 in tmpdf.index) and (0 in tmpdf.index) and (1 in tmpdf.index) ): 
                diffPosGlob.append(tmpdf.loc[4,'T_POSITION']+x_shift - (tmpdf.loc[4,'m']*z+tmpdf.loc[4,'q']))
                
            if( 1 in tmpdf.index):     
                diffPosLoc.append(tmpdf.loc[4,'T_POSITION']+x_shift - (tmpdf.loc[1,'m']*z+tmpdf.loc[1,'q']))
                
    return diffPosGlob, diffPosLoc

#*****************************************************************************
# ANGLE RESOLUTION
#*****************************************************************************

def getAngleRes( runLists, cfg ):
    diffAngleGlob = []
    diffAngleLoc  = []

    tmpdf = pd.DataFrame()

    for resultsList in runLists:
        for event in resultsList:
            if( event is None): 
                continue
            else:
                tmpdf = event.set_index('SL')

            if( (2 in tmpdf.index) and (0 in tmpdf.index) and (1 in tmpdf.index) ): 
                diffAngleGlob.append(arctan(tmpdf.loc[4,'m'])*10**3- arctan(tmpdf.loc[2,'T_ANGLE'])*10**3)
                
            if( 1 in tmpdf.index):    
                diffAngleLoc.append(arctan(tmpdf.loc[1,'m'])*10**3- arctan(tmpdf.loc[1,'T_ANGLE'])*10**3)
                
    return diffAngleGlob, diffAngleLoc
