import pandas as pd
import numpy as np
from modules.constants import DURATION_BX, DURATION_TDC, ORBIT_BX, XCELL, ZCELL, Z_POS_SHIFT, X_POS_SHIFT, TM, TDRIFT_NS, VDRIFT    

def map_hit_position(hits_df_, local=False):
    """
    Compute the position of the hit.
    If local=True the positions left/right will 
    be computed in the local frame of reference
    (e.g. X_LEFT_LOC). Else, the positions will be
    computed in the global frame of reference
    """

    hits_df_['HIT_DRIFT_TIME'] = (hits_df_['BX_COUNTER']+hits_df_['TDC_MEAS']/30)*25-hits_df_['T0_NS']

    if local:
        ref = 'LOC'
    else:
        ref = 'GLOB'

    hits_df_[f'X_LEFT_{ref}'] = hits_df_[f'WIRE_X_{ref}'] - hits_df_['HIT_DRIFT_TIME']*VDRIFT
    hits_df_[f'X_RIGHT_{ref}'] = hits_df_[f'WIRE_X_{ref}'] + hits_df_['HIT_DRIFT_TIME']*VDRIFT

    return hits_df_

def getEvents(stream_df, cfg, runTimeShift, useTrigger):
    
    from modules.mappings import Mapping
    # create a dataframe with only valid hits ->
    # trigger words and scintillator hits are removed
    hits_df = stream_df[
        (stream_df.HEAD == cfg['headers']['valid_hit']) &
        (stream_df.TDC_CHANNEL <= 127)
        ].copy()

    # fix TDC_MEAS data type
    hits_df = hits_df.dropna()
    hits_df = hits_df.astype({'TDC_MEAS': 'int32'})

    # create mapping with the loaded configurations
    mapper = Mapping(cfg)

    # map hits
    hits_df = mapper.global_map(hits_df)
    
    # select all orbits with a trigger signal from
    # the scintillators coincidence
    trigger_df = stream_df[
        (stream_df['HEAD'] == cfg['scintillator']['head']) & 
        (stream_df['FPGA'] == cfg['scintillator']['fpga']) & 
        (stream_df['TDC_CHANNEL'] == cfg['scintillator']['tdc_ch'])
    ].copy()
    
    # create a T0 column (in ns)
    trigger_df['T0'] = (trigger_df['BX_COUNTER'] + trigger_df['TDC_MEAS'] / 30)
    
    # select only hits in the same orbit of a scint trigger signal
    hits_df_ = pd.merge(
        hits_df, trigger_df[['ORBIT_CNT','T0']],
        left_on='ORBIT_CNT', right_on='ORBIT_CNT',
        suffixes=(None, None)
    )
    print(len(hits_df))
    
    
    # if true consider only the hits in the same orbit of a scint trigger AND FPGA trigger signals
    if (useTrigger == True):
        
        # select all orbits with a trigger signal from
        # the FPGA , saving also trigger angle and position
        trigger_df_ = stream_df[
            (stream_df['HEAD'] == cfg['headers']['t0_trg']) & 
            (stream_df['FPGA'] == 0)
        ].copy()

        triggerAngle_df = stream_df[
            (stream_df['HEAD'] == 4) & 
            (stream_df['FPGA'] == 0)
        ].copy()

        triggerPosition_df = stream_df[
            (stream_df['HEAD'] == 5) & 
            (stream_df['FPGA'] == 0)
        ].copy()
        
        trigger_df_['T_T0'] = (trigger_df_['BX_COUNTER']*25 + trigger_df_['TDC_MEAS']*25 / 30)
        
        triggerAngle_df.rename(columns={'TDC_MEAS': 'T_ANGLE'}, inplace=True)
        trigger_df_ = pd.merge(
            trigger_df_, triggerAngle_df[['ORBIT_CNT','T_ANGLE']],
            left_on='ORBIT_CNT', right_on='ORBIT_CNT',
            suffixes=(None, None)
        )

        triggerPosition_df.rename(columns={'TDC_MEAS': 'T_POSITION'}, inplace=True)
        trigger_df_ = pd.merge(
            trigger_df_, triggerPosition_df[['ORBIT_CNT','T_POSITION']],
            left_on='ORBIT_CNT', right_on='ORBIT_CNT',
            suffixes=(None, None)
        )

        hits_df_ = pd.merge(
            hits_df_, trigger_df_[['ORBIT_CNT','T_T0', 'T_ANGLE', 'T_POSITION']],
            left_on='ORBIT_CNT', right_on='ORBIT_CNT',
            suffixes=(None, None)
        )
        
    # TIME SHIFTING
    
    
    # apply scintillator calibration -> they shoul be computed
    # for each run! For this example run, values are provided in
    # the configurations file.

    # create a time column in NS
    hits_df_['T0_NS'] = hits_df_['T0'] * DURATION_BX

    for sl, offset_sl in cfg['time_offset_sl'].items():
        # correction is in the form:
        # coarse offset valid for all SL + fine tuning 
        hits_df_.loc[
            hits_df_['SL'] == sl, 'T0_NS'
        ] -= cfg['time_offset_scint'] -runTimeShift + offset_sl # 6 for 123, 8 for 0
    
    # having the time allows us to compute hit position 
    # with left / right ambiguity

    hits_df_ = map_hit_position(hits_df_, local=False)
    
    # creat a list of event ( an event is a group of hits all in the same orbit
    # the probability of multiples hits in the same orbit is very low)
    events=[]
    for x in pd.unique(hits_df_.ORBIT_CNT) :
        events.append(hits_df_[hits_df_['ORBIT_CNT']==x])
        events[-1] = events[-1].reset_index(drop = True)
    return events