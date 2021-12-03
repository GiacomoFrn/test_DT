"""
File containing usefull constants
"""

# times
DURATION_BX = 25        # BX duration in ns
DURATION_TDC = 25./30.  # TDC duration in ns
ORBIT_BX = 3564         # number of BX in an orbit

# geometry
XCELL = 42. #  cell width in mm
ZCELL = 13. # cell height in mm

Z_POS_SHIFT = [ -9999, -ZCELL*1.5,   -ZCELL*0.5,     ZCELL*0.5,     ZCELL*1.5  ]
X_POS_SHIFT = [ -9999, -XCELL*7,       -XCELL*7.5,    -XCELL*7,     -XCELL*7.5 ] 

# max drift time in BX
TM = 15.6

# max drift time in ns
TDRIFT_NS = TM*DURATION_BX

# drift velocity
VDRIFT  = XCELL*0.5 / TDRIFT_NS 