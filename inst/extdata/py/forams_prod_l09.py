# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:59:25 2016
Last modified on Tue Oct 10 20:28:00 CEST 2017
@author: Didier M. Roche
"""

#  !! changes from version 0.2: modified the print name of pachy_d to incompta
#       following (Darling et al., 2006)
#       DARLING, K. F., KUCERA, M., KROON, D. and WADE, C. M., 2006.
#       A resolution for the coiling direction paradox in Neogloboquadrina pachyderma.
#       Paleoceanography 21: PA2011, doi:10.1029/2005PA001189.
#

__version__ = "1.0"

global l09_cnsts_dic

# Order of the coeffictients from Lombard et al., 2009
#  µ(T1) alias mut, TA, TL, TH, TAL, TAH ... last column below is r^2

# O. universa          0.23 (±0.09) 626 (±3851)   289.8 (±3.4) 304.0 (± 0.5) 30953 (±13,119)    249412 (±389,004) 0.75
# G. sacculifer        0.30 (±0.01) 2155 (±686)   289.3 (±0.6) 304.4 (± 0.3) 94385 (±23,740)    171209 (±48,514)  0.95
# G. siphonifera       0.17 (±0.02) 7105 (±1861)  284.9 (±0.3) 301.8 (± 0.5) 349991 (±283,836)  130852 (±66,186)  0.77
# G. ruber             0.24 (±0.11) 1086 (± 5772) 292.3 (±2.3) 303.4 (± 0.5) 53765 (±21,035)    165490 (±102,487) 0.90
# N. dutertrei         0.10 (±0.02) 6876 (±4805)  280.1 (±1.6) 298.8 (± 3.0) 210746 (± 602,206)  57855 (±39,102)  0.63
# G. bulloides         0.25 (±0.02) 6482 (±1256)  281.1 (±0.3) 298.5 (± 0.3) 295374 (±266,947)  159618 (±44,756)  0.95
# N. pachyderma (dex.) 0.13 (±0.01) 6584 (± 1080) 277.8 (±0.3) 293.7 (± 0.5) 300239 (±249,645)  110583 (±35,218)  0.96
# N. pachyderma (sin.) 0.37 (±0.05) 6584          999999.0     279.8 (± 0.7) 999999.0            59491 (±15,098)  0.90



l09_cnsts_dic = {
    'universa'    :[0.23,  626., 289.8, 304.0,  30953.0, 249412.0],
    'sacculifer'  :[0.30, 2155., 289.3, 304.4,  94385.0, 171209.0],
    'siphonifera' :[0.17, 7105., 284.9, 301.8, 349991.0, 130852.0],
    'ruber'       :[0.24, 1086., 292.3, 303.4,  53765.0, 165490.0],
    'dutertrei'   :[0.10, 6876., 280.1, 298.8, 210746.0,  57855.0],
    'bulloides'   :[0.25, 6482., 280.0, 298.5, 295374.0, 159618.0],
    'pachy_d'     :[0.13, 6584., 277.8, 293.7, 300239.0, 110583.0],
    'pachy_s'     :[0.37, 6584., 999.0, 279.8,      0.0,  59491.0]
                }

global claire_depth_dic

# Depth of happy living foram, in meters from the surface ocean
#
# universa, sacculifer, siphonifera, dutertrei I made up to 100 meters

# Optimized version from error analysis

claire_depth_dic = {
    'universa'    :[-100],
    'sacculifer'  :[-100],
    'siphonifera' :[-100],
    'ruber'       :[ -10],
    'dutertrei'   :[-100],
    'bulloides'   :[-400],
    'pachy_d':[ -65],
    'pachy_s':[-550]
                }

foram_names_dic = {
    'ruber':'Globigerinoides ruber',
    'bulloides':'Globigerina bulloides',
    'pachy_s':'Neogloboquadrina pachyderma',
    'pachy_d':'Neogloboquadrina incompta',
    'sacculifer':'Globigerinoides sacculifer'
              }

# Need the maximum growth rates computed from the function given
# Computed once numerically with scipy.optimize functions

global l09_maxgrowth_dic

l09_maxgrowth_dic = {
    'ruber'       : [0.260441142387609],
    'sacculifer'  : [0.3693161596433823],
    'bulloides'   : [0.3170973153524813],
    'pachy_s': [0.04886423011338358],
    'dutertrei'   : [0.10951633092705156],
    'pachy_d': [0.10944195438579478],
    'universa'    : [0.24084903555860263],
    'siphonifera' : [0.2791874260946192]
                    }

def growth_rate_l09(nm_foram, input_T, norm=0):

# Growth rate in day^-1 ...
# Inputs:
#   foram name as string
#   input_T temperature in °K

# Optional input norm is a normalization factor

    import math

    mut, TA, TL, TH, TAL, TAH = l09_cnsts_dic[nm_foram]

    if (norm == 0):
      if (TAL > 0.0):
        return ( mut*math.exp(TA/293.0-TA/input_T)
                 /(1+math.exp(TAL/input_T-TAL/TL)+math.exp(TAH/TH-TAH/input_T))
               )
      else:
        return ( mut*math.exp(TA/293.0-TA/input_T)
                 /(2+math.exp(TAH/TH-TAH/input_T))
               )
    else:
      max_func = l09_maxgrowth_dic[nm_foram][0]

      return ( mut*math.exp(TA/293.0-TA/input_T)
               /(1+math.exp(TAL/input_T-TAL/TL)+math.exp(TAH/TH-TAH/input_T))
             ) / max_func

#enddef growth_rate_l09

def growth_rate_l09_array(nm_foram, input_T, norm=0):

# Growth rate in day^-1 ...
# Inputs:
#   foram name as string
#   input_T temperature in °K

# Optional input norm is a normalization factor

    from numpy import ma
    from numpy import asarray

    mut, TA, TL, TH, TAL, TAH = l09_cnsts_dic[nm_foram]
    
    input_T = asarray(input_T)

    if (norm == 0):
      if (TAL > 0.0):
        return ( mut*ma.exp(TA/293.0-TA/input_T)
                 /(1+ma.exp(TAL/input_T-TAL/TL)+ma.exp(TAH/TH-TAH/input_T))
               )
      else:
        return ( mut*ma.exp(TA/293.0-TA/input_T)
                 /(2+ma.exp(TAH/TH-TAH/input_T))
               )
    else:
      max_func = l09_maxgrowth_dic[nm_foram][0]
      return ( mut*ma.exp(TA/293.0-TA/input_T)
                 /(1+ma.exp(TAL/input_T-TAL/TL)+ma.exp(TAH/TH-TAH/input_T))
             ) / max_func

#enddef growth_rate_l09

def get_living_depth(nm_foram):

    return claire_depth_dic[nm_foram]

#enddef get_living_depth


def get_foram_fullname(short_name):

    return foram_names_dic[short_name]
#enddef get_foram_fullname

#
# The End of All Things (op. cit.)
#
