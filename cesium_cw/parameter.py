#!/usr/bin/python2
from __future__ import division
import sys
import numpy as np
import string
from atom import Atom
from math import pow
from parameter_common import *

# D1 line Energy level diagram
#    _   F=4 9 levels -
# _< _B               } group1
# |  _C  F=3 7 levels -
# |
# |A
# |
# |  _   F=4 9 levels group 2
# _< _D
#    _E  F=3 7 levels group 3
# total 32 levels
# in Hz
A=335.116048807e12
B=510.860e6
C=656.820e6
D=4.021776399375e9
E=5.170855370625e9


if __name__ == '__main__':
    filename =  sys.argv[1]
    para = Parameter((4,3),(4,3),0.01,1,
                     (((1,3),(0,4)),((1,3),(0,3)),((1,4),(0,4)),((1,4),(0,3))),
                     (E+A+B,E+A-C,E+D,0),
#                     (20000+B+C,20000,10000,0),
                     {#'nu': [3.35120562842e14-9192411945.14,3.35120562842e14-220526.367456],
                      'nup': [((1,4,0),(0,4,-1)),((1,4,0),(0,3,1))],
                      'e_amp': [(1,(1,)), (1,(-1,))],
                      'sweep_profile':[1,-1.5E5,1.5E5,5000],
                      },
                     filename)
    para.write()
    


