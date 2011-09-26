#!/usr/bin/python2
from __future__ import division
import sys
from parameter_common import Parameter

# D2 line Energy level diagram
# refer to page 25 of D line data IMPORTANT the graph underneath is not accurate
#     _F=5
# _<B _F=4 C
# |   _F=3 D  32 levels
# |   _F=2 E
# |A
# |
# |  _F=4 F  9 levels
# _< _
#    _F=3 G  7 levels
# total 48 levels
# in Hz
A=351.72571850e12
B=12.79851e6
C=263.8906e6
D=188.4885e6
E=399.7128e6
F=4.021776399375e9
G=5.170855370625e9

if __name__ == '__main__':
    filename =  sys.argv[1]
    arg_dict = {
        'l1f':(5,4,3,2),
        'l0f':(4,3),
        'B':0.0,
        'd1':0,
        'gamma':5.0,
        'egpair':(((1,5),(0,4)),((1,5),(0,3)),((1,4),(0,4)),((1,4),(0,3)),((1,3),(0,4)),((1,3),(0,3)),((1,2),(0,4)),((1,2),(0,3))),
        'omega_list':(G+A+C,G+A+B,G+A-D,G+A-E,G+F,0),
        'parameter':{'nup': [((1,3,0),(0,4,-1)),((1,3,0),(0,3,1))],
                      'e_amp': [(5,(-1,1)), (5,(-1,1))],
                      'sweep_profile':[0,-4E5,4E5,500],
                      },
        'filename':filename
        }
    para = Parameter(**arg_dict)
    para.write()
