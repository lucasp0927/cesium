#!/usr/bin/python2
from __future__ import division
import sys
import numpy as np
import string
from atom import Atom
from math import pow

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


def omega(parameter):
    parameter['omega']=[]
    for i in range(11):
        parameter['omega'].append(G+A+C)
    for i in range(9):
        parameter['omega'].append(G+A+B)
    for i in range(7):
        parameter['omega'].append(G+A-D)
    for i in range(5):
        parameter['omega'].append(G+A-E)
    for i in range(9):
        parameter['omega'].append(G+F)
    for i in range(7):
        parameter['omega'].append(0)
    return parameter

def index2lfm(n):
    """
    11   9     7     5     9     7
    0~10 11~19 20~26 27~31 32~40 41~47
    """
    if n < 32:
        l = 1
        if n < 11:
            f = 5
            m = n-5
        elif n > 10 and n < 20:
            f = 4
            m = n - 11 - 4
        elif n > 19 and n < 27:
            f = 3
            m = n - 20 - 3
        elif n > 26 and n < 32:
            f = 2
            m = n - 27 - 2
    else:
        l = 0
        if n > 31 and n < 41:
            f = 4
            m = n - 32 - 4
        if n > 40:
            f = 3
            m = n - 41 - 3
    return l,f,m

def lfm2index(l,f,m):
    if l == 0:
        index = 32.0
        if f == 3.0:
            index += 9.0 + m + f
        else:
            index += m + f
    else:
        index = 0.0
        if f == 5.0:
            index += m + f
        elif f == 4.0:
            index += 11 + m + f
        elif f == 3.0:
            index += 20 + m + f
        elif f == 2.0:
            index += 27 + m + f
    return index

def dipole(parameter):
    n=parameter['n']
    parameter['dipole'] = [[0 for i in range(n)] for j in range(n)]
    cs = Atom()
    for i in range(n):
        for j in range(n):
            d1 = index2lfm(i)
            d2 = index2lfm(j)
            if d1[0] == 0 and d2[0] == 1:
                coef = {'q':0,
                         'L1':0,
                         'L2':1,
                         'F1':d1[1],
                         'F2':d2[1],
                         'mf1':d1[2],
                         'mf2':d2[2],
                         'J1':1.0/2.0,
                         'J2':3.0/2.0,
                         'I':7.0/2.0}
                parameter['dipole'][i][j] = cs.dipole_element(**coef)
            elif d2[0] == 0 and d1[0] == 1:
                coef = {'q':0,
                         'L1':0,
                         'L2':1,
                         'F1':d2[1],
                         'F2':d1[1],
                         'mf1':d2[2],
                         'mf2':d1[2],
                         'J1':1.0/2.0,
                         'J2':3.0/2.0,
                         'I':7.0/2.0}
                parameter['dipole'][i][j] = cs.dipole_element(**coef)
            else:
                parameter['dipole'][i][j] = 0.0
    return parameter

def decoherence(parameter):
    Gamma = 2*np.pi*5.234e6 #this is parameter for D2 line
    n=parameter['n']
    parameter['decoherence_matrix'] = [[[] for i in range(n)] for j in range(n)]
    cs = Atom()
    """
    e l=1 f=3 g l=0 f=4
    e l=1 f=3 g l=0 f=3
    """
    egpair=(((1,5),(0,4)),((1,5),(0,3)),((1,4),(0,4)),((1,4),(0,3)),((1,3),(0,4)),((1,3),(0,3)),((1,2),(0,4)),((1,2),(0,3)))
#    egpair=(((1,5),(0,4)),)
    
    for pair in egpair:
        for i in range(n):
            for j in range(i,n):
                d1 = index2lfm(i)
                d2 = index2lfm(j)
                if d1[0:2] == pair[0] and d2[0:2] == pair[0] and i != j:
                     parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma])
                if d1[0:2] == pair[0] and d2[0:2] == pair[1]:
                    parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
                elif d1[0:2] == pair[1] and d2[0:2] == pair[0]:
                    parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
                elif d1[0:2] == pair[1] and d2[0:2] == pair[1]:
                    for q in (-1.0,0.0,1.0):
                        f1 = pair[0][1]
                        if (d1[2]+q <= f1 and d1[2]+q >= -1*f1) and (d2[2]+q <= f1 and d2[2]+q >= -1*f1):
                            coef1 = {'q':q,
                                     'L1':0,
                                     'L2':1,
                                     'F1':pair[1][1],
                                     'F2':pair[0][1],
                                     'mf1':d1[2],
                                     'mf2':d1[2]+q,
                                     'J1':1.0/2.0,
                                     'J2':3.0/2.0,
                                     'I':7.0/2.0}
                            coef2 = {'q':q,
                                     'L1':0,
                                     'L2':1,
                                     'F1':pair[1][1],
                                     'F2':pair[0][1],
                                     'mf1':d2[2],
                                     'mf2':d2[2]+q,
                                     'J1':1.0/2.0,
                                     'J2':3.0/2.0,
                                     'I':7.0/2.0}
                            tmp = Gamma*cs.cg_coef(**coef1)*cs.cg_coef(**coef2)
                            if tmp != 0.0:
                                ii = lfm2index(pair[0][0],pair[0][1],d1[2]+q)
                                jj = lfm2index(pair[0][0],pair[0][1],d2[2]+q)
                                parameter['decoherence_matrix'][i][j].append([ii,jj,tmp])
                                if ii == jj:
                                    parameter['decoherence_matrix'][int(ii)][int(jj)].append([ii,jj,-1*tmp])
    return parameter

if __name__ == '__main__':
    filename =  sys.argv[1]
    parameter={'nu': [A-F,A+G],
               'level_group': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31], [32, 33, 34, 35, 36, 37, 38, 39, 40],  [41, 42, 43, 44, 45, 46, 47]],
               'e_amp': [100, 100],
               'n': 48,
               'sweep_profile':[0,-1E10,1E10,400]
                }
    omega(parameter)
    dipole(parameter)
    decoherence(parameter)
    print "      L   F   M|"
    print "----------------------------------------"
    sum = 0.0
    psum = 0.0
    for i in range(parameter['n']):
        iter = parameter['decoherence_matrix'][i][i]
        for j in iter:
            sum += j[2]
            psum += j[2]
        print "{:3d}:{:3d} {:3d} {:3d}| {:<100} |Sum is: {:>20f}".format(i,index2lfm(i)[0],index2lfm(i)[1],index2lfm(i)[2],iter,psum)
        psum = 0.0        

    print "the sum is %f" %sum
    txtf = open(filename+'.txt','w')
    txtf.write(str(parameter))
    txtf.close()

