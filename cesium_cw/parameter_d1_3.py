#!/usr/bin/python2
from __future__ import division
import sys
import numpy as np
import string
from atom import Atom
from math import pow

# D1 line Energy level diagram
# _C     F=3 7 levels group 1
# |
# |F
# |
# |  _   F=4 9 levels group 2
# _< _D
#    _E  F=3 7 levels group 3
# in Hz
# total 23 levels

A=335.116048807e12
B=510.860e6
C=656.820e6
D=4.021776399375e9
E=5.170855370625e9
F=A-C

def omega(parameter):
    parameter['omega']=[]
    for i in range(7):
        parameter['omega'].append(F+E)
    for i in range(9):
        parameter['omega'].append(E+D)
    for i in range(7):
        parameter['omega'].append(0)
    return parameter

def index2lfm(n):
    #7 9 7
    #0~6 7~15 16~22
    # note the order of m
    if n < 7:
        l = 1
    else:
        l = 0
    if n<7:
        f = 3
        m = n-3
    elif n>6 and n<16:
        f = 4
        m = n - 7 - 4
    else:
        f=3
        m = n - 16 - 3
    return l,f,m

def lfm2index(l,f,m):
    if l==0:
        index = 7
        if f == 4:
            index += m + 4
        elif f ==3:
            index += m + 3 + 9
    else:
        index = 0
        index += m + 3
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
                         'J2':1.0/2.0,
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
                         'J2':1.0/2.0,
                         'I':7.0/2.0}
                parameter['dipole'][i][j] = cs.dipole_element(**coef)
            else:
                parameter['dipole'][i][j] = 0.0
    return parameter

def decoherence(parameter):
    Gamma = 2*np.pi*4.575e6 #this is parameter for D1 line
    n=parameter['n']
    parameter['decoherence_matrix'] = [[[] for i in range(n)] for j in range(n)]
    cs = Atom()
    """
    e l=1 f=3 g l=0 f=4
    e l=1 f=3 g l=0 f=3
    """
    # for f in (3.0,4.0):
    #     for m in range(-1*f,f+1):
    #         for m2 in range(-3,4)
    #             for q in (-1.0,0.0,1.0):
    #                   coef1 = {'q':q,
    #                            'L1':0,
    #                            'L2':1,
    #                            'F1':f,
    #                            'F2':3,
    #                            'mf1':m,
    #                            'mf2':m2,
    #                            'J1':1.0/2.0,
    #                            'J2':1.0/2.0,
    #                            'I':7.0/2.0}
    #                   tmp = Gamma * pow(cs.cg_coef(**coef1),2)
    #                   parameter['decoherence_matrix'][i][j].append([lfm2index(pair[0][0],pair[0][1],d1[2]+q),lfm2index(pair[0][0],pair[0][1],d2[2]+q),tmp])                  
    #egpair=(((1,3),(0,4)),((1,3),(0,3)))
    egpair=(((1,3),(0,4)),)
    for pair in egpair:
        for i in range(n):
            for j in range(i,n):
                d1 = index2lfm(i)
                d2 = index2lfm(j)
                if d1[0:2] == pair[0] and d2[0:2] == pair[0]:
                    parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma])
                elif d1[0:2] == pair[0] and d2[0:2] == pair[1]:
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
                                     'J2':1.0/2.0,
                                     'I':7.0/2.0}
                            coef2 = {'q':q,
                                     'L1':0,
                                     'L2':1,
                                     'F1':pair[1][1],
                                     'F2':pair[0][1],
                                     'mf1':d2[2],
                                     'mf2':d2[2]+q,
                                     'J1':1.0/2.0,
                                     'J2':1.0/2.0,
                                     'I':7.0/2.0}
                            tmp = Gamma*cs.cg_coef(**coef1)*cs.cg_coef(**coef2)
                            if tmp != 0.0:
                                parameter['decoherence_matrix'][i][j].append([lfm2index(pair[0][0],pair[0][1],d1[2]+q),lfm2index(pair[0][0],pair[0][1],d2[2]+q),tmp])
    return parameter

if __name__ == '__main__':
    filename =  sys.argv[1]
    parameter={'nu': [F-D,F+E],
               'level_group': [[0,1,2,3,4,5,6],[7,8,9,10,11,12,13,14,15],[16,17,18,19,20,21,22]],
               'e_amp': [100, 100],
               'n': 23,
               'sweep_profile':[0,-1E11,1E11,400]
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
