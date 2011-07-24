#!/usr/bin/python2
import sys
import numpy as np
from atom import Atom
filename =  sys.argv[1]
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
# in Hz
A=335.116048807e12
B=510.860e6
C=656.820e6
D=4.021776399375e9
E=5.170855370625e9

def omega(parameter):
    parameter['omega']=[]
    for i in range(9):
        parameter['omega'].append(E+A+B)
    for i in range(7):
        parameter['omega'].append(E+A-C)
    for i in range(9):
        parameter['omega'].append(E+D)
    for i in range(7):
        parameter['omega'].append(0)
    return parameter

def index2lfm(n):
    #9 7 9 7
    #0~8 9~15 16~24 25~31
    # note the order of m
    if n < 16:
        l = 1
    else:
        l = 0
    if n<9:
        f = 4
        m = n-4
    elif n>8 and n< 16:
        f = 3
        m = n - 9 - 3
    elif n>15 and n<25:
        f = 4
        m = n - 16 - 4
    else:
        f=3
        m = n - 25 - 3
    return l,f,m

def lfm2index(l,f,m):
    if l==0:
        index = 16.0
    else:
        index = 0.0
    if f == 3.0:
        index += 9.0
    index += m+f
    return index

def dipole(parameter):
    n=parameter['n']
    #parameter['dipole'] = np.zeros([n,n])
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
    for i in range(n):
        for j in range(n):
            d1 = index2lfm(i)
            d2 = index2lfm(j)
#            if d1[0] != d2[0]:
                # defind ground state and excited state
            if d1[0] == 0: #excited state is the one with l=1
                g = d1[1] 
                e = d2[1]
            else:
                g = d2[1]
                e = d1[1]
            if d1[1] == e and d2[1] == e and d1[0] == 1 and d1[0] == 1 :
                parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma])
            elif d1[1] == e and d1[0] == 1 and d2[1] == g and d2[0] == 0 :
                parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
            elif d1[1] == g and d1[0] == 0 and d2[1] == e and d2[0] == 1 :
                parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
            elif d1[1] == g and d2[1] == g:
                for q in [-1.0,0.0,1.0]:
                    if (d1[2]+q <= d1[1] and d1[2]+q >= -1*d1[1]) and (d2[2]+q <= d2[1] and d2[2]+q >= -1*d2[1]):
                        coef1 = {'q':q,
                                 'L1':0,
                                 'L2':1,
                                 'F1':g,
                                 'F2':e,
                                 'mf1':d1[2],
                                 'mf2':d1[2]+q,
                                 'J1':1.0/2.0,
                                 'J2':1.0/2.0,
                                 'I':7.0/2.0}
                        coef2 = {'q':q,
                                 'L1':0,
                                 'L2':1,
                                 'F1':g,
                                 'F2':e,
                                 'mf1':d2[2],
                                 'mf2':d2[2]+q,
                                 'J1':1.0/2.0,
                                 'J2':1.0/2.0,
                                 'I':7.0/2.0}
                        tmp = Gamma*cs.cg_coef(**coef1)*cs.cg_coef(**coef2)
                        if tmp != 0.0:
                            parameter['decoherence_matrix'][i][j].append([lfm2index(1.0,e,d1[2]+q),lfm2index(1.0,e,d2[2]+q),tmp])
    return parameter

if __name__ == '__main__':
    parameter={'nu': [A+E,A-D],
               'level_group': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], [16, 17, 18, 19, 20, 21, 22, 23, 24],  [25, 26, 27, 28, 29, 30, 31]],
               'e_amp': [100, 100],
               'n': 32,
               'sweep_profile':[0,-1E11,1E11,400]
                }
    omega(parameter)
    dipole(parameter)
    decoherence(parameter)
    txtf = open(filename+'.txt','w')
    txtf.write(str(parameter))
    txtf.close()

