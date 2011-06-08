#!/usr/bin/python2
import sys
import os
from density_matrix import System
import numpy as np

def plot(n):
    f=open('./tmp.gp','w')
    f.write('set terminal png\nset output \'graph.png\'\n')
    tmp_str = 'plot \'%s\' using 1:2 with lines'%filename
    for i in range(n-1):
        tmp_str += ', \'%s\' using 1:%d with lines'%(filename,i+3)
    tmp_str += '\n'
    f.write(tmp_str)
    f.close()
    os.system('gnuplot tmp.gp')
    os.remove('tmp.gp')

if __name__ ==  '__main__':
    # n=3
    # #decoherence
    # Gamma1 = 5000000
    # Gamma12 = 2500000
    # Gamma13 = 2500000
    # gamma1 = 10000
    # gamma2 = 10000
    # #in conjugate form
    # #decoherence
    # decoherence_matrix = [[[[0,0,-1*Gamma1]],[[0,1,-0.5*Gamma1]],[[0,2,-0.5*Gamma1]]],
    #                       [[],[[0,0,Gamma12],[1,1,-1*gamma1],[2,2,gamma1]],[[1,2,-1*gamma2]]],
    #                       [[],[],[[0,0,Gamma13],[2,2,-1*gamma2],[1,1,gamma2]]]]


    # parameter = {'n':n,
    #              'omega':[105E10,9E9,0],
    #              'dipole':[[0,1000000,1000000],
    #                        [1000000,0,0],
    #                        [1000000,0,0]],
    #              'nu':[105E10-9E9,105E10], # on resonence
    #              'e_amp':[1,1],#amplitude of electric field
    #              'level_group':[[0],[1],[2]],
    #              'decoherence_matrix': decoherence_matrix
    #              }


    n=4
    #decoherence
    Gamma1 = 9000000
    Gamma12 = 3000000
    Gamma13 = 3000000
    gamma1 = 10000
    gamma2 = 10000
    #in conjugate form
    #decoherence
    decoherence_matrix = [[[[0,0,Gamma12]],[[0,1,-1*Gamma1]],[[0,0,0]],[[0,0,0]]],
                          [[],[[1,1,-1*Gamma1]],[[1,2,-0.5*Gamma1]],[[1,3,-0.5*Gamma1]]],
                          [[],[],[[1,1,Gamma12],[2,2,-1*gamma1],[3,3,gamma1]],[[2,3,-1*gamma2]]],
                          [[],[],[],[[1,1,Gamma13],[3,3,-1*gamma2],[2,2,gamma2]]]]

    # decoherence_matrix = [[[[0,0,-1*Gamma1]],[[0,1,-0.5*Gamma1]],[[0,2,-0.5*Gamma1]]],
    #                       [[],[[0,0,Gamma12],[1,1,-1*gamma1],[2,2,gamma1]],[[1,2,-1*gamma2]]],
    #                       [[],[],[[0,0,Gamma13],[2,2,-1*gamma2],[1,1,gamma2]]]]

    parameter = {'n':n,
                 'omega':[1E9+1E6,1E9,1E6,0],
                 'dipole':[[0,1000000,0,0],
                           [1000000,0,1000000,1000000],
                           [0,1000000,0,0],
                           [0,1000000,0,0]],
                 'nu':[1E9,1E6], # on resonence
                 'e_amp':[1,1],#amplitude of electric field
                 'level_group':[[0,1],[2],[3]],
                 'decoherence_matrix': decoherence_matrix
                 }
    



    filename = './test.txt'
#    system = System(n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    system = System(parameter)
    system.sweep(-1E8,1E8,400,'./test.txt')#TODO: add file name
    plot(n)
