#!/usr/bin/python2
import sys
import os
from density_matrix import System

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
    n=3
    omega = [105E10,9E9,0]
    dipole=[[0,1000000,1000000],
            [1000000,0,0],
            [1000000,0,0]]
    #remember to /2
    nu = [105E10-9E9,105E10] # on resonence
    e_amp = [1,1] 
    level_group = [[0],[1],[2]]
    #decoherence
    Gamma1 = 5000000
    Gamma12 = 2500000
    Gamma13 = 2500000
    gamma1 = 10000
    gamma2 = 10000
    '''
    n=4
    omega = [105E10+1E6,105E10,9E9,0]
    dipole=[[0,1000000,0,0],
            [1000000,0,1000000,1000000],
            [0,1000000,0,0],
            [0,1000000,0,0]]
    #remember to /2
    nu = [105E10-9E9,105E10] # on resonence
    e_amp = [1,1] 
    level_group = [[0,1],[2],[3]]
    #decoherence
    Gamma1 = 9000000
    Gamma12 = 3000000
    Gamma13 = 3000000
    gamma1 = 10000
    gamma2 = 10000
    '''
    filename = './test.txt'
    system = System(n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    system.sweep(-1E7,1E7,400,'./test.txt')#TODO: add file name
    plot(n)
