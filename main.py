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
    omega = [5E10,9E9,0]
    dipole=[[0,100000,100000],
            [100000,0,0],
            [100000,0,0]]
    #remember to /2
    nu = [5E10-9E9,5E10] # on resonence
    e_amp = [10,10] 
    level_group = [[0],[1],[2]]
    #decoherence
    Gamma1 = 1000000
    Gamma12 = 500000
    Gamma13 = 500000
    gamma1 = 1000
    gamma2 = 1000    
    filename = './test.txt'
    system = System(n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    system.sweep(-1E7,1E7,400,'./test.txt')#TODO: add file name
    plot(n)
