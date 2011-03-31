#!/usr/bin/python2
import sys
import os
from density_matrix import System

def plot():
    f=open('./tmp.gp','w')
    f.write('set terminal png\nset output \'graph.png\'\n')
    f.write('plot \'%s\' using 1:2 with lines, \'%s\' using 1:3 with lines, \'%s\' using 1:4 with lines \n' %(filename,filename,filename ) )
    f.close()
    os.system('gnuplot tmp.gp')
    os.remove('tmp.gp')

if __name__ ==  '__main__':
    n=3
    omega = [351E12,9E9,0]
    dipole=[[0,100000,100000],
            [100000,0,0],
            [100000,0,0]]
    #remember to /2
    nu=[351E12-9E9,351E12] # on resonence
    up = [0]
    low1 = [1]
    low2 = [2]
    #decoherence
    Gamma1 = 1000000
    Gamma12 = 500000
    Gamma13 = 500000
    gamma1 = 10000
    gamma2 = 10000
    filename = './test.txt'
    
    system = System(n,omega,dipole,nu,up,low1,low2,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    system.sweep(-1E7,1E7,400,'./test.txt')#TODO: add file name
    
    plot()
