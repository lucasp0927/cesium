#!/usr/bin/python2
from density_matrix import System
import sys
import os
if __name__ ==  '__main__':
    n=3
    omega = [351E12,9E9,0]
    dipole=[[0,1000000,1000000],
            [1000000,0,0],
            [1000000,0,0]]
    #remember to /2
    nu=[351E12-9E9,351E12] # on resonence
    up = [0]
    low1 = [1]
    low2 = [2]
    #decoherence
    Gamma1 = 5000000
    Gamma12 = 2500000
    Gamma13 = 2500000
    gamma1 = 10000
    gamma2 = 10000
    filename = './test.txt'
    system = System(n,omega,dipole,nu,up,low1,low2,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    system.sweep(-1E7,1E7,1000,'./test.txt')#TODO: add file name

    f=open('./tmp.gp','w')
    f.write('set terminal png\nset output \'graph.png\'\n')
    f.write('plot \'%s\' using 1:2 with lines, \'%s\' using 1:3 with lines, \'%s\' using 1:4 with lines \n' %(filename,filename,filename ) )
    f.close()
    os.system('gnuplot tmp.gp')
    os.remove('tmp.gp')
