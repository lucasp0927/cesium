#!/usr/bin/python2
import sys
import os
from density_matrix import System
import numpy as np
import pickle #pickle is not very safe
import matplotlib.pyplot as plt

def count_line():
    """
    count lines in file_out
    """
    f = open(file_out,'r')
    line = 0
    for aline in f:
        line += 1
    f.close()
    return line

def plot_plt(n):
    line = count_line()
    f = open(file_out,'r')
    data = np.zeros([n+1,line])
    for aline in enumerate(f):
        data[:,aline[0]] = map(float,aline[1].split())
    f.close()
    plt.figure()
    plt.xlabel("frequency (Hz)")
    for i in range(n):
        plt.plot(data[0,:],data[i+1,:],label = str(i))
    plt.legend()
    plt.savefig("%s.png"%(file_out) )
    plt.show()
    

def plot_gnuplot(n):
    """
    Plot using gnuplot. 
    """
    f=open('./tmp.gp','w')
    f.write('set terminal png\nset output \'%s.png\'\n'%(file_out) )
    tmp_str = 'plot \'%s\' using 1:2 with lines'%file_out
    for i in range(n-1):
        tmp_str += ', \'%s\' using 1:%d with lines'%(file_out,i+3)
    tmp_str += '\n'
    f.write(tmp_str)
    f.close()
    os.system('gnuplot tmp.gp')
    os.remove('tmp.gp')

if __name__ ==  '__main__':
    file_in = sys.argv[1]
    file_out = sys.argv[2]
    
    # dictf = open('three_level','w')
    # pickle.dump(parameter,dictf)
    # dictf.close()

    dictf = open(file_in,'r')
    parameter = eval(dictf.read())
#    parameter = pickle.load(dictf)
    dictf.close()
    
    # txtf = open(file_in+'.txt','w')
    # txtf.write(str(parameter))
    # txtf.close()

    system = System(parameter)
    parameter['sweep_profile'].append(file_out)
    system.sweep(*parameter['sweep_profile'])# can parameter add after unpack array?
    # plot_gnuplot(parameter['n'])
    plot_plt(parameter['n'])
