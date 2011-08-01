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

def plot_plt_d1(n):
    
    def l1f4():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=4')
    def l1f3():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=3')
    def l0f4():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=0,f=4')
    def l0f3():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=0,f=3')                    
    option = {1:l1f4,
              2:l1f3,
              3:l0f4,
              4:l0f3
        }

    line = count_line()
    f = open(file_out,'r')
    data = np.zeros([n+1,line])
    data_sum = np.zeros([5,line])    
    for aline in enumerate(f):
        data[:,aline[0]] = map(float,aline[1].split())
    f.close()
    for i in range(line):
        data_sum[0,i] = data[0,i]
        data_sum[1,i] = sum(data[1:10,i])
        data_sum[2,i] = sum(data[10:17,i])
        data_sum[3,i] = sum(data[17:26,i])
        data_sum[4,i] = sum(data[26:33,i])                        
    plt.figure()
    plt.xlabel("frequency (Hz)")
    for i in range(4):
        option[i+1]()
    plt.legend()
    plt.savefig("%s.png"%(file_out) )
    plt.show()

def plot_plt_d2(n):

    def l1f5():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=5')
    def l1f4():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=4')
    def l1f3():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=3')
    def l1f2():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=1,f=2')
    def l0f4():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=0,f=4')
    def l0f3():
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = 'l=0,f=3')                    
    option = {1:l1f5,
              2:l1f4,
              3:l1f3,
              4:l1f2,              
              5:l0f4,
              6:l0f3
        }    
    line = count_line()
    f = open(file_out,'r')
    data = np.zeros([n+1,line])
    data_sum = np.zeros([7,line])    
    for aline in enumerate(f):
        data[:,aline[0]] = map(float,aline[1].split())
    f.close()
    for i in range(line):
        data_sum[0,i] = data[0,i]
        data_sum[1,i] = sum(data[1:12,i])
        data_sum[2,i] = sum(data[12:21,i])
        data_sum[3,i] = sum(data[21:28,i])
        data_sum[4,i] = sum(data[28:33,i])
        data_sum[5,i] = sum(data[33:42,i])
        data_sum[6,i] = sum(data[42:49,i])                                
    plt.figure()
    plt.xlabel("frequency (Hz)")
    for i in range(6):
        option[i+1]()
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
    plot_plt_d2(parameter['n'])

