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

    option = {1:'l=1,f=4',
              2:'l=1,f=3',
              3:'l=0,f=4',
              4:'l=0,f=3'
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
        plt.subplot(4,1,i)
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = option[i+1])
        plt.legend()
    plt.savefig("%s.png"%(file_out) )
    plt.show()

def plot_plt_d2(n):

    option = {1:'l=1,f=5',
              2:'l=1,f=4',
              3:'l=1,f=3',
              4:'l=1,f=2',
              5:'l=0,f=4',
              6:'l=0,f=3'
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
        plt.subplot(6,1,i)        
        plt.plot(data_sum[0,:],data_sum[i+1,:],label = option[i+1])
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
    file_name = str.split(file_in,'.')[0]
    #factor = [6,7,8,9,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]
    factor = [1000,2000]

    for f in factor:

        dictf = open(file_in,'r')
        parameter = eval(dictf.read())
        parameter['e_amp'] = list(parameter['e_amp'])
    #    parameter = pickle.load(dictf)
        dictf.close()

        file_out  = file_name + str(f) + "_freq.dat"
        e_amp =  parameter['e_amp']

        """
        calculate power
        """        
        epsilon_0 =  8.854187817620e-12
        C = 2.99792458e8
        amp = np.sqrt(2.162964e-2 * f / (epsilon_0 * C))
        power = epsilon_0 * C / 2.0 *  amp**2 + epsilon_0 * C / 2.0 *  amp**2
        print power
        parameter['e_amp'][0][0] = amp
        parameter['e_amp'][1][0] = amp
        parameter['power'] = power

        """
        initialize
        """
        system = System(parameter)
        parameter['sweep_profile'].append(file_out)

        # try:
        if parameter['d1'] == 1: #d1
            plot_current = plot_plt_d1
        elif parameter['d1'] == 2:#d2
            plot_current = plot_plt_d2
        # except KeyError:#d1 is not specific in parameter
                        #plot_current = plot_plt
        system.sweep(*parameter['sweep_profile'])# can parameter add after unpack array?
        #plot_current(parameter['n'])
