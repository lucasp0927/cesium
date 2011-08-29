#!/usr/bin/python2
import sys
import os
from solver import Solver
import numpy as np
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


if __name__ ==  '__main__':
    file_in = sys.argv[1]
    file_out = sys.argv[2]

    #read parameter
    dictf = open(file_in,'r')
    parameter = eval(dictf.read())
    dictf.close()

#    system = System(parameter)
    """
    initiate Solver here
    """

    
