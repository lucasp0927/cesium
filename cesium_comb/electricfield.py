#!/usr/bin/python2
import numpy as np
import matplotlib.pyplot as plt

def comb_field(t,n):
    mu_c = 1e7 #carier frequency
    PSI = 1.0/2.0*np.pi #phase difference
    #envolope
    return np.cos(2*np.pi*mu_c*t - n*PSI)*envolope(t)
    
def envolope(t):
    E_0 = 1.0
    tao = 2e-6
    return E_0*np.exp(-1.0*np.power(t/tao,2))
    
if __name__ == '__main__':
    t = np.linspace(-10e-6,10e-6,1000)
    x = comb_field(t,1)
    plt.figure()
    plt.plot(t,x)
    plt.show()
