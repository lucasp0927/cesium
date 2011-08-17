#!/usr/bin/python2
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


class Electricfield():
    """
    """
    
    def __init__(self,para):
        """
        """
        [self.Tr,#repetition rate
        self.mu_c, #carrier freq
        self.PSI, #phase difference
        self.E_0, #electric field
        self.tao]=[para['Tr'],
                   para['mu_c'],
                   para['PSI'],
                   para['E_0'],
                   para['tao']]
        self.period = self.calculate_period()
        
    def comb_field(self,t,n):
        mu_c = 351.72571850e12 #carier frequency 852 nm
        PSI = 1.0/2.0*np.pi #phase difference
        return np.cos(2*np.pi*mu_c*t - n*PSI)*envolope(t)

    def envolope(self,t):
        #envolope
        E_0 = 1.0
        tao = 3e-14
        return E_0*np.exp(-1.0*np.power(t/tao,2))

    def calculate_period(self):
        """
        Return after how many pulse will PSI add up to 2 pi's multipier.
        The donominator of i/j*2pi can be up to 300.
        """
        # calculate the greatest common divisor
        def gcd(first, second):
            counter = 0
            if (first < second):
                first, second = second, first
            while np.abs(second) > 1e-10:
                counter += 1
                first, second =  first % second , second
                if (first < second):
                    first, second = second, first
                if counter > 50:
                    raise ValueError('PSI is not good. %e' %(second) )
            return first

        # calculate the least common multiple
        def lcm(first, second):
                gcd_num = gcd(first, second)
                return first * second / gcd_num

        a = self.PSI/np.pi/2
        gcd_num = gcd(a,1)
        period = int(np.round(lcm(a/gcd_num,1/gcd_num)*gcd_num/a))
        #check answer
        m = np.round(self.PSI*period/(2.0*np.pi))
        error = self.PSI*period - m*2.0*np.pi
        if  error > 1e-11:
            raise ValuError('Period is incorrect %e' %(second) )
        return period
        
if __name__ == '__main__':
    pass

