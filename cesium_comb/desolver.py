#!/usr/bin/python2
from __future__ import division
import numpy as np
from electricfield import Electricfield
from scipy.weave import converters
import scipy

class DESolver(object):
    """
    """
    def __init__(self,efield,step,matrix_static,matrix_electric):
        self.EF = efield
        self.step = step
        self.matrix_static = np.array(matrix_static)
        self.matrix_electric = np.array(matrix_electric)
        self.fine_step = 2*step
        self.t_arr = np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,self.fine_step+1)
        self.dt_fine = 2.0*self.EF.time_no_field/(self.fine_step)
        self.dt = 2.0*self.EF.time_no_field/(self.step)        
        self.E_arr = []
        for period in range(self.EF.period):
            self.E_arr.append([self.EF.comb_field(self.t_arr[i],period) for i in xrange(self.fine_step+1)])

    def solve(self, initial_state, period):
        # t_d =  np.dot((self.matrix_electric*E+self.matrix_static),state)
        current = initial_state
        for i in xrange(0,self.step,2):
            k1 =  np.dot((self.matrix_electric*self.E_arr[period][i]+self.matrix_static),current)*self.dt
            k2 =  np.dot((self.matrix_electric*self.E_arr[period][i+1]+self.matrix_static),current+k1*0.5)*self.dt
            k3 =  np.dot((self.matrix_electric*self.E_arr[period][i+1]+self.matrix_static),current+k2*0.5)*self.dt
            k4 =  np.dot((self.matrix_electric*self.E_arr[period][i+2]+self.matrix_static),current+k3)*self.dt                         
            current = current + 1.0/6.0* (k1+2.0*k2+2.0*k3+k4)
        return current
