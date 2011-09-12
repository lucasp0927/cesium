#!/usr/bin/python2
from __future__ import division
import numpy as np
from electricfield import Electricfield
from scipy.weave import converters
import scipy

class DESolver(object):
    """
    """
    def __init__(self,efield,step):
        self.EF = efield
        self.step = step
        self.fine_step = 2*step
        self.t_arr = np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,self.fine_step+1)
        self.dt_fine = 2.0*self.EF.time_no_field/(self.fine_step)
        self.E_arr = []
        for period in range(self.EF.period):
            self.E_arr.append([self.EF.comb_field(self.t_arr[i],period) for i in xrange(self.fine_step+1)])
