#!/usr/bin/python2
import numpy as np
from electricfield import Electricfield
from constant import *
from scipy.linalg.matfuncs import *

class Solver(object):
    """
    """

    def __init__(self,parameter,electric_field):
        """
        """
        print 'initializing...'
        self.EF = electric_field
        self.n = parameter['n'] #number of state
        self.dipole = parameter['dipole']
        self.omega = parameter['omega']
        self.decoherence_matrix = parameter['decoherence_matrix']
        self.N = self.n**2 #number of independent density matrix variable

        self.matrix_static = np.zeros([self.N,self.N],dtype = complex)#diagonal + decoherence
        self.matrix_electric = np.zeros([self.N,self.N],dtype = complex)
        print 'static part...'
        self.decoherence()
        self.diagonal_h()
        #print self.no_field_matrix()
        
    def index(self,i,j):
        if i==j:
            return (self.n + (self.n - i + 1)) * i - i
        elif j<i:
            #conjugate
            return (self.n + (self.n - j + 1)) * j - j + (i - j) * 2
        else:
            return -1-i*i+2*j+2*i*(-1+self.n)

    def reverse_index(self,n):
        for i in xrange(0,self.n):
            if n < 2*(self.n-i)-1:
                break
            else:
                n -= 2*(self.n-i)-1
        if n % 2 == 1:
            j = i+(n+1)/2
        elif n % 2 == 0:
            j = i
            i = j+n/2
        return i,j

    def decoherence(self):
        for i in range(0,self.n):
            for j in range(i,self.n):
                if i == j:
                    for iter in self.decoherence_matrix[i][j]:
                        k = iter[0]
                        l = iter[1]
                        s = iter[2]
                        self.matrix_static[self.index(i,j)][self.index(k,l)] += s
                else:
                    for iter in self.decoherence_matrix[i][j]:
                        k = iter[0]
                        l = iter[1]
                        s = iter[2]
                        self.matrix_static[self.index(i,j)][self.index(k,l)] += s
                        self.matrix_static[self.index(j,i)][self.index(l,k)] += s #need to add conjudate here if s is complex

    def diagonal_h(self):
        #[H,r]*(-i h)
        for i in range(0,self.n):
            for j in range(0,self.n):
                self.matrix_static[self.index(i,j)][self.index(i,j)]+=(self.omega[i]-self.omega[j])*(-1j)

    def no_field_matrix(self):
        return expm(np.matrix(self.matrix_static)*self.EF.zero_segment)
        
if __name__ == '__main__':
    pass
