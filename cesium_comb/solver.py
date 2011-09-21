#!/usr/bin/python2
from __future__ import division
import numpy as np
from electricfield import Electricfield
#from desolver import DESolver
from desolver_matrix import DESolver_matrix
from constant import *
from scipy.linalg.matfuncs import *
from scipy.weave import converters
from math import factorial
import scipy

class Solver(object):
    """
    """

    def __init__(self,parameter,electric_field,initial_state,end_time):
        """
        """
        self.end_time = end_time
        print 'SOLVER INITIALIZE'
        print '   initializing...'
        self.EF = electric_field
        self.n = parameter['n'] #number of state
        self.dipole = parameter['dipole']
        self.omega = parameter['omega']
        self.decoherence_matrix = parameter['decoherence_matrix']
        self.N = self.n**2 #number of independent density matrix variable

        self.matrix_static = np.zeros([self.N,self.N],dtype = complex)#diagonal + decoherence
        self.matrix_electric = np.zeros([self.N,self.N],dtype = complex)
        self.matrix_total = np.identity(self.N,dtype = complex)
        self.period_matrix = []
        if initial_state.size != self.N:
            raise ValueError('wrong initial state size')
        #may add normalize check here
        self.initial_state = np.matrix(initial_state).transpose()
        self.current_state = self.initial_state

        print '   static part...'
        self.decoherence()
        self.diagonal_h()
        print '   dipole....'
        self.matrix_dipole()
        print '   calculate no field matrix...'
        self.matrix_no_field = self.no_field_matrix()
        # print '   initialize DE solver...'
        # self.desolver = DESolver(efield = self.EF,step = 10000,matrix_static = self.matrix_static,matrix_electric = self.matrix_electric)
        print '   initialize DE solver (matrix version)...'
        self.desolver_matrix = DESolver_matrix(efield = self.EF,step = 10000,matrix_static = self.matrix_static,matrix_electric = self.matrix_electric,MATRIX_SIZE = self.N,order = 5)
        print '\n\n'

    def total_period(self):
        return int(np.round(self.end_time/self.EF.total_time))


    def time_differ(self,E,state):
        return np.dot((matrix_electric*E+matrix_static),state)

    def main_control_matrix(self):
        #prepare period matrix
        print 'prepare period matrix'
        for n in range(self.EF.period):
            print 'calculate period %d'%(n)
            self.period_matrix.append(self.desolver_matrix.solve(n,self.N))
        #calculate

    def main_control(self):
        time = 0.0
        for counter in xrange(self.total_period()):
            for period in xrange(self.EF.period):
                self.calculate_period(period)
            time += self.EF.total_time
            print counter
        print 'total time is ',time

    def calculate_period(self,n):
        self.current_state = np.dot(self.matrix_no_field,self.current_state)
        self.current_state = self.desolver.solve(self.current_state,n)
        self.current_state = np.dot(self.matrix_no_field,self.current_state)




    def calculate_period_old(self):
        matrix_nofield = self.no_field_matrix()
        ##DEBUG
        #matrix_nofield = np.matrix(np.eye(self.N,dtype = complex))
        for n in range(self.EF.period):
            print 'calculate period %d'%(n)
            self.matrix_total = matrix_nofield*self.matrix_total
            self.matrix_total = self.calculate_matrix_electric(n)*self.matrix_total
            self.matrix_total = matrix_nofield*self.matrix_total

    def calculate_matrix_electric(self,n):
        num = 10000
        result = np.matrix(np.zeros([self.N,self.N],dtype = complex))
        identity_matrix =np.matrix(np.identity(self.N,dtype = complex))
        t_arr = np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,num)
        dt = 2.0*self.EF.time_no_field/(num-1)
        E_arr = np.zeros(num)
        for i in xrange(num):
            E_arr[i] = self.EF.comb_field(t_arr[i],n)
        Hs = self.matrix_static * dt
        He = self.matrix_electric *dt
        print self.sum_matrix(Hs)
        print self.sum_matrix(He)
        print Hs
        for pass_n in range(1,2):
#            tmp = result = np.matrix(np.zeros([self.N,self.N],dtype = complex))
            dict = self.build_matrix_dict(pass_n)

            ###This is the slowest part
            c_tmp = 0
            E_x = np.zeros(pass_n)

#             for x in self.perm(pass_n,0,num):
#                 if c_tmp != x[0]:
#                     print x[0]
#                 c_tmp = x[0]

#                 for i in range(pass_n):
#                     E_x[i] = E_arr[x[i]]

#                 # patt = dict['pattern']
#                 # coef = dict['coef']

#                 # code="""
#                 #      double coef_tmp;
#                 #      for (int i = 0; i<2^pass_n; i++)
#                 #      {
#                 #         coef_tmp = 1.0;
#                 #         for (int p_id = 0; p_id<pass_n; p_id++)
#                 #         {
#                 #         if (patt(i,p_id) == 1)
#                 #            coef_tmp *= E_x(p_id);
#                 #         }
#                 #         coef(i) += coef_tmp;
#                 #      }
#                 #      return_val = coef;
#                 #      """
#                 # coef = scipy.weave.inline(code,
#                 #                           ['patt','E_x','pass_n','coef'],
#                 #                           type_converters=converters.blitz,
#                 #                           compiler = 'gcc')

#                 for index in range(2**pass_n):
#  #                    coef_tmp = np.prod(E_x**dict['pattern'][index])
#                      patt = dict['pattern'][index]
#                      code="""
#                           double coef_tmp;
#                           coef_tmp = 1.0;
#                           for (int p_id = 0; p_id<pass_n; p_id++)
#                           {
#                             if (patt(p_id) == 1)
#                                coef_tmp *= E_x(p_id);
#                           }
#                           return_val = coef_tmp;
#                           """
#                      dict['coef'][index]+=scipy.weave.inline(code,
#                                                        ['patt','E_x','pass_n'],
#                                                        type_converters=converters.blitz,
#                                                        compiler = 'gcc')
# #                    dict['coef'][index] += 1
            dict['coef'][0] = (factorial(num) / (factorial(pass_n)*factorial(num-pass_n)))
            print dict


            for index in range(2**pass_n):
                tmp = np.matrix(np.identity(self.N,dtype = complex))
                print dict['pattern'][index]
                for p_id in range(pass_n):
                    if dict['pattern'][index][p_id] == 1:
                        tmp = np.dot(He,tmp)
                    else:
                        tmp = np.dot(Hs,tmp)
                tmp_0 = tmp*dict['coef'][index]
                print 'average :', self.sum_matrix(tmp_0)
                result = tmp_0 + result
        result = np.matrix(np.identity(self.N,dtype = complex))+result
        print "test"
        print np.dot(Hs,np.dot(Hs,Hs))
        return result

        # simple version, slow. Not accurate.
        # m1 = identity_matrix + self.matrix_static * dt
        # print len(m1)
        # for t in t_arr:
        #     print(t)
        #     E = self.EF.comb_field(t,n) * dt
        #     tmp = np.dot((m1+ self.matrix_electric * E),tmp)
        # return tmp

    def sum_matrix(self,A):
        sum = np.sum(np.abs(A))
        return sum/(self.N**2)

    def build_matrix_dict(self,n):
        if n > 20:
            print 'n to big'
            raise ValueError
        dict = {}
        dict['pattern'] = np.zeros([2**n,n])
#        dict['coef'] = [0 for i in range(2**n)]
        dict['coef'] = np.zeros(2**n)
        for i in range(2**n):
#            dict['pattern'][i] = ("{:0"+str(n)+"b}").format(i)#todo:find better way
            tmp = list(("{:0"+str(n)+"b}").format(i))
            for j in range(len(tmp)):
                tmp[j] = int(tmp[j])
            tmp = np.array(tmp)
            dict['pattern'][i] = tmp
        return dict


    def perm(self,n,start,end):
        #return array of permutation of start to end - 1
        if n == 1:
            for i in xrange(start,end):
                yield [i]
        else:
            result = []
            for i in xrange(start,end):
                result = [i]
                for j in self.perm(n-1,i+1,end):
                    result.extend(j)
                    yield(result)
                    result = [i]

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
        #[H,r]*(-i/h)
        for i in range(0,self.n):
            for j in range(0,self.n):
                self.matrix_static[self.index(i,j)][self.index(i,j)]+=(self.omega[i]-self.omega[j])*(-1j)

    def matrix_dipole(self):
        """
        need to time the electric field to get real time differential.
        """
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    self.matrix_electric[self.index(i,j)][self.index(k,j)]+=self.dipole_h(i,k)*(-1j/HBAR())
                    self.matrix_electric[self.index(i,j)][self.index(i,k)]-=self.dipole_h(k,j)*(-1j/HBAR())

    def dipole_h(self,i,j):
        if i <= j:
            return self.dipole[0][i][j]
        else:
            return self.dipole[0][j][i]

    def no_field_matrix(self):
            return expm(np.matrix(self.matrix_static)*self.EF.zero_segment)

if __name__ == '__main__':
    pass
