#!/usr/bin/python2
from solver import Solver
from electricfield import Electricfield
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    '''
    S = Solver()
    S.n = 10
    A = np.zeros((S.n,S.n))
    for i in xrange(S.n):
        for j in xrange(S.n):
            A[i][j] = S.index(i,j)
    print(np.matrix(A))
    for i in xrange(S.n):
        tmp = ''
        for j in xrange(S.n):
            tmp += ' %d,%d ' %(S.reverse_index(A[i][j]))
        print(tmp)
    '''
    para = {}
    para['Tr'] = 1.0/91.9262177e6
    para['mu_c'] = 351.72571850e12
    para['PSI'] = 3.0/4.0*2.0*np.pi
    para['E_0'] = 1.0
    para['tao'] = 3e-14
    EF = Electricfield(para)
    print 'period is ',EF.period
    print 'zero_segment ',EF.zero_segment
    
    file_in = 'setting/three_level.txt'
    dictf = open(file_in,'r')
    parameter = eval(dictf.read())
    dictf.close()    
    S = Solver(parameter,EF)
    #print np.matrix(S.matrix_static)
    # print np.matrix(S.matrix_electric)
    
    ### test_free decay
    init = np.zeros([9,1],dtype = complex)
    init[0][0] = 1.0+0.0j
    init = np.matrix(init)
    print init
    step = np.matrix(S.no_field_matrix())
    step_i = step
    a=[]
    b=[]
    c=[]
    for i in range(100):
        tmp = step_i * init
        a.append(tmp[0,0].real)
        b.append(tmp[5,0].real)
        c.append(tmp[8,0].real)
        step_i = step*step_i
    plt.figure()
    plt.plot(np.arange(100)*EF.zero_segment,a)
    plt.plot(np.arange(100)*EF.zero_segment,b)
    plt.plot(np.arange(100)*EF.zero_segment,c)    
    plt.show()

        
