#!/usr/bin/python2
from solver import Solver
from electricfield import Electricfield
import numpy as np
import matplotlib.pyplot as plt
from math import factorial

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
    
    """
    Set numpy printing setting
    """
    np.set_printoptions(precision = 2)   # print np array for human reading 
    # #np.set_printoptions(suppress=True)  
    np.set_printoptions(linewidth= 200)

    """
    stuff for initialization of solver object.
    """
    para = {}
    # para['Tr'] = 1.0/91.9262177e6
    # para['mu_c'] = 351.72571850e12
    # para['PSI'] = 3.0/4.0*2.0*np.pi
    # para['E_0'] = 1.0
    # para['tao'] = 3e-14

    
#    para['Tr'] = 1.0/91.9262177e6
    para['Tr'] = 1.0/91.9262177e6
    para['mu_c'] = 351.72571850e12
#    para['PSI'] = 3.0/4.0*2.0*np.pi
    para['PSI'] = 2.0*np.pi
    para['E_0'] = 0.0
    para['tao'] = 3e-14

    
    EF = Electricfield(para)

    print '---INFORMATION OF ELECTRICFIELD---'
    print 'period is ',EF.period
    print 'total time is ',EF.total_time
    print 'zero_segment ',EF.zero_segment
    print 'wave in packet', para['mu_c']*EF.time_no_field*2
    print '---END---\n\n'
    
#    file_in = 'setting/three_level.txt'
    file_in = 'setting/d2_2.txt'
    
    dictf = open(file_in,'r')
    parameter = eval(dictf.read())
    dictf.close()

    initial_state = np.zeros(parameter['n']**2,dtype = complex)
    initial_state[0] = 1.0+0.0j

    def new_test():
        S = Solver(parameter,EF,initial_state,1e-3)
        print 'going to simulate',S.total_period(),'total periods.'
#        S.main_control()
        S.main_control_matrix()
        print 'CPU'
        print S.matrix_no_field
        print 'GPU'        
        print S.period_matrix[0]
        print S.EF.zero_segment
        print S.EF.time_no_field *2.0
        print np.linalg.norm(S.period_matrix[0]- S.matrix_no_field)        
        
    new_test()


    
    def test1():
        S = Solver(parameter,EF)
        print S.build_matrix_dict(3)
        A = S.calculate_matrix_electric(0)
        print 'middle'
        print A
        print 'start'
        print S.no_field_matrix()
        print S.matrix_static

    def test2():
        S = Solver(parameter,EF)
        dt = 2.0*S.EF.time_no_field

        Hs =  np.matrix(S.matrix_static*dt)
        result = Hs
        print np.linalg.eig(Hs)
        for i in range(2,70):
            tmp = Hs
            for j in range(i-1):
                tmp = np.dot(Hs,tmp)
            result = np.matrix(result + tmp/float(factorial(i)))
        print result
        print "\n---CORRECT ANSWER---"
        print S.no_field_matrix()-np.identity(S.N,dtype = complex)
        

    # print S.build_matrix_dict(1)    
    #print np.matrix(S.matrix_static)
#    print np.matrix(S.matrix_total)
    
    # # ### test_free decay
    # init = np.zeros([9,1],dtype = complex)
    # init[0][0] = 1.0+0.0j
    # init = np.matrix(init)
    # print init
    # step = np.matrix(S.matrix_total)
    # step_i = step
    # a=[]
    # b=[]
    # c=[]
    # for i in range(100):
    #     tmp = step_i * init
    #     a.append(tmp[0,0].real)
    #     b.append(tmp[5,0].real)
    #     c.append(tmp[8,0].real)
    #     step_i = step*step_i
    # plt.figure()
    # plt.plot(np.arange(100)*EF.zero_segment,a)
    # plt.plot(np.arange(100)*EF.zero_segment,b)
    # plt.plot(np.arange(100)*EF.zero_segment,c)    
    # plt.show()

        
