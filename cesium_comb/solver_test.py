#!/usr/bin/python2
from solver import Solver
import numpy as np

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
    file_in = 'setting/three_level.txt'
    dictf = open(file_in,'r')
    parameter = eval(dictf.read())
    dictf.close()    
    S = Solver(parameter)
    print np.matrix(S.matrix_static)
    print np.matrix(S.matrix_electric)
