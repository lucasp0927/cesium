#!/usr/bin/python2
from ctypes import *
import time
import numpy as np
import copy
if __name__ == '__main__':
    libmatrixMul = CDLL("./obj/libmatrixMul.so")
    #libmatrixMul.runTest()
    #    tenDouble = c_double * 10
    #    data = tenDouble(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0)
    A = np.array([[1.0+1.0j,2.0+2.0j,3.0+3.0j,4.0+4.0j,5.0+5.0j],[6.0+6.0j,7.0+7.0j,8.0+8.0j,9.0+9.0j,10.+10.0j]])
    A = A.ravel()
    libmatrixMul.deviceVerify()
    Ar = copy.copy(A.real)
    print A
    print Ar
    Ai = A.imag    
    libmatrixMul.doubleMatrixPrint(Ar.ctypes.data_as(POINTER(c_double)),10)
    for i in range(100):
        time.sleep(0.1)
        libmatrixMul.progressBar(100)
