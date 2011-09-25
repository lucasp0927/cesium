#!/usr/bin/python2
from ctypes import *
import time
if __name__ == '__main__':
    libmatrixMul = CDLL("./obj/libmatrixMul.so")
    #libmatrixMul.runTest()
    tenDouble = c_double * 10
    data = tenDouble(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0)
    libmatrixMul.deviceVerify()
    libmatrixMul.doubleMatrixPrint(byref(data),10)
    for i in range(100):
        time.sleep(0.1)
        libmatrixMul.progressBar(100)
