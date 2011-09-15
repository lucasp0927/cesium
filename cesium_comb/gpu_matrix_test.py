#!/usr/bin/python2

from gpu_matrix import GPU_Matrix
import numpy as np
import time 
MATRIX_SIZE = 1024
gpu = GPU_Matrix(MATRIX_SIZE)
a = np.random.randn(MATRIX_SIZE, MATRIX_SIZE).astype(np.float64)
b = np.random.randn(MATRIX_SIZE, MATRIX_SIZE).astype(np.float64)
c = np.random.randn(MATRIX_SIZE, MATRIX_SIZE).astype(np.float64)
d = np.random.randn(MATRIX_SIZE, MATRIX_SIZE).astype(np.float64)
A = a+b*1j
B = c+d*1j
t = time.time()
C = gpu.matrix_mul(a,b,c,d)
print time.time() - t


