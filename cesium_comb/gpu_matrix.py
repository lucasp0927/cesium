#!/usr/bin/python2
import numpy as np
from numpy import linalg as la
from pycuda import driver, compiler, gpuarray, tools
import time
# -- initialize the device
import pycuda.autoinit

class GPU_Matrix(object):
    """ 
    Multiples two square matrices together using multiple blocks and shared memory. 
    Each thread block is assigned a "tile" of the resulting matrix and is responsible
    for generating the elements in that tile.  Each thread in a block computes one element 
    of the tile.
    """
    
    def __init__(self,MATRIX_SIZE):
        """
            
        """
        self.MATRIX_SIZE = MATRIX_SIZE
        # define size of blocks and tiles sub-matrix 
        # (we assume that the block size is same as tile size)
        self.TILE_SIZE = 2
        self.BLOCK_SIZE = self.TILE_SIZE        
        kernel_code_template = """
        __global__ void MatrixMulKernel(double *A, double *B, double *C, double *D, double *E, double *F)
        {

          const uint wA = %(MATRIX_SIZE)s;
          const uint wB = %(MATRIX_SIZE)s;  

          // Block index
          const uint bx = blockIdx.x;
          const uint by = blockIdx.y;

          // Thread index
          const uint tx = threadIdx.x;
          const uint ty = threadIdx.y;

          // Index of the first sub-matrix of A processed by the block
          const uint aBegin = wA * %(BLOCK_SIZE)s * by;
          // Index of the last sub-matrix of A processed by the block
          const uint aEnd = aBegin + wA - 1;
          // Step size used to iterate through the sub-matrices of A
          const uint aStep = %(BLOCK_SIZE)s;

          // Index of the first sub-matrix of B processed by the block
          const uint bBegin = %(BLOCK_SIZE)s * bx;
          // Step size used to iterate through the sub-matrices of B
          const uint bStep = %(BLOCK_SIZE)s * wB;

          // The element of the block sub-matrix that is computed
          // by the thread
          double Esub = 0;
          double Fsub = 0;  
          // Loop over all the sub-matrices of A and B required to
          // compute the block sub-matrix
          for (int a = aBegin, b = bBegin;
               a <= aEnd;
               a += aStep, b += bStep) 
            {
              // Shared memory for the sub-matrix of A
              __shared__ double As[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];
              __shared__ double Cs[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];      
              // Shared memory for the sub-matrix of B
              __shared__ double Bs[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];
              __shared__ double Ds[%(BLOCK_SIZE)s][%(BLOCK_SIZE)s];      

              // Load the matrices from global memory to shared memory;
              // each thread loads one element of each matrix
              As[ty][tx] = A[a + wA * ty + tx];
              Bs[ty][tx] = B[a + wA * ty + tx];
              Cs[ty][tx] = C[b + wB * ty + tx];
              Ds[ty][tx] = D[b + wB * ty + tx];      
              // Synchronize to make sure the matrices are loaded
              __syncthreads();

              // Multiply the two matrices together;
              // each thread computes one element
              // of the block sub-matrix
              for (int k = 0; k < %(BLOCK_SIZE)s; ++k)
              {
                Esub += As[ty][k] * Cs[k][tx] - Bs[ty][k] * Ds[k][tx];
                Fsub += As[ty][k] * Ds[k][tx] + Bs[ty][k] * Cs[k][tx];        
              }
              // Synchronize to make sure that the preceding
              // computation is done before loading two new
              // sub-matrices of A and B in the next iteration
              __syncthreads();
            }

          // Write the block sub-matrix to global memory;
          // each thread writes one element
          const uint c = wB * %(BLOCK_SIZE)s * by + %(BLOCK_SIZE)s * bx;
          E[c + wB * ty + tx] = Esub;
          F[c + wB * ty + tx] = Fsub;
        }
        
        """
        # get the kernel code from the template 
        # by specifying the constants MATRIX_SIZE and BLOCK_SIZE
        kernel_code = kernel_code_template % { 
            'MATRIX_SIZE': self.MATRIX_SIZE,
            'BLOCK_SIZE': self.BLOCK_SIZE,
            }

        # compile the kernel code
        mod = compiler.SourceModule(kernel_code)

        # get the kernel function from the compiled module
        self.matrixmul = mod.get_function("MatrixMulKernel")


    def matrix_mul(self,a_gpu,b_gpu,c_gpu,d_gpu,e_gpu,f_gpu):
        # a_gpu = gpuarray.to_gpu(Ar)
        # b_gpu = gpuarray.to_gpu(Ai)
        # c_gpu = gpuarray.to_gpu(Br) 
        # d_gpu = gpuarray.to_gpu(Bi)
        # e_gpu = gpuarray.empty((self.MATRIX_SIZE, self.MATRIX_SIZE), np.float64)
        # f_gpu = gpuarray.empty((self.MATRIX_SIZE, self.MATRIX_SIZE), np.float64)

        self.matrixmul(
            # inputs
            a_gpu, b_gpu, c_gpu ,d_gpu,
            # output
            e_gpu, f_gpu, 
            # grid of multiple blocks
            grid = (self.MATRIX_SIZE / self.TILE_SIZE, self.MATRIX_SIZE / self.TILE_SIZE),
            # block of multiple threads
            block = (self.TILE_SIZE, self.TILE_SIZE, 1), 
            )
#        return e_gpu.get(Cr),f_gpu.get(Ci)
