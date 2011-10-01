#ifndef _MATRIXMUL_H_
#define _MATRIXMUL_H_
// Utilities and system includes
#include <stdio.h>
#include <cublas_v2.h>

//define forward
extern "C"
void doubleMatrixPrint(double* data,unsigned int size);
extern "C"
void returnMatrixPointer(cuDoubleComplex* data,char ri,int size,double* result);
extern "C"
cuDoubleComplex* complexIdentityMatrix(unsigned int size);
extern "C"
void solve(cuDoubleComplex* Hs, cuDoubleComplex* He,cuDoubleComplex* result, double* E_arr,double dt,int N, int finestep);
extern "C"  
void progressBar(unsigned int full);
extern "C"
cuDoubleComplex* complexMatrixCreate(double* real,double* imag,unsigned int size);
void complexTest();
extern "C"
void deviceVerify();
extern "C"
void runTest();
extern "C"
void computeGold(cuDoubleComplex* C, const cuDoubleComplex* A, const cuDoubleComplex* B, unsigned int hA, unsigned int wA, unsigned int wB);
extern "C"
void inline checkError(cublasStatus_t status, const char* msg);
extern "C"
void randomInit(cuDoubleComplex* data, int size);
extern "C"
void printMatrix(cuDoubleComplex* data, int size);
extern "C"
void printDiff(cuDoubleComplex *data1, cuDoubleComplex *data2, int width, int height, int iListLength, float fListTol);

#endif //#ifndef _MATRIXMUL_H_
