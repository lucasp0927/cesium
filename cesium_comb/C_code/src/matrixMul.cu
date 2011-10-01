#include "matrixMul.h"
#define SCALE(vector,num)\
{const cuDoubleComplex alpha = make_cuDoubleComplex(num,0.0);  \
  cublasZscal(handle,N*N,&alpha,vector,1);\
}\

#define MUL(x,y,result)                               \
  {             \
  const cuDoubleComplex alpha = make_cuDoubleComplex(1.0,0.0);\
  const cuDoubleComplex beta = make_cuDoubleComplex(0.0,0.0);\
  cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha, x, N, y, N, &beta, result, N);\
  }             \


#define COPY(from,to)  cublasZcopy(handle,N*N,from,1,to,1)
// y = ax + y
#define AXPY(y,x,a)\
{\
  const cuDoubleComplex alpha = make_cuDoubleComplex(a,0.0);\
  cublasZaxpy (handle,N*N, &alpha , x , 1 , y , 1);             \
}\



void progressBar(unsigned int full)
{
  static unsigned int now = 1;
  printf("%d\% \r",now*100/full);
  fflush(stdout);
  now++;

}

void doubleMatrixPrint(double* data,unsigned int size)
{
  for (int i = 0; i < size; ++i)
    {
      printf("%d:\t %f \n",i,data[i]);
    }
}

cuDoubleComplex* complexMatrixCreate(double* real,double* imag,unsigned int size)
{
    cuDoubleComplex* A  = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex)*size);
    for (int i = 0; i < size; ++i)
      {
        A[i] = make_cuDoubleComplex(real[i],imag[i]);
      }
    return A;
}

cuDoubleComplex* complexIdentityMatrix(unsigned int size)
{
  cuDoubleComplex* A  = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex)*(size*size));
  for (int i = 0; i < size*size; ++i)
    {
      A[i] = make_cuDoubleComplex(0.0,0.0);
    }
  for (int i = 0; i < size; ++i)
    {
      A[i*size+i] = make_cuDoubleComplex(1.0,0.0);
    }
  return A;
}

void deviceVerify ()
{
  int devID;
  cudaDeviceProp props;
  // get number of SMs on this GPU
  cudaGetDevice(&devID);
  cudaGetDeviceProperties(&props, devID);
  printf("Device %d: \"%s\" with Compute %d.%d capability\n", devID, props.name, props.major, props.minor);
}

void complexTest()
{
  /*
    http://bccd-ng.cluster.earlham.edu/svn/bccd-ng/branches/lemanal-devel/trees/software/cluster/software/cuda/include/cuComplex.h
    usage of cuda complex.
    cuCreal
    cuCimag
    cuCadd
    cuCsub
    cuCmul
    cuCdiv
    cuCabs
   */
  cuDoubleComplex* A  = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex)*10);
  cuDoubleComplex* B  = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex)*10);
  randomInit(A,10);
  randomInit(B,10);
  printMatrix(A,10);
  printMatrix(B,10);
  printDiff(A,B,3,3,100,1e-6);
  free(A);
  free(B);
}


//cpu reference
void computeGold(cuDoubleComplex* C, const cuDoubleComplex* A, const cuDoubleComplex* B, unsigned int hA, unsigned int wA, unsigned int wB)
{
  for (unsigned int i = 0; i < hA; ++i)
    for (unsigned int j = 0; j < wB; ++j) {
      cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
      for (unsigned int k = 0; k < wA; ++k) {
        cuDoubleComplex a = A[i * wA + k];
        cuDoubleComplex b = B[k * wB + j];
        sum = cuCadd(sum,cuCmul(a,b));
      }
      C[i * wB + j] = sum;
    }
}


void inline checkError(cublasStatus_t status, const char* msg)
{
  if(status != CUBLAS_STATUS_SUCCESS){
    printf(msg);
    exit(-1);
  }
}

void randomInit(cuDoubleComplex* data, int size)
{
  // set seed for rand()
  srand(2006);
  for (int i = 0; i < size; ++i)
    {
      cuDoubleComplex tmp;
      tmp = make_cuDoubleComplex(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
      data[i] = tmp;
    }

}

void returnMatrixPointer(cuDoubleComplex* data,char ri,int size, double* result)
{
  for (int i = 0; i < size; ++i)
    {
     if (ri == 'r')
       {
         result[i] = cuCreal(data[i]);
       }
     else
       {
         result[i] = cuCimag(data[i]);
       }
    }
}

void printMatrix(cuDoubleComplex* data, int size)
{
    for (int i = 0; i < size; ++i)
    {
      printf("%d:  %f + %f j \n",i,cuCreal(data[i]),cuCimag(data[i]));
    }
    printf("\n \n");

}

void printDiff(cuDoubleComplex *data1, cuDoubleComplex *data2, int width, int height, int iListLength, float fListTol)
{
  printf("Listing first %d Differences > %.6f...\n", iListLength, fListTol);
  int i,j,k;
  int error_count=0;
  for (j = 0; j < height; j++)
    {
      if (error_count < iListLength)
        {
          printf("\n  Row %d:\n", j);
        }
      for (i = 0; i < width; i++)
        {
          k = j * width + i;
          float fDiff = cuCabs(cuCsub(data1[k],data2[k]));
          if (fDiff > fListTol)
            {
              if (error_count < iListLength)
                {
                  printf("    Loc(%d,%d)\tCPU=%.5f+%.5f j\tGPU=%.5f+%.5f j\tDiff=%.6f\n", i, j, cuCreal(data1[k]),cuCimag(data1[k]),cuCreal(data2[k]),cuCimag(data2[k]), fDiff);
                }
              error_count++;
            }
        }
    }
  printf(" \n  Total Errors = %d\n\n", error_count);
}

void solve(cuDoubleComplex* Hs, cuDoubleComplex* He,cuDoubleComplex* result, double* E_arr,double dt,int N, int finestep)
{
  printf ("\nstart gpu solving...\n");
  printf ("matrix size:%d\n",N*N);
  printf ("step size:%d\n\n",finestep/6);

  unsigned int mem_size = sizeof(cuDoubleComplex) * (N*N);

  // allocate device memory
  cuDoubleComplex* d_Hs, *d_He, *d_result,* d_result_tmp, *d_k1, *d_k2, *d_k3, *d_k4, *d_k5, *d_k6, *d_tmp,*d_tmp2, *d_I;
  cudaMalloc((void**) &d_Hs, mem_size);
  cudaMalloc((void**) &d_He, mem_size);
  cudaMalloc((void**) &d_result, mem_size);
  cudaMalloc((void**) &d_result_tmp, mem_size);  
  cudaMalloc((void**) &d_k1, mem_size);
  cudaMalloc((void**) &d_k2, mem_size);
  cudaMalloc((void**) &d_k3, mem_size);
  cudaMalloc((void**) &d_k4, mem_size);
  cudaMalloc((void**) &d_k5, mem_size);
  cudaMalloc((void**) &d_k6, mem_size);
  cudaMalloc((void**) &d_tmp, mem_size);
  cudaMalloc((void**) &d_tmp2, mem_size);
  cudaMalloc((void**) &d_I, mem_size);

  cuDoubleComplex* I = complexIdentityMatrix(N);

  // copy host memory to device
  cudaMemcpy(d_Hs, Hs, mem_size, cudaMemcpyHostToDevice) ;
  cudaMemcpy(d_He, He, mem_size, cudaMemcpyHostToDevice) ;
  cudaMemcpy(d_result, result, mem_size, cudaMemcpyHostToDevice) ;
  cudaMemcpy(d_I, I, mem_size, cudaMemcpyHostToDevice) ;

  //cublas_v2
  cublasHandle_t handle;
  checkError(cublasCreate(&handle), "cublasCreate() error!\n");

  //test: scale result  y = ax + y

  /* COPY(d_result,d_tmp); */
  /* AXPY(d_tmp,d_result,2.0); */
  /* MUL(d_result,d_Hs,d_tmp); */
  /* COPY(d_tmp,d_result); */

  for (int i = 0; i < finestep; i=i+6)
    {
      progressBar(finestep/6);
      /* k1 */
      COPY(d_He,d_k1);
      SCALE(d_k1,E_arr[i]);
      AXPY(d_k1,d_Hs,1.0);
      /* tmp */
      COPY(d_He,d_tmp);
      SCALE(d_tmp,E_arr[i+1]);
      AXPY(d_tmp,d_Hs,1.0);
      /* tmp2 */
      COPY(d_I,d_tmp2);
      AXPY(d_tmp2,d_k1,0.25*dt);
      /* k2 */
      MUL(d_tmp,d_tmp2,d_k2);
      /* tmp */
      COPY(d_He,d_tmp);
      SCALE(d_tmp,E_arr[i+2]);
      AXPY(d_tmp,d_Hs,1.0);
      /* tmp2 */
      COPY(d_I,d_tmp2);
      AXPY(d_tmp2,d_k1,3.0/32.0*dt);
      AXPY(d_tmp2,d_k2,9.0/32.0*dt);
      /* k3 */
      MUL(d_tmp,d_tmp2,d_k3);
      /* tmp */
      COPY(d_He,d_tmp);
      SCALE(d_tmp,E_arr[i+3]);
      AXPY(d_tmp,d_Hs,1.0);
      /* tmp2 */
      COPY(d_I,d_tmp2);
      AXPY(d_tmp2,d_k1,1932.0/2197.0*dt);
      AXPY(d_tmp2,d_k2,-1.0*7200.0/2197.0*dt);
      AXPY(d_tmp2,d_k3,7296.0/2197.0*dt);
      /* k4 */
      MUL(d_tmp,d_tmp2,d_k4);
      /* tmp */
      COPY(d_He,d_tmp);
      SCALE(d_tmp,E_arr[i+4]);
      AXPY(d_tmp,d_Hs,1.0);
      /* tmp2 */
      COPY(d_I,d_tmp2);
      AXPY(d_tmp2,d_k1,439.0/216.0*dt);
      AXPY(d_tmp2,d_k2,-1.0*8.0*dt);
      AXPY(d_tmp2,d_k3,3680.0/513.0*dt);
      AXPY(d_tmp2,d_k4,-1.0*845.0/4104.0*dt);
      /* k5 */
      MUL(d_tmp,d_tmp2,d_k5);

      /* tmp */
      COPY(d_He,d_tmp);
      SCALE(d_tmp,E_arr[i+5]);
      AXPY(d_tmp,d_Hs,1.0);
      /* tmp2 */
      COPY(d_I,d_tmp2);
      AXPY(d_tmp2,d_k1,-1.0*8.0/27.0*dt);
      AXPY(d_tmp2,d_k2,-1.0*2.0*dt);
      AXPY(d_tmp2,d_k3,-1.0*3544.0/2565.0*dt);
      AXPY(d_tmp2,d_k4,1859.0/4104.0*dt);
      AXPY(d_tmp2,d_k5,-1.0*11.0/40.0*dt);      
      /* k6 */
      MUL(d_tmp,d_tmp2,d_k6);
      /* tmp */
      COPY(d_I,d_tmp);
      AXPY(d_tmp,d_k1,16.0/135.0*dt);
      AXPY(d_tmp,d_k3,6656.0/12825.0*dt);
      AXPY(d_tmp,d_k4,28561.0/56430.0*dt);
      AXPY(d_tmp,d_k5,-1.0*9.0/50.0*dt);
      AXPY(d_tmp,d_k6,2.0/55.0*dt);
      /* result */
      MUL(d_tmp,d_result,d_result_tmp);
      COPY(d_result_tmp,d_result);
    }

  //retrieve result
  cudaMemcpy(result, d_result, mem_size, cudaMemcpyDeviceToHost) ;
  checkError(cublasDestroy(handle), "cublasDestroy() error!\n");
  //free memory
  cudaFree(d_Hs);
  cudaFree(d_He);
  cudaFree(d_result);
  cudaFree(d_k1);
  cudaFree(d_k2);
  cudaFree(d_k3);
  cudaFree(d_k4);
  cudaFree(d_k5);
  cudaFree(d_k6);
}

void runTest()
{
  unsigned int uiWA, uiHA, uiWB, uiHB, uiWC, uiHC;
    uiWA = 2;
    uiHA = 2;
    uiWB = 2;
    uiHB = 2;
    uiWC = 2;
    uiHC = 2;
  printf("\nUsing Matrix Sizes: A(%u x %u), B(%u x %u), C(%u x %u)\n\n",
         uiWA, uiHA, uiWB, uiHB, uiWC, uiHC);

  // allocate host memory for matrices A and B
  unsigned int size_A = uiWA * uiHA;
  unsigned int mem_size_A = sizeof(cuDoubleComplex) * size_A;
  cuDoubleComplex* h_A = (cuDoubleComplex*)malloc(mem_size_A);
  unsigned int size_B = uiWB * uiHB;
  unsigned int mem_size_B = sizeof(cuDoubleComplex) * size_B;
  cuDoubleComplex* h_B = (cuDoubleComplex*)malloc(mem_size_B);

  // initialize host memory
  randomInit(h_A, size_A);
  randomInit(h_B, size_B);
  printMatrix(h_A, size_A);
  printMatrix(h_B, size_B);

  // allocate device memory
  cuDoubleComplex* d_A, *d_B, *d_C;
  unsigned int size_C = uiWC * uiHC;
  unsigned int mem_size_C = sizeof(cuDoubleComplex) * size_C;

  // allocate host memory for the result
  cuDoubleComplex* h_C      = (cuDoubleComplex*) malloc(mem_size_C);
  cuDoubleComplex* h_CUBLAS = (cuDoubleComplex*) malloc(mem_size_C);

  cudaMalloc((void**) &d_A, mem_size_A);
  cudaMalloc((void**) &d_B, mem_size_B);

  // copy host memory to device
  cudaMemcpy(d_A, h_A, mem_size_A, cudaMemcpyHostToDevice) ;
  cudaMemcpy(d_B, h_B, mem_size_B, cudaMemcpyHostToDevice) ;

  cudaMalloc((void**) &d_C, mem_size_C);


  //cublas_v2
  cublasHandle_t handle;
  checkError(cublasCreate(&handle), "cublasCreate() error!\n");
  const cuDoubleComplex alpha = make_cuDoubleComplex(1.0,0.0);
  const cuDoubleComplex beta = make_cuDoubleComplex(0.0,0.0);

  //note cublas is column primary!
  //need to transpose the order
  cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, uiWB, uiHA, uiWA, &alpha, d_B, uiWB, d_A, uiWA, &beta, d_C, uiWA);

  // copy result from device to host
  cudaMemcpy(h_CUBLAS, d_C, mem_size_C, cudaMemcpyDeviceToHost) ;
  checkError(cublasDestroy(handle), "cublasDestroy() error!\n");

  // compute reference solution
  cuDoubleComplex* reference = (cuDoubleComplex*)malloc(mem_size_C);
  computeGold(reference, h_A, h_B, uiHA, uiWA, uiWB);

  // check result (CUBLAS)
  printf("Comparing CUBLAS & Host results\n");
  printDiff(reference, h_CUBLAS, uiWC, uiHC, 100, 1.0e-5f);
  printMatrix(h_CUBLAS, size_C);

  // clean up memory
  free(h_A);
  free(h_B);
  free(h_C);
  free(reference);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
}
