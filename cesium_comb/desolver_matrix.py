#!/usr/bin/python2
import numpy as np
from electricfield import Electricfield
from scipy.weave import converters
from gpu_matrix import GPU_Matrix
import scipy
import time
from pycuda import driver, compiler, gpuarray, tools
import pycuda.autoinit

class DESolver_matrix(object):
    """
    """
    def __init__(self,efield,step,matrix_static,matrix_electric,MATRIX_SIZE):
        self.EF = efield
        self.step = step
        self.matrix_static = np.array(matrix_static)
        self.matrix_electric = np.array(matrix_electric)
        self.fine_step = 2*step
        self.t_arr = np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,self.fine_step+1)
        self.dt_fine = 2.0*self.EF.time_no_field/(self.fine_step)
        self.dt = 2.0*self.EF.time_no_field/(self.step)        
        self.E_arr = []
        for period in range(self.EF.period):
            self.E_arr.append([self.EF.comb_field(self.t_arr[i],period) for i in xrange(self.fine_step+1)])
	self.N = MATRIX_SIZE
	self.gpu = GPU_Matrix(self.N)
	    
    def solve(self, period,N):
	resultr = np.identity(N,dtype = np.float64)
	resulti = np.zeros([N,N],dtype = np.float64)	
	I = np.identity(N,dtype = np.float64)	
	# E_arr = np.array(self.E_arr)
	# matrix_electric = self.matrix_electric
	# matrix_static = self.matrix_static
	Her = self.matrix_electric.real
	Hei = self.matrix_electric.imag
	Hsr = self.matrix_static.real
	Hsi = self.matrix_static.imag
	#python version
	k2r = np.zeros([self.N,self.N],dtype = np.float64)
	k2i = np.zeros([self.N,self.N],dtype = np.float64)
	k3r = np.zeros([self.N,self.N],dtype = np.float64)
	k3i = np.zeros([self.N,self.N],dtype = np.float64)
	k4r = np.zeros([self.N,self.N],dtype = np.float64)
	k4i = np.zeros([self.N,self.N],dtype = np.float64)			
        for i in xrange(0,self.fine_step,2):
	    print i
	    t = time.time()
	    k1r = (Her*self.E_arr[period][i] + Hsr)*self.dt
	    k1i = (Hei*self.E_arr[period][i] + Hsi)*self.dt
	    tmpr = Her*self.E_arr[period][i+1]+Hsr
	    tmpi = Hei*self.E_arr[period][i+1]+Hsi
            self.gpu.matrix_mul(tmpr,tmpi,I+k1r*0.5,k1i*0.5,k2r,k2i)
	    k2r *= self.dt
	    k2i *= self.dt
            self.gpu.matrix_mul(tmpr,tmpi,I+k2r*0.5,k2i*0.5,k3r,k3i)
	    k3r *= self.dt
	    k3i *= self.dt	    
            self.gpu.matrix_mul((Her*self.E_arr[period][i+2]+Hsr),(Hei*self.E_arr[period][i+2]+Hsi),I+k3r,k3i,k4r,k4i)
	    k4r *= self.dt
	    k4i *= self.dt	    
	    self.gpu.matrix_mul(I+1.0/6.0*(k1r+2.0*k2r+2.0*k3r+k4r),1.0/6.0*(k1i+2.0*k2i+2.0*k3i+k4i),resultr,resulti,resultr,resulti)
	    print time.time() - t
        return resultr + resulti*1j


        # for i in xrange(0,self.fine_step,2):
	#     t = time.time()
	#     print i
        #     k1 =  (self.*self.E_arr[period][i]+self.matrix_static)*self.dt
        #     k2 =  self.gpu.matrix_mul((self.matrix_electric*self.E_arr[period][i+1]+self.matrix_static),I+k1*0.5)*self.dt
        #     k3 =  self.gpu.matrix_mul((self.matrix_electric*self.E_arr[period][i+1]+self.matrix_static),I+k2*0.5)*self.dt
        #     k4 =  self.gpu.matrix_mul((self.matrix_electric*self.E_arr[period][i+2]+self.matrix_static),I+k3)*self.dt
	#     result = self.gpu.matrix_mul(I+1.0/6.0*(k1+2.0*k2+2.0*k3+k4),result)
	#     print time.time() - t
        # return result
		
# 	test = np.ones([2,2],dtype = complex)*2.0
# 	step = self.fine_step
# 	dt = float(self.dt)
# 	support_code = """
# 		static inline void matrix_matrix_product(const  blitz::Array<std::complex<double>,2> & A, const blitz::Array<std::complex<double>,2> & B, blitz::Array<std::complex<double>,2>  & X, int N)
# 		{
#                 std::complex<double>  somme;
#                 for (int i=0;i < N;i++){
#                     for (int j=0;j < N;j++){
#                    	somme=0.0;
#                        	for (int k=0;k < N;k++){
#                        	  somme += A(i,k)*B(k,j);
#                             	}
#                          	X(i,j)=somme;
# 			   }
#                        }
#                  }
# 		       """
# 	code_old = """
# 	typedef blitz::Array<std::complex<double>,2> complex_array;
# 	complex_array k1(N,N);
# 	complex_array k2(N,N);
# 	complex_array k3(N,N);
# 	complex_array k4(N,N);
# 	complex_array temp(N,N);
# 	complex_array temp_current(N,N);
# 	complex_array temp_result(N,N);	    	    			   
# 	for (int i = 0; i < step;i += 2)
# 	{
# 	        printf("   %d \\n",i);
# 		printf("   %s \\n", "k1");
# 		k1 = (matrix_electric*E_arr(period,i)+matrix_static)*dt;
# 		temp = matrix_electric*E_arr(period,i+1)+matrix_static;
# 		temp_current = I+k1*0.5;
# 		printf("   %s \\n", "k2");
# 		matrix_matrix_product(temp,temp_current,k2,N);
# 		k2 = k2*dt;
# 		temp_current = I + k2*0.5;
# 		printf("   %s \\n", "k3");
# 		matrix_matrix_product(temp,temp_current,k3,N);
# 		k3 = k3*dt;
# 		temp = matrix_electric*E_arr(period,i+2)+matrix_static;
# 		temp_current = I + k3;
# 		printf("   %s \\n", "k4");
# 		matrix_matrix_product(temp,temp_current,k3,N);
# 		k3 = k3*dt;				     
# 		temp_current = I+1.0/6.0*(k1+2.0*k2+2.0*k3+k4);
# 		printf("   %s \\n", "result");
# 		matrix_matrix_product(temp_current,result,temp_result,N);				     
# 		printf("   %s \\n", "copy result");		
# 		result = temp_result;
# 	}


# 	"""

# 	code = """
# 	typedef blitz::Array<std::complex<double>,2> complex_array;
# 	complex_array k1(N,N);
# 	complex_array k2(N,N);
# 	complex_array k3(N,N);
# 	complex_array k4(N,N);
# 	complex_array temp(N,N);
# 	complex_array temp_current(N,N);
# 	complex_array temp_result(N,N);	    	    	
# 	for (int i = 0; i < step;i += 2)
# 	{
# 	        printf("   %d \\n",i);
# 		k1 = (matrix_electric*E_arr(period,i)+matrix_static)*dt;
# 		temp = matrix_electric*E_arr(period,i+1)+matrix_static;
# 		temp_current = I+k1*0.5;
# 		k2 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
# 		temp_current = I + k2*0.5;
# 		k3 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
# 		temp = matrix_electric*E_arr(period,i+2)+matrix_static;
# 		temp_current = I + k3;
# 		k4 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
# 		temp_current = I+1.0/6.0*(k1+2.0*k2+2.0*k3+k4);
# 		temp_result = sum(temp_current(blitz::tensor::i,blitz::tensor::k)*result(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k);
# 		result = temp_result;
# }


# 	"""	    	
#     	# scipy.weave.inline(code,
#     	# 		   ['step','N','period','result','I','matrix_electric','matrix_static','E_arr','test','dt'],
# 	# #		   support_code = support_code,
#     	# 		   type_converters=converters.blitz,
# 	# 		   compiler = 'gcc')
# 	# #		   extra_compile_args = ["-O3"])

