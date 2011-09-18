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
    def __init__(self,efield,step,matrix_static,matrix_electric,MATRIX_SIZE,order):
        self.EF = efield
        self.step = step
        self.matrix_static = np.array(matrix_static)
        self.matrix_electric = np.array(matrix_electric)
        if order == 4:
            self.fine_step = 2*step
            self.t_arr = np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,self.fine_step+1)
            self.dt_fine = 2.0*self.EF.time_no_field/(self.fine_step)
            self.dt = 2.0*self.EF.time_no_field/(self.step)
            self.E_arr = []
            for period in range(self.EF.period):
                self.E_arr.append([self.EF.comb_field(self.t_arr[i],period) for i in xrange(self.fine_step+1)])
            self.solve = self.solve_4
        elif order == 5:
            self.dt = 2.0*self.EF.time_no_field/(self.step)
            self.fine_step = 6*step
            self.t_arr = []
            self.E_arr = []            
            for t in np.linspace(-1.0*self.EF.time_no_field,self.EF.time_no_field,self.step+1):
                #0 1/4 3/8 12/13 1 1/2
                self.t_arr.append(t)
                self.t_arr.append(t+1.0/4.0*self.dt)
                self.t_arr.append(t+3.0/8.0*self.dt)
                self.t_arr.append(t+12.0/13.0*self.dt)
                self.t_arr.append(t+1.0*self.dt)
                self.t_arr.append(t+1.0/2.0*self.dt)
            for period in range(self.EF.period):
                self.E_arr.append([self.EF.comb_field(self.t_arr[i],period) for i in xrange(self.fine_step)])
            self.solve = self.solve_5
            
        self.N = MATRIX_SIZE
        self.gpu = GPU_Matrix(self.N)

    def solve_4(self, period,N):
        resultr = np.identity(N,dtype = np.float64)
        resulti = np.zeros([N,N],dtype = np.float64)
        I = np.identity(N,dtype = np.float64)
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
        return resultr + resulti*1j

    def solve_5_gpu(self, period,N):
        resultr = np.identity(N,dtype = np.float64)
        resulti = np.zeros([N,N],dtype = np.float64)
        I = np.identity(N,dtype = np.float64)
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
        k5r = np.zeros([self.N,self.N],dtype = np.float64)
        k5i = np.zeros([self.N,self.N],dtype = np.float64)
        k6r = np.zeros([self.N,self.N],dtype = np.float64)
        k6i = np.zeros([self.N,self.N],dtype = np.float64)        
        for i in xrange(0,self.fine_step,6):
            print i
            t = time.time()
            k1r = (Her*self.E_arr[period][i] + Hsr)
            k1i = (Hei*self.E_arr[period][i] + Hsi)
            tmpr = Her*self.E_arr[period][i+1]+Hsr
            tmpi = Hei*self.E_arr[period][i+1]+Hsi
            self.gpu.matrix_mul(tmpr,tmpi,I+k1r*0.25*self.dt,k1i*0.25*self.dt,k2r,k2i)
            tmpr = Her*self.E_arr[period][i+2]+Hsr
            tmpi = Hei*self.E_arr[period][i+2]+Hsi                        
            self.gpu.matrix_mul(tmpr,tmpi,I+(k1r*3.0/32.0+k2r*9.0/32.0)*self.dt,(k1i*3.0/32.0+k2i*9.0/32.0)*self.dt,k3r,k3i)
            tmpr = Her*self.E_arr[period][i+3]+Hsr
            tmpi = Hei*self.E_arr[period][i+3]+Hsi                                    
            self.gpu.matrix_mul(tmpr,tmpi,I+(k1r*1932.0-k2r*7200.0+k3r*7296.0)*self.dt/2197.0,(k1i*1932.0-k2i*7200.0+k3i*7296.0)*self.dt/2197.0,k4r,k4i)
            tmpr = Her*self.E_arr[period][i+4]+Hsr
            tmpi = Hei*self.E_arr[period][i+4]+Hsi                                    
            self.gpu.matrix_mul(tmpr,tmpi,I+(k1r*439.0/216.0-k2r*8.0+k3r*3680.0/513.0-k4r*845.0/4104.0)*self.dt,(k1i*439.0/216.0-k2i*8.0+k3i*3680.0/513.0-k4i*845.0/4104.0)*self.dt,k5r,k5i)            
            tmpr = Her*self.E_arr[period][i+5]+Hsr
            tmpi = Hei*self.E_arr[period][i+5]+Hsi                                    
            self.gpu.matrix_mul(tmpr,tmpi,I+(-1.0*k1r*8.0/27.0-k2r*2.0-k3r*3544.0/2565.0+k4r*1859.0/4104.0-k5r*11.0/40.0)*self.dt,(-1.0*k1i*8.0/27.0-k2i*2.0-k3i*3544.0/2565.0+k4i*1859.0/4104.0-k5i*11.0/40.0)*self.dt,k6r,k6i)            
            self.gpu.matrix_mul(I+(k1r*16.0/135.0+k3r*6656.0/12825.0+k4r*28561.0/56430.0-k5r*9.0/50.0+k6r*2.0/55.0)*self.dt,(k1i*16.0/135.0+k3i*6656.0/12825.0+k4i*28561.0/56430.0-k5i*9.0/50.0+k6i*2.0/55.0)*self.dt,resultr,resulti,resultr,resulti)
            print time.time() - t
        return resultr + resulti*1j

    def solve_5(self, period,N):
        result = np.identity(N,dtype = np.complex)
        I = np.identity(N,dtype = np.complex)
        He = self.matrix_electric
        Hs = self.matrix_static
        #python version
        for i in xrange(0,self.fine_step,6):
            print i
            t = time.time()
            k1 = (He*self.E_arr[period][i] + Hs)
            tmp = He*self.E_arr[period][i+1]+Hs
            k2 = np.dot(tmp,I+k1*0.25*self.dt)
            tmp = He*self.E_arr[period][i+2]+Hs
            k3 = np.dot(tmp,I+(k1*3.0/32.0+k2*9.0/32.0)*self.dt)
            tmp = He*self.E_arr[period][i+3]+Hs
            k4 = np.dot(tmp,I+(k1*1932.0-k2*7200.0+k3*7296.0)*self.dt/2197.0)
            tmp = He*self.E_arr[period][i+4]+Hs
            k5 = np.dot(tmp,I+(k1*439.0/216.0-k2*8.0+k3*3680.0/513.0-k4*845.0/4104.0)*self.dt)
            tmp = He*self.E_arr[period][i+5]+Hs
            k6 = np.dot(tmp,I+(-1.0*k1*8.0/27.0-k2*2.0-k3*3544.0/2565.0+k4*1859.0/4104.0-k5*11.0/40.0)*self.dt)
            result = np.dot(I+(k1*16.0/135.0+k3*6656.0/12825.0+k4*28561.0/56430.0-k5*9.0/50.0+k6*2.0/55.0)*self.dt,result)
            print time.time() - t
        return result


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

#       test = np.ones([2,2],dtype = complex)*2.0
#       step = self.fine_step
#       dt = float(self.dt)
#       support_code = """
#               static inline void matrix_matrix_product(const  blitz::Array<std::complex<double>,2> & A, const blitz::Array<std::complex<double>,2> & B, blitz::Array<std::complex<double>,2>  & X, int N)
#               {
#                 std::complex<double>  somme;
#                 for (int i=0;i < N;i++){
#                     for (int j=0;j < N;j++){
#                       somme=0.0;
#                               for (int k=0;k < N;k++){
#                                 somme += A(i,k)*B(k,j);
#                               }
#                               X(i,j)=somme;
#                          }
#                        }
#                  }
#                      """
#       code_old = """
#       typedef blitz::Array<std::complex<double>,2> complex_array;
#       complex_array k1(N,N);
#       complex_array k2(N,N);
#       complex_array k3(N,N);
#       complex_array k4(N,N);
#       complex_array temp(N,N);
#       complex_array temp_current(N,N);
#       complex_array temp_result(N,N);
#       for (int i = 0; i < step;i += 2)
#       {
#               printf("   %d \\n",i);
#               printf("   %s \\n", "k1");
#               k1 = (matrix_electric*E_arr(period,i)+matrix_static)*dt;
#               temp = matrix_electric*E_arr(period,i+1)+matrix_static;
#               temp_current = I+k1*0.5;
#               printf("   %s \\n", "k2");
#               matrix_matrix_product(temp,temp_current,k2,N);
#               k2 = k2*dt;
#               temp_current = I + k2*0.5;
#               printf("   %s \\n", "k3");
#               matrix_matrix_product(temp,temp_current,k3,N);
#               k3 = k3*dt;
#               temp = matrix_electric*E_arr(period,i+2)+matrix_static;
#               temp_current = I + k3;
#               printf("   %s \\n", "k4");
#               matrix_matrix_product(temp,temp_current,k3,N);
#               k3 = k3*dt;
#               temp_current = I+1.0/6.0*(k1+2.0*k2+2.0*k3+k4);
#               printf("   %s \\n", "result");
#               matrix_matrix_product(temp_current,result,temp_result,N);
#               printf("   %s \\n", "copy result");
#               result = temp_result;
#       }


#       """

#       code = """
#       typedef blitz::Array<std::complex<double>,2> complex_array;
#       complex_array k1(N,N);
#       complex_array k2(N,N);
#       complex_array k3(N,N);
#       complex_array k4(N,N);
#       complex_array temp(N,N);
#       complex_array temp_current(N,N);
#       complex_array temp_result(N,N);
#       for (int i = 0; i < step;i += 2)
#       {
#               printf("   %d \\n",i);
#               k1 = (matrix_electric*E_arr(period,i)+matrix_static)*dt;
#               temp = matrix_electric*E_arr(period,i+1)+matrix_static;
#               temp_current = I+k1*0.5;
#               k2 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
#               temp_current = I + k2*0.5;
#               k3 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
#               temp = matrix_electric*E_arr(period,i+2)+matrix_static;
#               temp_current = I + k3;
#               k4 = sum(temp(blitz::tensor::i,blitz::tensor::k)*temp_current(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k)*dt;
#               temp_current = I+1.0/6.0*(k1+2.0*k2+2.0*k3+k4);
#               temp_result = sum(temp_current(blitz::tensor::i,blitz::tensor::k)*result(blitz::tensor::k,blitz::tensor::j),blitz::tensor::k);
#               result = temp_result;
# }


#       """
#       # scipy.weave.inline(code,
#       #                  ['step','N','period','result','I','matrix_electric','matrix_static','E_arr','test','dt'],
#       # #                support_code = support_code,
#       #                  type_converters=converters.blitz,
#       #                  compiler = 'gcc')
#       # #                extra_compile_args = ["-O3"])
