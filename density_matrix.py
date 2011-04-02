#!/usr/bin/python2
from __future__ import division
import numpy as np
import sys
from progress_bar import ProgressBar


class System:
    """
    """
    def density_length(self,n):
        """
        Calculate number of independant variable in the density matrix.
        """
        return n + (n * (n - 1))
        """
        n+2(n-1+1)(n-1)/2 because non diagonal terms are complex
        """
        
    def hamiltonian(self,i,j):
        if i==j:
            return self.omega[i]
        else:
            return self.dipole[i][j]
        #rabi frequency
        
    def density_index(self,i,j):
        if j<i:
            raise IOError('only upper diagonal part')
        if i==j:
            return (self.n + (self.n - i + 1)) * i - i
        else:
            return (self.n + (self.n - i + 1)) * i - i + (j - i) * 2 - 1

    def normalize(self,matrix):
        """
        No need to calculate bloch equation of rho_nn.
        rho_11+rho_22....rho_nn = 1
        """
        for i in range(self.n):
            matrix[self.N - 1][self.density_index(i,i)] = 1
        return matrix

    def same_group(self,i,j):
        """
        If levels are in same group?
        """
        for g in self.level_group:
            if (i in g) and (j in g):
                return True
        else:
            return False

    def interaction_freq(self,i,j):
        if self.same_group(i,j):
            return 0
        else:
            if i == 0:
                return -1 * self.nu[j-1]
            elif j == 0:
                return self.nu[i-1]
            else:
                return self.nu[i-1]-self.nu[j-1]
            """
            if (i in self.up) and (j in self.low1):
                return -1 * self.nu[0]
            elif ((i in self.up) and (j in self.low2)):
                return -1 * self.nu[1]
            elif ((i in self.low1) and (j in self.low2)):
                return self.nu[0] - self.nu[1]
            elif (j in self.up) and (i in self.low1):
                return self.nu[0]
            elif ((j in self.up) and (i in self.low2)):
                return self.nu[1]
            elif ((j in self.low1) and (i in self.low2)):
                return self.nu[1] - self.nu[0]
            
            else:
                print 'interaction_freq error'
                print i,j
            """
            
    def rotating_wave_approx(self,freq):
        for i in range(self.nu.__len__()):
            if ((self.nu[i] + freq  == 0) or (self.nu[i] - freq == 0)):
                return [True,i]
        else:
            print 'rotating wave approximation,'
            print self.nu + freq, self.nu - freq
            return [False,0]
        
    def diagonal_part(self,system,i,j):
        complex_tmp = np.zeros(self.N,dtype='complex')
        """
        plus part
        """
        for k in range(self.n):
            if k==j: # diagonal hamiltonian
                complex_tmp[self.density_index(i,k)]+=1j*self.hamiltonian(k,j)
            else:
                df = self.rotating_wave_approx(self.interaction_freq(i,k))
                if k>i:
                    if df[0]:
                        complex_tmp[self.density_index(i,k)] += 1j *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                        complex_tmp[self.density_index(i,k)+1] += -1 *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                else:
                    if df[0]:
                        complex_tmp[self.density_index(k,i)] += 1j *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                        complex_tmp[self.density_index(k,i)+1] += self.e_amp[df[1]]*self.hamiltonian(k,j)/2
        """
        minus part
        """
        for k in range(self.n):
            if k==i: # diagonal hamiltonian
                complex_tmp[self.density_index(k,j)]-=1j*self.hamiltonian(i,k)
            else:
                df = self.rotating_wave_approx(self.interaction_freq(k,j))
                if j>k:
                    if df[0]:
                        complex_tmp[self.density_index(k,j)] -= 1j *self.e_amp[df[1]] * self.hamiltonian(i,k)/2
                        complex_tmp[self.density_index(k,j)+1] -= -1 * self.e_amp[df[1]] * self.hamiltonian(i,k)/2
                else:
                    if df[0]:
                        complex_tmp[self.density_index(j,k)] -= 1j *self.e_amp[df[1]] * self.hamiltonian(i,k)/2
                        complex_tmp[self.density_index(j,k)+1] -= self.e_amp[df[1]] *self.hamiltonian(i,k)/2
        for num in complex_tmp:
            if num.imag != 0:
                raise IOError('complex number in density matrix diagonal!')
        for s in range(self.N):
            system[self.density_index(i,j)][s] = complex_tmp[s].real
        return system

    def non_diagonal_part(self,system,i,j):
        complex_tmp = np.zeros(self.N,dtype='complex')
        """
        plus part
        """
        for k in range(self.n):
            if k==j: # diagonal hamiltonian
                complex_tmp[self.density_index(i,k)] += 1j*self.hamiltonian(k,j)
                complex_tmp[self.density_index(i,k)+1] += -1*self.hamiltonian(k,j)
            else:
                df=self.rotating_wave_approx(self.interaction_freq(i,k)-1*self.interaction_freq(i,j))
                if i==k: #rho_ii
                    if df[0]:
                        complex_tmp[self.density_index(i,k)] += 1j * self.e_amp[df[1]]*self.hamiltonian(k,j)/2
                elif k>i: #rho_ik, upper diagonal
                    if df[0]:
                        complex_tmp[self.density_index(i,k)] += 1j *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                        complex_tmp[self.density_index(i,k)+1] += -1 *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                else: # lower diagonal
                    if df[0]:
                        complex_tmp[self.density_index(k,i)] += 1j *self.e_amp[df[1]]* self.hamiltonian(k,j)/2
                        complex_tmp[self.density_index(k,i)+1] +=  self.e_amp[df[1]]*self.hamiltonian(k,j)/2
        """
        minus part
        """
        for k in range(self.n):
            if k==i: # diagonal hamiltonian
                complex_tmp[self.density_index(k,j)] -= 1j*self.hamiltonian(i,k)
                complex_tmp[self.density_index(k,j)+1] -= -1*self.hamiltonian(i,k)
            else:
                df = self.rotating_wave_approx(self.interaction_freq(k,j)-1*self.interaction_freq(i,j))
                if j==k:
                    if df[0]:
                        complex_tmp[self.density_index(k,j)] -= 1j * self.e_amp[df[1]]*self.hamiltonian(i,k)/2
                elif j>k:
                    if df[0]:
                        complex_tmp[self.density_index(k,j)] -=1j *self.e_amp[df[1]]* self.hamiltonian(i,k)/2
                        complex_tmp[self.density_index(k,j)+1] -= -1 * self.e_amp[df[1]]*self.hamiltonian(i,k)/2
                else:
                    if df[0]:
                        complex_tmp[self.density_index(j,k)] -= 1j * self.e_amp[df[1]]*self.hamiltonian(i,k)/2
                        complex_tmp[self.density_index(j,k)+1] -= self.e_amp[df[1]]*self.hamiltonian(i,k)/2
        for s in range(self.N):
            system[self.density_index(i,j)][s] = complex_tmp[s].real
            system[self.density_index(i,j)+1][s] = complex_tmp[s].imag
        return system

    def decoherence(self,system):
        system[self.density_index(0,0)][self.density_index(0,0)] -= self.Gamma1
        system[self.density_index(1,1)][self.density_index(0,0)] += self.Gamma12
        system[self.density_index(1,1)][self.density_index(1,1)] -= self.gamma1
        system[self.density_index(1,1)][self.density_index(2,2)] += self.gamma2
        system[self.density_index(0,1)][self.density_index(0,1)] -= 0.5*self.Gamma1
        system[self.density_index(0,1)+1][self.density_index(0,1)+1] -= 0.5*self.Gamma1
        system[self.density_index(0,2)][self.density_index(0,2)] -= 0.5*self.Gamma1
        system[self.density_index(0,2)+1][self.density_index(0,2)+1] -= 0.5*self.Gamma1
        system[self.density_index(1,2)][self.density_index(1,2)] -= self.gamma2
        system[self.density_index(1,2)+1][self.density_index(1,2)+1] -= self.gamma2
        return system

    def add_freq(self,system):
        for i in range(self.n):
            for j in range(self.n):
                if j>i:
                    system[self.density_index(i,j)][self.density_index(i,j)+1] += self.interaction_freq(i,j)
                    system[self.density_index(i,j) + 1][self.density_index(i,j)] -= self.interaction_freq(i,j)

        return system

    def sweep(self,e_num,start,end,points,filename='./test.txt'):
        counter = 0 # progress bar's counter
        f=open(filename,'w')# w option will overwrite file if file exist
        prog = ProgressBar(counter, points, 50, mode='fixed', char='#')
        
        a = np.zeros(self.N)
        a[self.N-1] = 1 #1 because rho_11+rho_22 ... =1
        a = np.matrix(a)
        a = a.T
        
        for freq in np.linspace(start,end,points):
            counter +=1
            prog.increment_amount()
            print prog, '\r',
            sys.stdout.flush()
            system_sweep = self.system.copy()
            """
            keep self.system independant of frequency,
            only do frequency dependent operation on system_sweep
            """
            self.nu[e_num]=self.nu2[e_num]+freq
            system_sweep = self.add_freq(system_sweep)
            system_sweep=np.matrix(system_sweep)
            solution = np.linalg.solve(system_sweep,a)
            f.write('%.0f %.8f %.8f %.8f\n'%(freq,solution[0,0],solution[5,0],solution[8,0]))


    def von_neumann(self,system):
        for i in range(self.n):
            for j in range(self.n):
                """
                for every density matrix element
                """
                if j>i and ([i,j] != [self.n-1,self.n-1]):
                    self.non_diagonal_part(system,i,j)
                """
                diagonal element of density matrix
                """
                if j==i and ([i,j] != [self.n-1,self.n-1]):
                    self.diagonal_part(system,i,j)
                else:
                    pass
        return system

    def __init__(self,n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2 ):
        """
        """
        print 'initializing...'
        self.n = n #number of state
        self.omega = omega
        self.dipole = dipole
        self.nu = np.array(nu)
        self.e_amp = e_amp
        self.nu2 = self.nu.copy()
        self.level_group = level_group
        self.Gamma1 = Gamma1
        self.Gamma12 = Gamma12
        self.Gamma13 = Gamma13
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        #self.efreq = np.array([-1*nu[0],nu[0],-1*nu[1],nu[1]]) # frequencies of electric field
        self.N = self.density_length(n) #number of independent density matrix variable
        
        self.system = np.zeros([self.N,self.N])
        self.system = self.normalize(self.system)
        print 'von_neumann'
        self.system = self.von_neumann(self.system)
        self.system = self.decoherence(self.system)

if __name__ ==  '__main__':
    pass
