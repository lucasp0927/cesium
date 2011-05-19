#!/usr/bin/python2
from __future__ import division
from Expo import Expo
from Expolist import Expolist
import numpy as np
import sys
from progress_bar import ProgressBar

class System:
    """
    """
    def hamiltonian(self,i,j):
        """
        now only for real dipole moment
        """
        if i==j:
            return Expolist([Expo(self.omega[i],0)])
        else:
            terms = []
            for k in range(2):
                for l in range(2):
                    mag = self.dipole[i][j]*self.e_amp[k]
                    freq = (-1)**l*self.nu[k]
                    terms.append(Expo(mag,freq))
            return Expolist(terms)
        #rabi frequency

    def density_index(self,i,j):
        if i==j:
            return (self.n + (self.n - i + 1)) * i - i
        elif j<i:
            #conjugate
            return (self.n + (self.n - j + 1)) * j - j + (i - j) * 2
        else:
            return (self.n + (self.n - i + 1)) * i - i + (j - i) * 2 - 1

    def normalize(self,matrix):
        """
        No need to calculate bloch equation of rho_nn.
        rho_11+rho_22....rho_nn = 1
        """
        for i in range(self.N):
            matrix[self.N-1][i] = 0
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

    def group_number(self,i):
        for k in range(self.level_group.__len__()):
            if i in self.level_group[k]:
                return k
        else:
            raise IOError('no level found in group')

    def interaction_freq(self,i,j):
        if self.same_group(i,j):
            return 0
        else:
            if self.group_number(i) == 0:
                return self.nu[self.group_number(j)-1]
            elif self.group_number(i) == 1  and self.group_number(j) == 2:
                return self.nu[1] - self.nu[0]
            elif self.group_number(j) == 1  and self.group_number(i) == 2:
                return self.nu[0] - self.nu[1]
            elif self.group_number(j) == 0:
                return -1*self.nu[self.group_number(i)-1]
            else:
                print 'interaction_freq error'
                print i,j

    def rotating_wave_approx(self,freq):
        for i in range(self.nu.__len__()):
            if ((self.nu[i] + freq  == 0) or (self.nu[i] - freq == 0)):
                return [True,i]
        else:
            #print 'rotating wave approximation,'
            #print self.nu + freq, self.nu - freq
            return [False,0]


    def decoherence(self,system):
        """
        dont put terms in self.density_index(n,n), since it is used to normalize.
        """

        system[self.density_index(0,0)][self.density_index(0,0)] -= self.Gamma1
        system[self.density_index(1,1)][self.density_index(0,0)] += self.Gamma12
#        system[self.density_index(2,2)][self.density_index(0,0)] += self.Gamma13
        system[self.density_index(0,1)][self.density_index(0,1)] -= 0.5*self.Gamma1
        system[self.density_index(0,1)+1][self.density_index(0,1)+1] -= 0.5*self.Gamma1
        system[self.density_index(0,2)][self.density_index(0,2)] -= 0.5*self.Gamma1
        system[self.density_index(0,2)+1][self.density_index(0,2)+1] -= 0.5*self.Gamma1
        system[self.density_index(1,1)][self.density_index(1,1)] -= self.gamma1
        system[self.density_index(1,1)][self.density_index(2,2)] += self.gamma1
#        system[self.density_index(2,2)][self.density_index(2,2)] -= self.gamma2
#        system[self.density_index(2,2)][self.density_index(1,1)] += self.gamma2
        system[self.density_index(1,2)][self.density_index(1,2)] -= self.gamma2
        system[self.density_index(1,2)+1][self.density_index(1,2)+1] -= self.gamma2

        '''
        system[self.density_index(1,1)][self.density_index(1,1)] -= self.Gamma1
        system[self.density_index(2,2)][self.density_index(1,1)] += self.Gamma12
        system[self.density_index(0,0)][self.density_index(1,1)] += self.Gamma12
#        system[self.density_index(3,3)][self.density_index(1,1)] += self.Gamma13
        system[self.density_index(1,2)][self.density_index(1,2)] -= 0.5*self.Gamma1
        system[self.density_index(1,2)+2][self.density_index(1,2)+1] -= 0.5*self.Gamma1
        system[self.density_index(1,3)][self.density_index(1,3)] -= 0.5*self.Gamma1
        system[self.density_index(1,3)+2][self.density_index(1,3)+1] -= 0.5*self.Gamma1
        system[self.density_index(2,2)][self.density_index(2,2)] -= self.gamma1
        system[self.density_index(2,2)][self.density_index(3,3)] += self.gamma1
        system[self.density_index(0,0)][self.density_index(2,2)] -= self.gamma1
        system[self.density_index(0,0)][self.density_index(1,1)] += self.gamma1
#        system[self.density_index(3,3)][self.density_index(3,3)] -= self.gamma2
#        system[self.density_index(3,3)][self.density_index(2,2)] += self.gamma3
        system[self.density_index(2,3)][self.density_index(2,3)] -= self.gamma2
        system[self.density_index(2,3)+1][self.density_index(2,3)+1] -= self.gamma2
        '''
        return system

    def add_freq(self,system):
        for i in range(self.n):
            for j in range(self.n):
                if j>i:
                    system[self.density_index(i,j)][self.density_index(i,j)+1] -= self.interaction_freq(i,j)
                    system[self.density_index(i,j) + 1][self.density_index(i,j)] += self.interaction_freq(i,j)

        return system

    def sweep(self,start,end,points,filename='./test.txt'):
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
            self.nu[0]=self.nu2[0]+freq
            system_sweep = self.add_freq(system_sweep)
            system_sweep=np.matrix(system_sweep)
            #system_sweep = self.normalize(system_sweep)
            solution = np.linalg.solve(system_sweep,a)
            # print all diagonal element to file
            tmp_str = '%.0f'%freq
            for i in range(self.n):
                tmp_str += ' %.8f'%solution[self.density_index(i,i),0]
            tmp_str += '\n'
            f.write(tmp_str)

    def von_neumann(self,system):
        exposystem=[[Expolist() for x in range(self.N)]for y in range(self.N)]

        for i in range(self.n):
            for j in range(self.n):
                #conjugate
                pivot = self.density_index(i,j)
                for l in range(self.n):
                    tmp = Expolist([Expo(1,self.interaction_freq(l,j))])
                    exposystem[pivot][self.density_index(l,j)] += self.hamiltonian(i,l)*tmp
                    tmp = Expolist([Expo(1,self.interaction_freq(i,l))])
                    exposystem[pivot][self.density_index(i,l)] -= self.hamiltonian(l,j)*tmp
        for i in range(self.n):
            for j in range(self.n):
                pivot = self.density_index(i,j)
                tmp = Expolist([Expo(1/1j,self.interaction_freq(i,j))])#[H,rho]/i<- is here
                for l in range(self.N):
                    exposystem[pivot][l] *= tmp
        print exposystem[0][1]
        print exposystem[0][2]
        '''
        change to RE IM
        '''

        for i in range(self.n):
            for j in range(i+1,self.n):
                # for all upper diagnal element
                index = self.density_index(i,j)
                for l in exposystem:
                    l[index],l[index+1] = l[index]+l[index+1],(l[index]-l[index+1])*1j

        for i in range(self.n):
            for j in range(i+1,self.n):
                index = self.density_index(i,j)
                for l in range(self.N):
                    exposystem[index][l],exposystem[index+1][l] = (exposystem[index][l]+exposystem[index+1][l])*0.5,(exposystem[index+1][l]-exposystem[index][l])*0.5j

        '''
        RWA
        '''
        for i in exposystem:
            for j in i:
                j.RWA()

        for i in range(self.N):
            for j in range(self.N):
                if exposystem[i][j].mag().imag != 0:
                    print 'non real!'
                system[i][j] = exposystem[i][j].mag().real
        '''
        assert all terms in system is real
        '''

        return system

    def __init__(self,n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2 ):
        """
        """
        print 'initializing...'
        self.n = n #number of state
        self.dipole,self.omega = dipole,omega
        self.nu = np.array(nu)
        self.e_amp = e_amp
        self.nu2 = self.nu.copy()
        self.level_group = level_group
        self.gamma2,self.gamma1,self.Gamma13,self.Gamma12,self.Gamma1 = gamma2,gamma1,Gamma13,Gamma12,Gamma1
        self.N = self.n*self.n #number of independent density matrix variable
        self.system = np.zeros([self.N,self.N])
        print 'von_neumann...'
        self.system = self.von_neumann(self.system)
        self.system = self.decoherence(self.system)
        self.system = self.normalize(self.system)        
        print 'system\n', self.system


if __name__ ==  '__main__':
    n=3
    omega = [105E10,9E9,0]
    dipole=[[0,1000000,1000000],
            [1000000,0,0],
            [1000000,0,0]]
    #remember to /2
    nu = [105E10-9E9,105E10] # on resonence
    e_amp = [1,1]
    level_group = [[0],[1],[2]]
    #decoherence
    Gamma1 = 5000000
    Gamma12 = 2500000
    Gamma13 = 2500000
    gamma1 = 10000
    gamma2 = 10000

    filename = './test.txt'
    system = System(n,omega,dipole,nu,e_amp,level_group,Gamma1,Gamma12,Gamma13,gamma1,gamma2)
    #system.sweep(-1E7,1E7,400,'./test.txt')#TODO: add file name
    #plot(n)
