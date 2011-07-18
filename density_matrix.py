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
        produce hamiltanian matrix
        now only accept real amplitute
        """
        if i==j:
            return Expolist([Expo(self.omega[i],0)])
        else:
            terms = []
            for k in range(2):# two laser freq
                for l in range(2):
                    mag = -1*self.dipole[i][j]*self.e_amp[k]/2.0
                    freq = (-1)**l*self.nu[k]
                    terms.append(Expo(mag,freq))
            return Expolist(terms)

    def index(self,i,j):
        if i==j:
            return (self.n + (self.n - i + 1)) * i - i
        elif j<i:
            #conjugate
            return (self.n + (self.n - j + 1)) * j - j + (i - j) * 2
        else:
            return -1-i*i+2*j+2*i*(-1+self.n)

    def normalize(self,matrix):
        """
        No need to calculate bloch equation of rho_nn.
        rho_11+rho_22....rho_nn = 1
        """
        for i in range(self.N):
            matrix[self.N-1][i] = 0
        for i in range(self.n):
            matrix[self.N - 1][self.index(i,i)] = 1
        return matrix

    def same_group(self,i,j):
        """
        Are levels in the same group?
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

    def interaction(self,i,j):
        """
        Return the interaction frequency of level i and j
        """
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

    def decoherence(self,system):
        """
        Use this function before to_ri
        """
        for i in range(self.n):
            for j in range(i,self.n):
                for item in self.decoherence_matrix[i][j]:
                    tmp=Expolist([Expo(item[2],0)])
                    t = int(self.index(item[0],item[1]))
                    system[int(self.index(i,j))][t]+=tmp
        return system

    def add_freq(self,system):
        """
        Add frequency dependent terms to system matrix.
        """
        for i in range(self.n):
            for j in range(i+1,self.n):
                system[self.index(i,j)][self.index(i,j)+1] -= self.interaction(i,j)
                system[self.index(i,j) + 1][self.index(i,j)] += self.interaction(i,j)
        return system

    def sweep(self,sweep_n,start,end,points,filename='./test.txt'):
        """
        nu[sweep_n] is sweeped.
        Sweep the frequency and output the result to filename.
        """
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
            self.nu[sweep_n]=self.nu2[sweep_n]+freq
            system_sweep = self.add_freq(system_sweep)
            system_sweep=np.matrix(system_sweep)
            #system_sweep = self.normalize(system_sweep)
            solution = np.linalg.solve(system_sweep,a)
            # print all diagonal element to file
            tmp_str = '%.0f'%freq
            for i in range(self.n):
                tmp_str += ' %.8f'%solution[self.index(i,i),0]
            tmp_str += '\n'
            f.write(tmp_str)

    def von_neumann(self,system):
        for i in range(self.n):
            for j in range(i,self.n):
                #for every upper diagonal terms
                if i==j:
                    for k in range(self.n):
                        tmp = self.hamilton[i][k]*Expolist([Expo(1,-1*self.interaction(k,i))])
                        system[self.index(i,i)][self.index(k,i)] += tmp*(-1*1j)
                        tmp = self.hamilton[k][i]*Expolist([Expo(1,-1*self.interaction(i,k))])
                        system[self.index(i,i)][self.index(i,k)] -= tmp*(-1*1j)
                else:
                    for k in range(self.n):                    
                        tmp = self.hamilton[i][k]*Expolist([Expo(1,-1*self.interaction(k,j))])
                        tmp *= Expolist([Expo(1,self.interaction(i,j))])
                        system[self.index(i,j)][self.index(k,j)] += tmp*(-1*1j)
                        tmp = self.hamilton[k][j]*Expolist([Expo(1,-1*self.interaction(i,k))])
                        tmp *= Expolist([Expo(1,self.interaction(i,j))])                        
                        system[self.index(i,j)][self.index(i,k)] -= tmp*(-1*1j)
        return system

    def to_ri(self,system):
        """
        Convert conjugate terms in system to real and imag terms. 
        """
        # rotating wave approx
        for i in system:
            for j in i:
                j.RWA()
        #change system to number
        for i in range(self.N):
            for j in range(self.N):
                system[i][j] = system[i][j].mag()
        #transform conjugate to real and imag
        for i in range(self.n):
            for j in range(i,self.n):
                for k in range(self.n):
                    for l in range(k,self.n):
                        if (i==j and k!=l):
                            id1 = self.index(i,j)
                            id2 = self.index(k,l)
                            tmp1 = (system[id1][id2]+system[id1][id2+1]).real
                            tmp2 = (system[id1][id2]-system[id1][id2+1]).imag          
                            system[id1][id2],system[id1][id2+1] = tmp1,-1*tmp2
                        elif i!=j:
                            id1 = self.index(i,j)
                            id2 = self.index(k,l)
                            if k==l:
                                tmp = system[id1][id2]
                                system[id1][id2],system[id1+1][id2]=tmp.real,tmp.imag
                            else:
                                id1 = self.index(i,j)
                                id2 = self.index(k,l)
                                tmp1 = system[id1][id2]+system[id1][id2+1]
                                tmp2 = system[id1][id2]-system[id1][id2+1]
                                system[id1][id2],system[id1][id2+1] = tmp1.real,-1*tmp2.imag
                                system[id1+1][id2],system[id1+1][id2+1] = tmp1.imag,tmp2.real
        return system

    def __init__(self,parameter):
        """
        """
        #initialize variables
        print 'initializing...'
        self.n = parameter['n'] #number of state
        self.dipole = parameter['dipole']
        self.omega = parameter['omega']
        self.nu = np.array(parameter['nu'])
        self.e_amp = parameter['e_amp']
        self.nu2 = self.nu.copy()
        self.level_group = parameter['level_group']
        self.decoherence_matrix = parameter['decoherence_matrix']

        self.N = self.n**2 #number of independent density matrix variable
        self.system = [[Expolist() for i in range(self.N)] for j in range(self.N)]
        self.hamilton = np.array([[self.hamiltonian(i,j) for j in range(self.n)]for i in range(self.n)])
        #Start calculate
        print 'von_neumann...'
        self.system = self.von_neumann(self.system) #(H rho-rho H)/i
        self.system = self.decoherence(self.system) #Add decoherence terms
        self.system = self.to_ri(self.system) 
        self.system = self.normalize(self.system) #change last row of system matrix to normalize condition.
        self.system = np.array(self.system) #convert system to np.array, so that it can be solved using scipy.

if __name__ ==  '__main__':
    pass
