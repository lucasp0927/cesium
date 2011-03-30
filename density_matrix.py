#!/usr/bin/python2
from __future__ import division
import numpy as np

n=3 #number of states
N=0

omega = [351E12,9E9,0]

dipole=[[0,500000,500000],
        [500000,0,0],
        [500000,0,0]]
#remember to /2

nu=[351E12-9E9,351E12]
nu2=[351E12-9E9,351E12]

efreq = np.array([-1*nu[0],nu[0],-1*nu[1],nu[1]]) # frequencies of electric field

up = [0]
low1 = [1]
low2 = [2]
groups = [up,low1,low2]

#decoherence
Gamma1 = 5000000
Gamma12 = 2500000
Gamma13 = 2500000
gamma1 = 10000
gamma2 = 10000


def density_length(n):
    """
    Calculate number of independant variable in the density matrix.
    """
    return n+(n*(n-1)) #n+2(n-1+1)(n-1)/2 because non diagonal terms are complex

def hamiltonian(i,j):
    if i==j:
        return omega[i]
    else:
        return dipole[i][j]

def density_index(i,j):
    if j<i:
        raise IOError('only upper diagonal part')
    if i==j:
        return (n+(n-i+1))*i-i
    else:
        return (n+(n-i+1))*i-i+(j-i)*2-1

#def reverse_density_index(k):


def normalize(matrix):
    """
    no need to calculate bloch equation of rho_nn
    """
    for i in range(n):
        matrix[N-1][density_index(i,i)]=1

def same_group(i,j):
    for g in groups:
        if (i in g) and (j in g):
            return True
    else:
        return False

def interaction_freq(i,j):
    if same_group(i,j):
        return 0
    else:
        if (i in up) and (j in low1):
            return -1*nu[0]
        elif ((i in up) and (j in low2)):
            return -1*nu[1]
        elif ((i in low1) and (j in low2)):
            return nu[0]-nu[1]
        elif (j in up) and (i in low1):
            return nu[0]
        elif ((j in up) and (i in low2)):
            return nu[1]
        elif ((j in low1) and (i in low2)):
            return nu[1]-nu[0]
        else:
            print 'interaction_freq error\n'
            print i,j

def diagonal_part(system,i,j):
    complex_tmp = np.zeros(N,dtype='complex')
    """
    plus part
    """
    for k in range(n):
        if k==j: # diagonal hamiltonian
            complex_tmp[density_index(i,k)]+=1j*hamiltonian(k,j)
        else:
            if k>i:
                if 0 in efreq+interaction_freq(i,k):
                    complex_tmp[density_index(i,k)]+=1j*hamiltonian(k,j)
                    complex_tmp[density_index(i,k)+1]+= -1*hamiltonian(k,j)
            else:
                if 0 in efreq+interaction_freq(i,k):
                    complex_tmp[density_index(k,i)]+=1j*hamiltonian(k,j)
                    complex_tmp[density_index(k,i)+1]+= hamiltonian(k,j)
    """
    minus part
    """
    for k in range(n):
        if k==i: # diagonal hamiltonian
            complex_tmp[density_index(k,j)]-=1j*hamiltonian(i,k)
        else:
            if j>k:
                if 0 in efreq+interaction_freq(k,j):
                    complex_tmp[density_index(k,j)]-=1j*hamiltonian(i,k)
                    complex_tmp[density_index(k,j)+1]-= -1*hamiltonian(i,k)
            else:
                if 0 in efreq+interaction_freq(k,j):
                    complex_tmp[density_index(j,k)]-=1j*hamiltonian(i,k)
                    complex_tmp[density_index(j,k)+1]-= hamiltonian(i,k)
    for num in complex_tmp:
        if num.imag !=0:
            raise IOError('complex number in density matrix diagonal!')
    for s in range(N):
        system[density_index(i,j)][s] = complex_tmp[s].real

def non_diagonal_part(system,i,j):
    complex_tmp = np.zeros(N,dtype='complex')
    """
    plus part
    """
    for k in range(n):
        if k==j: # diagonal hamiltonian
            if interaction_freq(i,k)-interaction_freq(i,j) == 0: #always same can be delete
                complex_tmp[density_index(i,k)] += 1j*hamiltonian(k,j)
                complex_tmp[density_index(i,k)+1] += -1*hamiltonian(k,j)
        else:
            if i==k:
                if 0 in efreq+interaction_freq(i,k)-interaction_freq(i,j):
                    complex_tmp[density_index(i,k)]+=1j*hamiltonian(k,j)
            elif k>i:
                if 0 in efreq+interaction_freq(i,k)-interaction_freq(i,j):
                    complex_tmp[density_index(i,k)]+=1j*hamiltonian(k,j)
                    complex_tmp[density_index(i,k)+1]+= -1*hamiltonian(k,j)
            else:
                if 0 in efreq+interaction_freq(i,k)-interaction_freq(i,j):
                    complex_tmp[density_index(k,i)]+=1j*hamiltonian(k,j)
                    complex_tmp[density_index(k,i)+1]+= hamiltonian(k,j)
    """
    minus part
    """
    for k in range(n):
        if k==i: # diagonal hamiltonian
            if interaction_freq(k,j)-interaction_freq(i,j) == 0: #always same can be delete
                complex_tmp[density_index(k,j)] -= 1j*hamiltonian(i,k)
                complex_tmp[density_index(k,j)+1] -= -1*hamiltonian(i,k)
        else:
            if j==k:
                if 0 in efreq+interaction_freq(k,j)-interaction_freq(i,j):
                    complex_tmp[density_index(k,j)]-= 1j*hamiltonian(i,k)
            elif j>k:
                if 0 in efreq+interaction_freq(k,j)-interaction_freq(i,j):
                    complex_tmp[density_index(k,j)]-=1j*hamiltonian(i,k)
                    complex_tmp[density_index(k,j)+1]-= -1*hamiltonian(i,k)
            else:
                if 0 in efreq+interaction_freq(k,j)-interaction_freq(i,j):
                    complex_tmp[density_index(j,k)]-=1j*hamiltonian(i,k)
                    complex_tmp[density_index(j,k)+1]-= hamiltonian(i,k)
    for s in range(N):
        system[density_index(i,j)][s] = complex_tmp[s].real
        system[density_index(i,j)+1][s] = complex_tmp[s].imag

def decoherence(system):
    system[density_index(0,0)][density_index(0,0)]-=Gamma1
    system[density_index(1,1)][density_index(0,0)]+=Gamma12
    system[density_index(1,1)][density_index(1,1)]-=gamma1
    system[density_index(1,1)][density_index(2,2)]+=gamma2
    system[density_index(0,1)][density_index(0,1)]-=0.5*Gamma1
    system[density_index(0,1)+1][density_index(0,1)+1]-=0.5*Gamma1
    system[density_index(0,2)][density_index(0,2)]-=0.5*Gamma1
    system[density_index(0,2)+1][density_index(0,2)+1]-=0.5*Gamma1
    system[density_index(1,2)][density_index(1,2)]-=gamma2
    system[density_index(1,2)+1][density_index(1,2)+1]-=gamma2

def sweep_freq(system):
    for i in range(n):
        for j in range(n):
            if j>i:
                system[density_index(i,j)][density_index(i,j)+1]+=interaction_freq(i,j)
                system[density_index(i,j)+1][density_index(i,j)]-=interaction_freq(i,j)

def von_neumann(system):
    for i in range(n):
        for j in range(n):
            #for every density matrix element
            if j>i and ([i,j] != [n-1,n-1]):
                non_diagonal_part(system,i,j)
            """
            diagonal element of density matrix
            """
            if j==i and ([i,j] != [n-1,n-1]):
                diagonal_part(system,i,j)
            else:
                pass

if __name__ ==  '__main__':
    N=density_length(n)
    a=np.zeros(N)
    a[N-1]=1
    a=np.matrix(a)
    a=a.T
    system=np.zeros([N,N])
    normalize(system)
    von_neumann(system)
    decoherence(system)
    for freq in range(int(-1E7),int(1E7),int(1E4)):
        nu[0]=nu2[0]+freq
        system_sweep=system.copy()
        sweep_freq(system_sweep)
        system_sweep=np.matrix(system_sweep)
        print freq,np.linalg.solve(system_sweep,a)[0,0]
