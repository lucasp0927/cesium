#!/usr/bin/python2
from __future__ import division
from Expo import Expo
from Expolist import Expolist
from sweep_thread import Sweep_Thread
from sweep_multi import *
import multiprocessing
import numpy as np
import sys
from progress_bar import ProgressBar
import time #record calculate time
from constant import *

NOTHING = object()
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
            for k in range(len(self.e_amp)):
                for l in range(2):
                    mag = (-1*self.dipole[k][i][j]*self.e_amp[k][0]/2.0)/HBAR()
                    freq = (-1)**l*self.nu[k]
                    terms.append(Expo(mag,freq))#*(1.0/HBAR()))
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
    
    def interaction(self,i,j,nu = NOTHING):
        """
        Return the interaction frequency of level i and j
        """
        if nu is NOTHING:
            nu = self.nu
        
        if self.same_group(i,j):
            return 0
        else:
            if self.group_number(i) == 0:
                return nu[self.group_number(j)-1]
            elif self.group_number(j) == 0:
                return -1*nu[self.group_number(i)-1]
            else:
                return nu[self.group_number(j)-1]-nu[self.group_number(i)-1]
            
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

    def add_freq(self,system,nu = NOTHING):
        """
        Add frequency dependent terms to system matrix.
        """
        if nu is NOTHING: #because can't use self.xxx as default 
            nu = self.nu
        for i in range(self.n):
            for j in range(i+1,self.n):
                system[self.index(i,j)][self.index(i,j)+1] -= self.interaction(i,j,nu)
                system[self.index(i,j) + 1][self.index(i,j)] += self.interaction(i,j,nu)
        return system

    def sweep_threading(self,sweep_n,start,end,points,filename='./test.txt'):
        """
        nu[sweep_n] is sweeped.
        Sweep the frequency and output the result to filename.
        """
        ###############################
        ##multithread preparation
        ##############################
        threads = 8
        points = points//threads*threads # points per thread
        self.result = [[0.0 for i in range(self.n+1)]for j in range(points)]#this is the matrix which store the result, it will be saved to file later.
        job = self.allocate_job(start,end,points,threads)

        
        ################################
        ##This are codes for progress bar
        ###############################
        prog = ProgressBar(0, points, 50, mode='fixed', char='#')
        ##the linear algebra start here
        a = np.zeros(self.N)
        a[self.N-1] = 1 #1 because rho_11+rho_22 ... =1
        a = np.matrix(a)
        a = a.T

        thread_list = []
        for x in range(threads):
            thread_list.append(Sweep_Thread(self.result,job[x],prog,self.system,self.nu2,a,self.add_freq,self.index,sweep_n,self.n))

        tStart = time.time()            
        for t in thread_list:
            t.start()

        for t in thread_list:
            t.join()
        tStop = time.time()
        print"spend",(tStop - tStart),"second"
            
        self.sweep_save_file(filename,points)

    def sweep_multiprocessing(self,sweep_n,start,end,points,filename='./test.txt'):
        """
        nu[sweep_n] is sweeped.
        Sweep the frequency and output the result to filename.
        """
        ###############################
        ##multiprocessing preparation
        ##############################
        core = 10
        points = points//core*core # points per thread
        self.result = [[0.0 for i in range(self.n+1)]for j in range(points)]#this is the matrix which store the result, it will be saved to file later.
        job = self.allocate_job(start,end,points,core)

        
        ################################
        ##This are codes for progress bar
        ###############################
        prog = ProgressBar(0, points, 50, mode='fixed', char='#')
        ##the linear algebra start here
        a = np.zeros(self.N)
        a[self.N-1] = 1 #1 because rho_11+rho_22 ... =1
        a = np.matrix(a)
        a = a.T

        done_queue = multiprocessing.Queue()
        process_list = []
        for x in range(core):
            process_list.append(multiprocessing.Process(target = sweep_mp,args = (job[x],self.system,self.nu2,a,self.add_freq,self.index,sweep_n,self.n,done_queue)))

        tStart = time.time()
        print 'start'
        for p in process_list:
            p.start()

        stop_num = 0
        while stop_num != core:
            a = done_queue.get()
            if a == 'STOP':
                stop_num += 1
            else:
                self.result[a[0]] = a[1]
                prog.increment_amount()
                print prog, '\r',
                sys.stdout.flush()

        print '\n'
        for p in process_list:
            p.join()
            print "%s.exitcode = %s" %(p.name, p.exitcode)

        tStop = time.time()
        print"spend",(tStop - tStart),"second"
            
        self.sweep_save_file(filename,points)

    def allocate_job(self,start,end,points,threads):
        points = points//threads*threads # points per thread
        ppt = points//threads #points per threads
        space = zip(range(points),np.linspace(float(start),float(end),points))
        job = []
        for i in range(threads):
            job.append(space[i*ppt:(i+1)*ppt])
        return job
    
    def sweep_save_file(self,filename,points):
        print "save file"
        f=open(filename,'w')# w option will overwrite file if file exist            
        tmp_str = "#average power: %fw/m^2\n" %(self.power)
        f.write(tmp_str)
        for i in range(points):
            tmp_str = '%.10f '%(self.result[i][0]/(2.0*np.pi)) # convert to hz
            if self.d1 == 1:
                tmp_str += '%.10f %.10f %.10f' %(np.sum(self.result[i][1:17]),np.sum(self.result[i][17:26]),np.sum(self.result[i][26:33]))
            elif self.d1 == 2:
                tmp_str += '%.10f %.10f %.10f' %(np.sum(self.result[i][1:33]),np.sum(self.result[i][33:42]),np.sum(self.result[i][42:49]))                
            # for j in range(self.n):
            #     tmp_str += ' %.10f'%self.result[i][j+1]
            tmp_str += '\n'
            f.write(tmp_str)
        f.close()
        
    def von_neumann(self,system):
        for i in range(self.n):
            for j in range(i,self.n):
                for k in range(self.n):                    
                        tmp = self.hamilton[i][k]*Expolist([Expo(1,-1*self.interaction(k,j))])
                        tmp *= Expolist([Expo(1,self.interaction(i,j))])#result from differential in left hand side
                        system[self.index(i,j)][self.index(k,j)] += tmp*(-1j)

                        #minus part
                        tmp = self.hamilton[k][j]*Expolist([Expo(1,-1*self.interaction(i,k))])
                        tmp *= Expolist([Expo(1,self.interaction(i,j))])                        
                        system[self.index(i,j)][self.index(i,k)] -= tmp*(-1j)                
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
#        print self.system
        #transform conjugate to real and imag
        for i in range(self.n):
            for j in range(i,self.n):
                id1 = self.index(i,j)                
                if i == j:
                    for k in range(self.n):
                        for l in range(k,self.n):
                            if k != l:
                                id2 = self.index(k,l)
                                tmp1 = (system[id1][id2]+system[id1][id2+1]).real
                                tmp2 = (system[id1][id2]-system[id1][id2+1]).imag          
                                system[id1][id2],system[id1][id2+1] = tmp1,-1*tmp2
                else:
                    for k in range(self.n):
                        for l in range(k,self.n):
                            id2 = self.index(k,l)
                            if k==l:
                                tmp = system[id1][id2]
                                system[id1][id2],system[id1+1][id2]=tmp.real,tmp.imag
                            else:
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
        self.d1 = parameter['d1']
        self.power = parameter['power']
        self.N = self.n**2 #number of independent density matrix variable
        self.result = []
        print '  initializing system matrix...'        
        self.system = [[Expolist() for i in xrange(self.N)] for j in xrange(self.N)]
        print '  initializing hamiltonian...'
        self.hamilton = np.array([[self.hamiltonian(i,j) for j in range(self.n)]for i in range(self.n)])
        #Start calculate
        print 'von_neumann...'
        self.system = self.von_neumann(self.system) #(H rho-rho H)/i
        print 'decoherence...'
        self.system = self.decoherence(self.system) #Add decoherence terms
        print 'to_ri'        
        self.system = self.to_ri(self.system)
#        print self.system
        self.system = self.normalize(self.system) #change last row of system matrix to normalize condition.
        self.system = np.array(self.system) #convert system to np.array, so that it can be solved using scipy.
        #choose threading or multiprocessing here, function pointer
        self.sweep = self.sweep_multiprocessing
        
if __name__ ==  '__main__':
    pass
