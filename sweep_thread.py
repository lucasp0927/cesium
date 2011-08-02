#!/usr/bin/python2
import threading
import numpy as np
import sys
from copy import copy
class Sweep_Thread(threading.Thread):
    def __init__(self,result,job,prog,system,nu2,a,add_freq,index,sweep_n,n):
        threading.Thread.__init__(self)        
        self.result = result
        self.job = job
        self.prog = prog
        self.system = system
        self.nu2 = nu2
        self.add_freq = add_freq
        self.index = index        
        self.a = a
        self.sweep_n = sweep_n
        self.n = n
        self.sub_result = [[0.0 for i in range(self.n+1)]for j in range(len(job))]
        
    def run(self):
        nu = self.nu2.copy()

        for freq in self.job:
            ## progress bar
            sys.stdout.flush()
            self.prog.increment_amount()
            print self.prog, '\r',
            sys.stdout.flush()

            ## copy system
            system_sweep = self.system.copy()
            """
            keep self.system independant of frequency,
            only do frequency dependent operation on system_sweep
            """
            nu[self.sweep_n]=self.nu2[self.sweep_n]+freq[1]
            system_sweep = self.add_freq(system_sweep,nu)#add freq dependent terms
            system_sweep=np.matrix(system_sweep)
            solution = np.linalg.solve(system_sweep,self.a)#solve linear algebra system
            self.sub_result[freq[0]-self.job[0][0]][0] = freq[1]
            # save all diagonal element to matrix
            for i in range(self.n):
                #self.result[freq[0]][i+1] = solution[self.index(i,i),0]
                self.sub_result[freq[0]-self.job[0][0]][i+1]=solution[self.index(i,i),0]
        for i in self.job:
            self.result[i[0]] = self.sub_result[i[0]-self.job[0][0]]
        
if __name__ == '__main__':
    pass
