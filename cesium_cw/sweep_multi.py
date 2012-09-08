#!/usr/bin/python2
import numpy as np
from copy import copy

def sweep_mp(job,system,nu2,a,add_freq,index,sweep_n,n,output):
    nu = nu2.copy()
    sub_result = [0.0 for i in range(n+1)]
    for freq in job:
        ## copy system
        system_sweep = system.copy()
        """
        keep self.system independant of frequency,
        only do frequency dependent operation on system_sweep
        """
        """
        sweep both frequency
        """
        #nu[sweep_n]=nu2[sweep_n]+freq[1]
        nu[0]=nu2[0]+freq[1]
        nu[1]=nu2[1]-freq[1]                
        system_sweep = add_freq(system_sweep,nu)#add freq dependent terms
        system_sweep = np.matrix(system_sweep)
        solution = np.linalg.solve(system_sweep,a)#solve linear algebra system
        sub_result[0] = freq[1]
        # save all diagonal element to matrix
        for i in range(n):
            sub_result[i+1]=solution[index(i,i),0]
        #output.put(((job[0][0],job[-1][0]),copy(sub_result)))
        output.put((freq[0],copy(sub_result)))
    output.put('STOP')
