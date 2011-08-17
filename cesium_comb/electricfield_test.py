#!/usr/bin/python2
from electricfield import Electricfield
import numpy as np

if __name__ == '__main__':
    para = {}
    para['Tr'] = 91.9262177e6
    para['mu_c'] = 351.72571850e12
    para['PSI'] = 3.0/4.0*2.0*np.pi
    para['E_0'] = 1.0
    para['tao'] = 3e-14
    print 'testing period calculation'
    for i in range(1,300):
        for j in range(1,300):
            para['PSI'] = float(i)/float(j)*2.0*np.pi
            try:
                EF = Electricfield(para)
            except ValueError as detail:
                print 'encounter value error!',i,j,detail
            m = np.round(para['PSI']*EF.period/(2.0*np.pi))
            error = para['PSI']*EF.period - m*2.0*np.pi
            if  error < 1e-11:
                result = 'OK %e' %(error)
            else:
                result = 'FAIL! %e' %(error)
                print i,j,EF.period,result

    
