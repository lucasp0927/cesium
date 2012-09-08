#!/usr/bin/python2
from electricfield import Electricfield
from matplotlib import pyplot as plt
import numpy as np

if __name__ == '__main__':
    para = {}
    para['Tr'] = 1.0/91.9262177e6
    para['mu_c'] = 351.72571850e12
    para['PSI'] = 3.0/4.0*2.0*np.pi
    para['E_0'] = 1.0
    para['tao'] = 3e-14
    def test_period(para):
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
    #test_period(para)

                    
    print 'plot'
    print para['Tr']
    EF = Electricfield(para)
    w =  EF.time_no_field
    print w
    print EF.period
    print EF.Tr
    print EF.full_length
#    t = np.linspace(4.5*para['Tr']-w,4.5*para['Tr']+w,100000)
    t = np.linspace(0,EF.full_length,1000000)
    x = []
    
    for t_i in t:
        x.append(EF.field(t_i))
    figure = plt.plot(t,x)
    plt.show()
    
