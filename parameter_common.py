#!/usr/bin/python2
from __future__ import division
import numpy as np
from atom import Atom
class Parameter(object):
    """
    """
    
    def __init__(self,l1f,l0f,d1,egpair,omega_list,parameter,filename):
        """
        if d1 = 1 then d1 else d2
        l1f = (5,4,3)
        l0f = (4,3)
        """
        self.d1 = d1
        self.egpair = egpair
        self.l1f = l1f
        self.l0f = l0f
        self.parameter = parameter
        self.omega_list = omega_list
        self.filename = filename
        n = 0
        for i in l1f:
            n += 2*i + 1
        for i in l0f:
            n += 2*i + 1
        self.parameter['n'] = n
        
    def level_group(self):
        level_group = []
        l1_levels = 0
        for i in self.l1f:
            l1_levels += 2*i+1
        level_group.append(range(l1_levels))
        for i in self.l0f:
            level_group.append(range(l1_levels,l1_levels+2*i+1))
            l1_levels += 2*i+1
        self.parameter['level_group']=level_group

    def omega(self):
        counter = 0
        self.parameter['omega']=[]
        for l in (self.l1f,self.l0f):
            for i in l:
                for j in range(2*i+1):
                    self.parameter['omega'].append(self.omega_list[counter])
                counter += 1

    def index2lfm(self,n):
        l = 1
        for i in self.l1f:
            if n < 2*i + 1:
                f = i
                m = n - f
                break
            else:
                n -= 2*i + 1
        else:
            l = 0
            for j in self.l0f:
                if n < 2*j + 1:
                    f = j
                    m = n - f
                    break
                else:
                    n -= 2*j + 1
        return l,f,m

    def lfm2index(self,l,f,m):
        index = 0
        if l == 1:
            for i in self.l1f:
                if f != i:
                    index += 2*i+1
                else:
                    index += m + f
                    break
        else:
            for i in self.l1f:
                index += 2*i+1
            for i in self.l0f:
                if f != i:
                    index += 2*i+1
                else:
                    index += m + f
                    break
        return index
    
    def dipole(self):
        if self.d1 == 1:
            j2 = 1.0/2.0
        else:
            j2 = 3.0/2.0
        n=self.parameter['n']
        self.parameter['dipole'] = [[0 for i in range(n)] for j in range(n)]
        cs = Atom()
        for i in range(n):
            for j in range(n):
                d1 = self.index2lfm(i)
                d2 = self.index2lfm(j)
                if d1[0] == 0 and d2[0] == 1:
                    coef = {'q':0,
                             'L1':0,
                             'L2':1,
                             'F1':d1[1],
                             'F2':d2[1],
                             'mf1':d1[2],
                             'mf2':d2[2],
                             'J1':1.0/2.0,
                             'J2':j2,
                             'I':7.0/2.0}
                    self.parameter['dipole'][i][j] = cs.dipole_element(**coef)
                elif d2[0] == 0 and d1[0] == 1:
                    coef = {'q':0,
                             'L1':0,
                             'L2':1,
                             'F1':d2[1],
                             'F2':d1[1],
                             'mf1':d2[2],
                             'mf2':d1[2],
                             'J1':1.0/2.0,
                             'J2':j2,
                             'I':7.0/2.0}
                    self.parameter['dipole'][i][j] = cs.dipole_element(**coef)
                else:
                    self.parameter['dipole'][i][j] = 0.0

    def decoherence(self):
        if self.d1 == 1:
            j2 = 1.0/2.0
        else:
            j2 = 3.0/2.0    
        Gamma = 2*np.pi*4.575e6 #this is parameter for D1 line
        n=self.parameter['n']
        self.parameter['decoherence_matrix'] = [[[] for i in range(n)] for j in range(n)]
        cs = Atom()
        """
        e l=1 f=3 g l=0 f=4
        e l=1 f=3 g l=0 f=3
        """
        for pair in self.egpair:
            for i in range(n):
                for j in range(i,n):
                    d1 = self.index2lfm(i)
                    d2 = self.index2lfm(j)
                    if d1[0:2] == pair[0] and d2[0:2] == pair[0] and i != j:
                         self.parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma])
                    if d1[0:2] == pair[0] and d2[0:2] == pair[1]:
                        self.parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
                    elif d1[0:2] == pair[1] and d2[0:2] == pair[0]:
                        self.parameter['decoherence_matrix'][i][j].append([i,j,-1.0*Gamma/2.0])
                    elif d1[0:2] == pair[1] and d2[0:2] == pair[1]:
                        for q in (-1.0,0.0,1.0):
                            f1 = pair[0][1]
                            if (d1[2]+q <= f1 and d1[2]+q >= -1*f1) and (d2[2]+q <= f1 and d2[2]+q >= -1*f1):
                                coef1 = {'q':q,
                                         'L1':0,
                                         'L2':1,
                                         'F1':pair[1][1],
                                         'F2':pair[0][1],
                                         'mf1':d1[2],
                                         'mf2':d1[2]+q,
                                         'J1':1.0/2.0,
                                         'J2':j2,
                                         'I':7.0/2.0}
                                coef2 = {'q':q,
                                         'L1':0,
                                         'L2':1,
                                         'F1':pair[1][1],
                                         'F2':pair[0][1],
                                         'mf1':d2[2],
                                         'mf2':d2[2]+q,
                                         'J1':1.0/2.0,
                                         'J2':j2,
                                         'I':7.0/2.0}
                                tmp = Gamma*cs.cg_coef(**coef1)*cs.cg_coef(**coef2)
                                if tmp != 0.0:
                                    ii = self.lfm2index(pair[0][0],pair[0][1],d1[2]+q)
                                    jj = self.lfm2index(pair[0][0],pair[0][1],d2[2]+q)
                                    self.parameter['decoherence_matrix'][i][j].append([ii,jj,tmp])
                                    if ii == jj:
                                        self.parameter['decoherence_matrix'][int(ii)][int(jj)].append([ii,jj,-1*tmp])                                
    def write(self):
        self.level_group()
        self.omega()
        self.dipole()
        self.decoherence()
        print "      L   F   M|"
        print "----------------------------------------"
        sum = 0.0
        psum = 0.0
        for i in range(self.parameter['n']):
            iter = self.parameter['decoherence_matrix'][i][i]
            for j in iter:
                sum += j[2]
                psum += j[2]
            print "{:3d}:{:3d} {:3d} {:3d}| {:<100} |Sum is: {:>20f}".format(i,self.index2lfm(i)[0],self.index2lfm(i)[1],self.index2lfm(i)[2],iter,psum)
            psum = 0.0        

        print "the sum is %f" %sum
        txtf = open(self.filename+'.txt','w')
        txtf.write(str(self.parameter))
        txtf.close()
        
if __name__ == '__main__':
    A=351.72571850e12
    B=12.79851e6
    C=263.8906e6
    D=188.4885e6
    E=399.7128e6
    F=4.021776399375e9
    G=5.170855370625e9
    para = Parameter((5,4,3,2),(4,3),0,
                     (((1,5),(0,4)),((1,5),(0,3)),((1,4),(0,4)),((1,4),(0,3)),((1,3),(0,4)),((1,3),(0,3)),((1,2),(0,4)),((1,2),(0,3))),
                     (G+A+C,G+A+B,G+A-D,G+A-E,G+F,0),
                     {'nu': [A-F,A+G],
                      'e_amp': [100, 100],
                      'sweep_profile':[0,-1E10,1E10,400]
                      },
                     'setting/d2')
    para.write()
