#!/usr/bin/python2
from __future__ import division
import numpy as np
from atom import Atom
import math
import pydot

class Parameter(object):
    """
    """
    
    def __init__(self,l1f,l0f,B,d1,egpair,omega_list,parameter,filename):
        """
        if d1 = 1 then d1 else d2
        l1f = (5,4,3)
        l0f = (4,3)
        """
        self.d1 = d1
        self.egpair = egpair
        self.l1f = l1f
        self.l0f = l0f
        self.B = B
        self.parameter = parameter
        self.omega_list = omega_list
        self.filename = filename
        self.graph = pydot.Dot('csgraph',graph_type = 'digraph')
        self.l1subg = {}
        self.l0subg = {}
        self.l1subn = {}
        self.l0subn = {}        
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
        mub = 9.27400915e-28/1.054571628e-34
        gs = 2.0023193043622
        gl = 0.99999587
        gjs = 2.00254032
        gjp12 = 0.665900
        gjp32 = 1.33400
        gi = -0.00039885395
        counter = 0
        I = 7.0/2.0
        self.parameter['omega']=[]
        for LL in enumerate((self.l1f,self.l0f)):
            if LL[0] == 0: #p
                L = 1
                if self.d1 == 1: #d1
                    gj = gjp12
                    J = 1.0/2.0
                elif self.d1 == 0: #d2
                    gj = gjp32
                    J = 3.0/2.0
            elif LL[0] == 1: #s
                gj = gjs
                L=0
                J = 1.0/2.0
                
            for F in LL[1]:#l
                gf = gj * (F*(F+1) - I*(I+1) + J*(J+1)) / (2*F*(F+1)) + gi * (F*(F+1)+ I*(I+1) - J*(J+1)) / (2*F*(F+1))
                for m in range(-F,F+1): #m
                    print L,F,m,self.omega_list[counter] + mub * self.B * m * gf
                    self.parameter['omega'].append(self.omega_list[counter] + mub * self.B * m * gf)
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
        self.parameter['dipole'] = []
        for k in range(len(self.parameter['e_amp'])):
            if self.d1 == 1:
                j2 = 1.0/2.0
            else:
                j2 = 3.0/2.0
            n=self.parameter['n']
            tmp = [[0 for i in range(n)] for j in range(n)]
            cs = Atom()
            for i in range(n):
                for j in range(n):
                    d1 = self.index2lfm(i)
                    d2 = self.index2lfm(j)
                    if d1[0] == 0 and d2[0] == 1:
                        q=self.parameter['e_amp'][k][1]
                        coef = {'q':q,
                                 'L1':0,
                                 'L2':1,
                                 'F1':d1[1],
                                 'F2':d2[1],
                                 'mf1':d1[2],
                                 'mf2':d2[2]+q,
                                 'J1':1.0/2.0,
                                 'J2':j2,
                                 'I':7.0/2.0}
                        tmp[i][j] = cs.dipole_element(**coef)
                    elif d2[0] == 0 and d1[0] == 1:
                        q=self.parameter['e_amp'][k][1]                        
                        coef = {'q':q,
                                 'L1':0,
                                 'L2':1,
                                 'F1':d2[1],
                                 'F2':d1[1],
                                 'mf1':d2[2],
                                 'mf2':d1[2]+q,
                                 'J1':1.0/2.0,
                                 'J2':j2,
                                 'I':7.0/2.0}
                        tmp[i][j] = cs.dipole_element(**coef)
                    else:
                        tmp[i][j] = 0.0
            self.parameter['dipole'].append(tmp)

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
                                #this correction coefficient (see equation 54) should be written to atom.py later
                                rev = (-1)**(pair[0][1]-pair[1][1]+q)*math.sqrt((2*pair[0][1]+1)/(2*pair[1][1]+1))
                                rev = rev**2
                                tmp = Gamma*cs.cg_coef(**coef1)*cs.cg_coef(**coef2)*rev
                                if tmp != 0.0:
                                    ii = self.lfm2index(pair[0][0],pair[0][1],d1[2]+q)
                                    jj = self.lfm2index(pair[0][0],pair[0][1],d2[2]+q)
                                    self.parameter['decoherence_matrix'][i][j].append([ii,jj,tmp])
                                    if ii == jj:
                                        self.parameter['decoherence_matrix'][int(ii)][int(jj)].append([ii,jj,-1*tmp])
                                        f1 = int(pair[0][1])
                                        f2 = int(pair[1][1])
                                        label = '%f'%tmp
                                        self.graph.add_edge(pydot.Edge(self.l1subn[f1][int(d1[2]+q+f1)],self.l0subn[f2][int(d1[2]+f2)],label = label))
                                        
    def write(self):
        self.prepare_graph()
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
        self.write_graph()
        
    def prepare_graph(self):
        for i in self.l1f:
            self.l1subg[i] = pydot.Subgraph('',rank = 'same')
            self.l1subn[i] = []
            for j in range(-1*i,i+1):
                name = 'L=1,F=%d,m=%d' %(i,j)
                self.l1subn[i].append(pydot.Node(name))
                self.l1subg[i].add_node(self.l1subn[i][j+i])
            self.graph.add_subgraph(self.l1subg[i])
        for i in self.l0f:
            self.l0subg[i] = pydot.Subgraph('',rank = 'same')
            self.l0subn[i] = []
            for j in range(-1*i,i+1):
                name = 'L=0,F=%d,m=%d' %(i,j)
                self.l0subn[i].append(pydot.Node(name))
                self.l0subg[i].add_node(self.l0subn[i][j+i])
            self.graph.add_subgraph(self.l0subg[i])        
    def write_graph(self):
        self.graph.write_png(self.filename+'png')
    
if __name__ == '__main__':
    pass
