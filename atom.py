#!/usr/bin/python2
from __future__ import division
import math
import numpy as np


class Atom():
    """
    Calculate 3j 6j and dipole matrix elements.
    """

    def __init__(self, ):
        """
        """

    def bad_value (self,j1,j2,j3,m1,m2,m3):
        """
        Check
        """
        if (j1<np.abs(j2-j3) or j1 > (j2+j3)):
            return True
        if (np.abs(m1)>j1 or np.abs(m2)>j2 or np.abs(m3)>j3): # add comment
            return True
        if (m1+m2+m3 != 0):
            return True
        return False

    def threej (self,j1,j2,j3,m1,m2,m3):
        if self.bad_value(j1,j2,j3,m1,m2,m3):
            return 0

        jphase = math.pow(-1.0,j1-j2-m3)
        fac=np.zeros(10)

        fac[0]=math.factorial(j1+j2-j3)
        fac[1]=math.factorial(j1-j2+j3)
        fac[2]=math.factorial(-j1+j2+j3)
        fac[3]=math.factorial(j1+m1)
        fac[4]=math.factorial(j1-m1)
        fac[5]=math.factorial(j2+m2)
        fac[6]=math.factorial(j2-m2)
        fac[7]=math.factorial(j3+m3)
        fac[8]=math.factorial(j3-m3)
        fac[9]=math.factorial(j1+j2+j3+1)

        jprodfac = math.sqrt(1.0/fac[9])

        for i in range(9):
            jprodfac*= math.sqrt(fac[i])

        kmax = np.min([j1+j2-j3 , j1-m1 , j2+m2]) # DIFFRENT
        kmin = np.max([0.0 , -(j3-j2+m1) , -(j3-j1-m2)])
        jsum = 0.0

        for k in range(int(kmin),int(kmax+1)):
            jsum += math.pow(-1,k)*1.0/(math.factorial(k)*math.factorial(j3-j2+m1+k)*math.factorial(j3-j1-m2+k)*math.factorial(j1+j2-j3-k)*math.factorial(j1-k-m1)*math.factorial(j2+m2-k))
        return jphase*jprodfac*jsum

    def bad_value6j (self,j1,j2,j3,l1,l2,l3):
        if (j1<np.abs(j2-j3) or j1>(j2+j3)):
            return True
        if (j1<np.abs(l2-l3) or l1>(l2+l3)):
            return True
        if (l1<np.abs(j2-l3) or l1>(j2+l3)):
            return True
        if (l1<np.abs(l2-j3) or l1>(l2+j3)):
            return True
        return False

    def delta(self,a,b,c):
        x=math.factorial(a+b-c)*math.factorial(a-b+c)*math.factorial(-a+b+c)
        return math.sqrt(float(x)/float(math.factorial(a+b+c+1)))

    def sixj (self,j1,j2,j3,l1,l2,l3):
        if self.bad_value6j(j1,j2,j3,l1,l2,l3):
            return 0

        jphase = math.pow(-1.0,j1+j2+l1+l2)
        proddelt=self.delta(j1,j2,j3)*self.delta(l1,l2,j3)*self.delta(l1,j2,l3)*self.delta(j1,l2,l3)

        kmax = np.min([j1+j2+l1+l2+1,j1+j2-j3,l1+l2-j3,j1+l2-l3,l1+j2-l3])
        kmin = np.max([0.0,j1+l1-j3-l3,j2+l2-j3-l3])
        jsum = 0
        jsfac = np.zeros(8)

        for k in range(int(kmin),int(kmax+1)):
            jsfac[0] = math.factorial(j1+j2+l1+l2+1-k)
            jsfac[1] = math.factorial(k)
            jsfac[2] = math.factorial(j1+j2-j3-k)
            jsfac[3] = math.factorial(l1+l2-j3-k)
            jsfac[4] = math.factorial(j1+l2-l3-k)
            jsfac[5] = math.factorial(l1+j2-l3-k)
            jsfac[6] = math.factorial(-j1-l1+j3+l3+k)
            jsfac[7] = math.factorial(-j2-l2+j3+l3+k)
            product = 1.0

            for j in range(1,8):
                product /= jsfac[j]
            jsum += math.pow(-1.0,k)*jsfac[0]*product
        return jphase*proddelt*jsum

    def dipole_element (self,q,L1,L2,F1,F2,mf1,mf2,J1,J2,I):
        """
        expressed as multiples of <j||er||j'>
        1 is S
        2 is P
        """
        #<J|er|J'>
        if J1 == 0.5 and J2 == 0.5:
            er = 2.7020e-29
        elif J1 == 0.5 and J2 == 1.5:
            er = 3.8014e-29
        elif J1 == 1.5 and J2 == 0.5:
            er = 3.8014e-29
        else:
            print 'J1 or J2 is incorrect'
        cg = self.cg_coef(q,L1,L2,F1,F2,mf1,mf2,J1,J2,I)
        # <F|er|f'>
        fer = math.pow(-1,F2+J1+1+I)*math.sqrt((2*F2+1)*(2*J1+1))*self.sixj(J1,J2,1,F2,F1,I)
        return (cg*fer*er if (L1 != L2) else 0) #this is the formula 34 on D line data

    def cg_coef (self,q,L1,L2,F1,F2,mf1,mf2,J1,J2,I):
        return (math.pow(-1,F2-1+mf1)*math.sqrt(2.0*F1+1.0)*self.threej(F2,1.0,F1,mf2,-q,-mf1) if (L1 != L2) else 0)

if __name__ == "__main__":
    a=Atom()
    for q in (-1,0,1):
        print "q=%i" %q
        for f in (4,3,2):
            ans = []
            print "f=%i"%f
            for m in range(-3,4):
                coef1 = {'q':q,
                 'L1':0,
                 'L2':1,
                 'F1':3,
                 'F2':f,
                 'mf1':m,
                 'mf2':m+q,
                 'J1':1.0/2.0,
                 'J2':3.0/2.0,
                 'I':7.0/2.0}
                ans.append(a.dipole_element(**coef1)/3.8014e-29)
            print "%s \n" %ans
        print "\n"
