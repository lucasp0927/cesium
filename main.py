#!/usr/bin/python2
from density_matrix import System

if __name__ ==  '__main__':
    an=3
    aomega = [351E12,9E9,0]
    adipole=[[0,1000000,1000000],
            [1000000,0,0],
            [1000000,0,0]]
    #remember to /2
    anu=[351E12-9E9,351E12] # on resonence
    aup = [0]
    alow1 = [1]
    alow2 = [2]
    #decoherence
    aGamma1 = 5000000
    aGamma12 = 2500000
    aGamma13 = 2500000
    agamma1 = 10000
    agamma2 = 10000

    system = System(an,aomega,adipole,anu,aup,alow1,alow2,aGamma1,aGamma12,aGamma13,agamma1,agamma2)
    system.sweep(-1E7,1E7,1000)
