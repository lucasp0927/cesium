#!/usr/bin/python2
class Expo(object):
    """
    mag*Exp(i freq t)
    """
    
    def __init__(self,mag=0,freq=0):
        """
        """
        self.mag,self.freq = mag,freq

    def __mul__(self,expo2):
        try:
            mag = self.mag*expo2.mag
            freq = self.freq+expo2.freq
        except AttributeError:
            #int or float
            mag = self.mag*expo2
            freq = self.freq
        result = Expo(mag,freq)
        return result
        
    def __repr__(self):
        return '%s*Exp(i*%s*t)'%(str(self.mag),str(self.freq))

if __name__ == '__main__':
    a = Expo(1.0+2.0j,1.0)
    b = Expo(2.0,2.0)
    d=[a,b]
    c=a*b
    print a
    print b

        
