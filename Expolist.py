#!/usr/bin/python2
from Expo import Expo
from copy import copy

class Expolist(object):
    """
    list of Expo object, that can be add, multiply...
    """
    def __init__(self,list=[Expo(0,0)]):
        """
        """
        self.terms = copy(list)
        self.simplify()

    def __repr__(self):
        tmp_str = '('
        for term in self.terms:
             tmp_str += term.__str__()+','
        tmp_str += ')'
        return tmp_str

    def __add__(self,expolist):
        result =  Expolist(self.terms+expolist.terms)
        result.simplify()
        return result

    def __sub__(self,expolist):
        expolist2 = expolist*-1
        result =  Expolist(self.terms+expolist2.terms)
        result.simplify()
        return result

    def __mul__(self,expolist):
        result_list = []
        try:
            for i in self.terms:
                for j in expolist.terms:
                    result_list.append(i*j)
        except AttributeError:
            #int or float
            for i in self.terms:
                result_list.append(i*expolist)
        result = Expolist(result_list)
        result.simplify()
        return result

    def simplify(self):
        #put all terms with same freq together
        freq_dict={}
        for i in range(len(self.terms)):
            if (self.terms[i].freq in freq_dict) == False:
                freq_dict[self.terms[i].freq]=[i]
            else:
                freq_dict[self.terms[i].freq].append(i)
        new_terms = []
        for freq in freq_dict:
            mag = 0
            for i in freq_dict[freq]:
                mag += self.terms[i].mag
            new_terms.append(Expo(mag,freq))
        #delete all terms with magnitude 0
        new_terms2 = filter(lambda x: x.mag!=0 and True or False,
                        new_terms)
        self.terms = copy(new_terms2)
        
    def RWA(self):
        """
        Rotating wave approximation
        """
        #remove terms with freq != 0
        new_terms = filter(lambda x: x.freq==0 and True or False,
                        self.terms)
        self.terms = copy(new_terms)

    def mag(self):
        self.simplify()
        if len(self.terms) == 0:
            return 0
        else:
            return self.terms[0].mag
        
if __name__ == '__main__':
    matrix = [Expolist() for i in range(2)]
    print Expolist([Expo(1,0)])
    matrix[0] = matrix[0]+Expolist([Expo(1,0)])
    print matrix
    a1 = Expolist([Expo(1,0)])
    a2 = Expolist()
    print a1
    print a2
    print a1*a2
