__doc__="""Simple storage class for modeling the 4C HMM parameters with I/O
support for Susan's format."""
import os
import numpy as n

class HmmParams:
    def __init__( self, fileName ):
        self.dark = n.zeros(4)
        self.stick = n.zeros(4)
        self.branch = n.zeros(4)
        self.miscall = n.zeros((4,4))
        self.__load(fileName)

    def __load( self, fileName ):
        infile = open( fileName, 'r' )
        i, j, state = 0, 0, ''
        for line in infile:
            if 'Dark' in line:
                state = 'dark'
                i = 0
                continue
            if 'Stick' in line:
                state = 'stick'
                i = 0
                continue
            if 'Branch' in line:
                i = 0
                state = 'branch'
                continue
            if 'Miscall' in line:
                state = 'miscall'
                i = 0
                j = 0
                continue
            line = line.rstrip()
            if state=='dark':
                self.dark[i] = float(line)
                i += 1
            elif state=='stick':
                self.stick[i] = float(line)
                i += 1
            elif state=='branch':
                self.branch[i] = float(line)
                i += 1
            elif state=='miscall':
                self.miscall[ j, i ] = float(line)
                i += 1
                if i==4:
                    j += 1
                    i = 0
        infile.close()

    def write( self, file ):
        print >>file, '>Dark'
        for v in self.dark:
            print >>file, v
        print >>file, '>Stick'
        for v in self.stick:
            print >>file, v
        print >>file, '>Branch'
        for v in self.branch:
            print >>file, v
        print >>file, '>Miscall'
        for v in self.miscall:
            for v2 in v:
                print >>file, v2

    def rescale( self, epsilon ):
        """Rescale the error represented by this model to give
        roughly epsilon worth of error (on average).  HMM error is
        non-additive, so we can only approximate this tranformation."""
        sumError = sum(self.dark) + sum(self.stick) + sum(self.branch)
        sumError += 4.0 - sum( [ self.miscall[i,i] for i in xrange(4) ] )
        rescale = epsilon / ( sumError / 4.0 )
        self.dark = rescale * self.dark
        self.stick = rescale * self.stick
        self.branch = rescale * self.branch
        for i in xrange(4):
            for j in xrange(4):
                if i==j: continue
                self.miscall[i,j] *= rescale
        for i in xrange(4):
            self.miscall[i,i] = 1.0 - rescale * ( 1.0 - self.miscall[i,i] )

