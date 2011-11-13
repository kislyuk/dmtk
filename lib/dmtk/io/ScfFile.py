import os
import sys
import copy

from struct import *

class ScfFile:
    def __init__(self, fileName):
        self.fileName = fileName
        self.numSamples = 0
        self.numBases = 0
        self.sampleSize = 2
        self.comments = [ "MACH=Astro" ]
        self.samples = []
        self.bases = []
        
    def write(self):
        outfile = open( self.fileName, 'wb' )
        self._writeHeaders(outfile)
        self._writeSamples(outfile)
        self._writeBases(outfile)
        self._writeComments(outfile)
        outfile.close()
        
    def _writeHeaders(self, outfile):
        #
        outfile.write( pack('>cccc','.','s','c','f') )
        outfile.write( pack('>III',self.numSamples,self._computeSampleOffset(),\
                            self.numBases) )
        # unused
        outfile.write( pack('>II',0,0) )
        #
        outfile.write( pack('>III',self._computeBasesOffset(), self._computeCommentsSize(),
                            self._computeCommentsOffset()) )
        # version format
        outfile.write( pack('>cccc','3','.','0','0') )
        outfile.write( pack('>IIII',self.sampleSize,0,0,self._computePrivateOffset()) )
        
        for i in xrange(18):
            outfile.write( pack('>I', 0) )
            
        print >>sys.stderr, "sample offset = ", self._computeSampleOffset()
        print >>sys.stderr, "bases offset = ", self._computeBasesOffset()
        print >>sys.stderr, "comments size = ", self._computeCommentsSize()
        print >>sys.stderr, "comments offset = ", self._computeCommentsOffset()
        print >>sys.stderr, "private offset = ", self._computePrivateOffset()
        
    def _writeSamples(self, outfile):
        print >>sys.stderr, "Writing samples, starting at offset=%d" % outfile.tell()
        for i in xrange(4):
            scratch = copy.copy( self.samples[i] )
            # 1st-level diff
            p_delta = 0
            for i in xrange(len(scratch)):
                p_sample = scratch[i]
                scratch[i] -= p_delta
                p_delta = p_sample
            # 2nd-level diff
            p_delta = 0
            for i in xrange(len(scratch)):
                p_sample = scratch[i]
                scratch[i] -= p_delta
                p_delta = p_sample
                
            for n in scratch:
                outfile.write( pack('>h',n) )
        
    def _writeBases(self, outfile):
        print >>sys.stderr, "Writing bases, starting at offset=%d" % outfile.tell()
        for base in self.bases:
            outfile.write( pack('>I',base.sampleIndex) )
        for i in xrange(4):
            for base in self.bases:
                outfile.write( pack('>b',base.qv[i]) )
        for base in self.bases:
            outfile.write( pack('>c',base.baseName) )
        # pad
        for i in xrange(3):
            for j in xrange(len(self.bases)):
                outfile.write( pack('>b',0) )
       
    def _writeComments(self, outfile):
        print >>sys.stderr, "Writing comments, starting at offset=%d" % outfile.tell()
        commentStr = "\n".join(self.comments) + '\n'
        outfile.write( pack('>%ds'%len(commentStr),commentStr) )
        
    def _computeSampleOffset(self):
        return 128
        
    def _computeBasesOffset(self):
        return self._computeSampleOffset() + self.numSamples * self.sampleSize * 4
    
    def _computeCommentsSize(self):
        cs = 0
        for comment in self.comments:
            cs += len(comment)+1
        return cs
    
    def _computeCommentsOffset(self):
        return self._computeBasesOffset() + self.numBases * 12
        
    def _computePrivateOffset(self):
        return self._computeCommentsOffset() + self._computeCommentsSize()
        
    def setTraces(self, channelData):
        self.samples = channelData
        self.numSamples = len(self.samples[0])
    
    def setBases(self, basecalls, basePositions):
        self.bases = []
        for b, bpos in zip(basecalls, basePositions):
            base = ScfBase(b,bpos)
            self.bases.append(base)
        self.numBases = len(self.bases)
    
class ScfBase:
    def __init__(self,baseName,bpos):
        self.sampleIndex = bpos
        self.qv = [0,0,0,0]
        self.baseName = baseName
        self._calcQv()
     
    def _calcQv(self):
        index = 'ACGT'.find(self.baseName)
        self.qv[index]=40
    
if __name__=='__main__':
    scf = ScfFile( sys.argv[1] )
    samples = [ \
               [ 5, 10, 15, 20, 15, 10, 5 ], \
               [ 6, 10, 15, 25, 15, 10, 6 ], \
               [ 7, 10, 15, 26, 15, 10, 7 ], \
               [ 8, 10, 15, 27, 15, 10, 8 ] ]
    basecalls = 'A'
    bpos = [ 3 ]
    scf.setTraces( samples )
    scf.setBases( basecalls, bpos )
    scf.write()