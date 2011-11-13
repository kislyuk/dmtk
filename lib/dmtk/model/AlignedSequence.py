__doc__=""" """
import sys

from dmtk.model.Range import Range

class AlignedSequence:
    GAP_CHAR = '-'
    



    def __init__(self, sequence=None, aligned=False, keepTermGaps=False, seqName=None, isReference=False ):
      
        self.seqName = seqName
        self.isReference = isReference
      
        if not aligned:
            self.sequence = sequence
        else:
            self.sequence = None
        # 0-based and [start,end) on the alignment
        self.alignedRange = Range()
        # gap positions are relative to start
        self.gaps = []
        if aligned:
            if keepTermGaps:
                self.setFromAlignment( sequence, start=0, end=len(sequence) )
            else:
                self.setFromAlignment( sequence )
            
    def clone(self):
        c = AlignedSequence()
        c.sequence = self.sequence
        c.alignedRange = Range(self.alignedRange.start,self.alignedRange.end)
        c.gaps = [ g for g in self.gaps ]
        return c
            
    def __getitem__( self, i ):
        """Returns the character at position i in the aligned sequence.
        Position 0 is the first character of this sequence.
        TODO: should be sped up with binary search."""
        d = 0
        oldg = -100
        for g in self.gaps:
            if g+d>=i:
                break
            d += 1
            oldg = g
        if oldg+d==i:
            return AlignedSequence.GAP_CHAR
        return self.sequence[ i-d ]
            
    def getAlignedStart(self):
        return self.alignedRange.start

    def getAlignedEnd(self):
        return self.alignedRange.end
        
    def setAlignedStart(self, value):
        self.alignedRange.start = value

    def setAlignedEnd(self, value):
        self.alignedRange.end = value
        
    def setFromAlignment(self, alignedSequence, start=-1, end=-1):
        if not self.sequence: 
            self.sequence = alignedSequence.replace(AlignedSequence.GAP_CHAR,'')
        if start<0 or end<0:
            self.alignedRange = self.__findAlignedStartStop( alignedSequence )
        else:
            self.alignedRange = Range( start, end )
        iUnalign = 0
        for i in self.alignedRange:
            ch = alignedSequence[i]
            if ch=='-':
                self.gaps.append( iUnalign-1 )
            else:
                iUnalign += 1
        
    def mergeAlignedSequence(self, alignedSequence):
        """ Adds gaps from an alternative alignment of this sequence. """
        iGap = 0
        newGaps = []
        diffs1, diffs2 = [], []
        # print 'diffing gaps=%s vs gaps2=%s' % ( str(self.gaps), str(alignedSequence.gaps) )
        for aGap in alignedSequence.gaps:
            truePos = aGap + alignedSequence.alignedRange.start
            if truePos>=self.alignedRange.start and truePos<self.alignedRange.end:
                myPos = truePos - self.alignedRange.start
                while iGap<len(self.gaps) and self.gaps[iGap]<myPos:
                    newGaps.append( self.gaps[iGap] )
                    diffs2.append( self.gaps[iGap] )
                    iGap += 1
                if iGap<len(self.gaps) and self.gaps[iGap]==myPos:
                    newGaps.append( self.gaps[iGap] )
                    iGap += 1
                    continue
                diffs1.append( myPos )
                newGaps.append( myPos )
            if truePos>=self.alignedRange.end:
                break
        while iGap<len(self.gaps):
            diffs2.append( self.gaps[iGap] )
            newGaps.append( self.gaps[iGap] )
            iGap += 1
        self.gaps = newGaps
        #print >>sys.stderr, 'diffs1 = ' + str(diffs1)
        #print >>sys.stderr, 'diffs2 = ' + str(diffs2)
        return diffs1, diffs2
    
    def toAlignedGaps(self, gaps, oldGaps=None):
        if len(gaps)==0: return []
        if not oldGaps:
            oldGaps = self.gaps
        j = 0
        delta = 0
        alignedGaps = []
        #print 'toAlignedGaps, gaps=%s, start=%d' % (str(gaps),self.start)
        #print '  oldGaps = ' , oldGaps
        #print >>sys.stderr, "gaps = " + str(gaps)
        #print >>sys.stderr, "oldGaps = " + str(oldGaps)
        for gap in gaps:
            while j<len(oldGaps) and oldGaps[j]<gap-self.alignedRange.start:
                delta += 1
                j += 1
            #print 'delta = %d' % delta
            aGap = gap + delta - self.alignedRange.start
            #print 'aGap = %d' % aGap
            alignedGaps.append( aGap )
            #if aGap>=0 and not old: delta += 1
        return alignedGaps
        
    def insertGaps(self, alignedGaps, offset=0):
        if len(alignedGaps)==0: return
        #print >>sys.stderr, "inserting " + str(alignedGaps)
        #
        # insert gaps which begin before we do
        # (these bump our start point up)
        #
        k = 0
        deltaStart = 0
        while k<len(alignedGaps) and alignedGaps[k]+offset<self.alignedRange.start:
            deltaStart += 1
            k += 1
        if k==len(alignedGaps):
            self.alignedRange.addDelta( deltaStart )
            return
        #print >>sys.stderr, "Now k=%d"%k
        #
        # consider old gaps in the aligned frame
        #
        myGaps = []
        oldAlignedGaps = []
        delta = 1
        for gap in self.gaps:
            oldAlignedGaps.append( gap + delta )
            delta += 1
        #print >>sys.stderr, "oldAlignedGaps = " + str(oldAlignedGaps)
        #
        # add new gaps, with bookkeeping for the old frame
        #
        j = 0
        delta = 0
        #print 'oldAlignedGaps = ', oldAlignedGaps
        #print 'alignedGaps = ', alignedGaps[k:]
        for gap in alignedGaps[k:]:
            aGap = gap - self.alignedRange.start + offset
            #print 'aGap = ', aGap
            while j<len(oldAlignedGaps) and oldAlignedGaps[j]<=aGap:
                delta += 1
                myGaps.append( oldAlignedGaps[j]-delta )
                j += 1
            #print 'delta = ', delta
            #print 'adding ', str(aGap-delta)
            myGaps.append( aGap-delta )
        #print >>sys.stderr, 'myGaps1 = %s' % str(myGaps)
        #
        # add in any old gaps at the end which we missed
        #
        while j<len(oldAlignedGaps):
            delta += 1
            myGaps.append( oldAlignedGaps[j]-delta )
            j += 1
        #print >>sys.stderr, 'myGaps2 = %s' % str(myGaps)
        self.alignedRange.addDelta( deltaStart )
        self.gaps = myGaps
        
    def __findAlignedStartStop( self, aSequence ):
        """ Returns Range for the beginning and end 
        of the sequence in this aligned sequence.
        End points are 0-based and [start,end)
        """
        start = len(aSequence)-1
        for i, c in enumerate(aSequence):
            if c!='-':
                start = i
                break
        end = 0
        for i in xrange(len(aSequence)):
            j = len(aSequence)-i-1
            c = aSequence[j]
            if c!='-':
                end = j
                break
        return Range( start, end+1 )
        
    def getFullString(self):
        return "-"*self.alignedRange.start + str(self)
        
    def __str__(self):
        return "".join( [ ch for ch in self ] )

    def __len__(self):
        "returns length of aligned sequence"
        return len(self.gaps)+len(self.sequence)
        
    def __iter__(self):
        iGap = 0
        for i, ch in enumerate(self.sequence):
            yield ch
            while iGap<len(self.gaps) and self.gaps[iGap]==i:
                yield AlignedSequence.GAP_CHAR
                iGap += 1
