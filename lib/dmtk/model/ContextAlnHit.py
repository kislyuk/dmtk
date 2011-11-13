#!/usr/bin/env python
"""The ContextAlnHit class extends CmpH5AlnHit, providing a simple interface for
analyzing data on a per-alignment-position basis. It is meant to be lightweight
and useful."""

from dmtk.io.cmph5.CmpH5AlnHit import CmpH5AlnHit

GAP = '-'

class ContextAlnHit( CmpH5AlnHit ):
    
    def __init__( self, sourceAlnHit=None ):
        """Can be initialized by another alnHit"""
        if sourceAlnHit != None:
            self.__dict__ = sourceAlnHit.__dict__.copy()
        else:
            CmpH5AlnHit.__init__( self )

        # Until Pat's Realigner bug is fixed...
        # self._favorBranching( )

    def _favorBranching( self ):
        """Slides insertions within homopolymers to the left"""
        ref = list(self.alignedTarget)
        read = self.alignedQuery
        changed = True
        while( changed ):
            changed = False
            for i in range(len(ref)):
                if ref[i] == GAP and ref[i-1] == read[i] and i-1>0:
                    ref[i] = ref[i-1]
                    ref[i-1] = GAP
                    changed = True

        self.alignedTarget = "".join( ref )

    def __iter__( self ):
        """Iterates over AlnPosn objects, which represent
        locations witin a pairwise alignment and provide
        relevant information."""
        for posn in range(len(self.alignedQuery)):
            yield AlnPosn( self, posn )

    @property
    def aln2ref( self ):
        """Caches a mapping of aln coordinate to reference position."""
        if not hasattr( self, "_aln2ref" ):
            self._aln2ref = {}
            rc = self.target_strand == "-"
            refPosn = self.target_end-1 if rc else self.target_start
            increment = -1 if rc else 1
            for i in range(len(self.alignedTarget)):
                self._aln2ref[ i ] = refPosn
                if self.alignedTarget[ i ] != GAP:
                    refPosn += increment
        return self._aln2ref

class AlnPosn( object ):
    """Provides a useful interface for single alignment positions"""
    
    def __init__( self, parent, posn ):
        """Takes in the parent AlignmentHit and the position"""
        self._parent    = parent
        self._posn      = posn
        self._reference = parent.alignedTarget
        self._read      = parent.alignedQuery

    def __getitem__( self, value ):
        """String indices return the corresponding pulse metric
        from the parent alignmentHit for this position."""
        if isinstance( value, str ):
            return self._parent[ value ][ self._posn ]

    @property
    def ref( self ):
        """Returns the reference base associated with this position."""
        return self._reference[ self._posn ]

    @property
    def read( self ):
        """Returns the read base associated with this position."""
        return self._read[ self._posn ]

    @property
    def refPosn( self ):
        """Returns the reference position of this AlnPosn."""
        return self._parent.aln2ref[ self._posn ]

    @property
    def errorType( self ):
        """Returns one of 'MM', 'Del', 'Ins', or None."""
        read, ref = self.read, self.ref
        if read == GAP:
            return 'Del'
        if ref == GAP:
            return 'Ins'
        if read != ref:
            return 'MM'
        return None

    def refContext( self, left=0, right=0 ):
        """Given the number of additional positions to include to the left
        and right of this position (with respect to the reference), 
        returns the local reference context."""
        if self.ref == GAP:
            right += 1
        basesLeft  = self._refBases( left, -1 )
        basesRight = self._refBases( right, 1 )
        if basesLeft == None or basesRight == None:
            return None
        return basesLeft + self.ref.replace(GAP,"") + basesRight
        
    def _refBases( self, number, direction ):
        """Given the number of desired bases and a direction
        to move along the alignment (-1,+1), returns the relevant
        reference bases."""
        curPosn, curBases = self._posn, 0
        while curBases < number:
            curPosn += direction
            if curPosn < 0 or curPosn >= len(self._reference):
                return None
            if self._reference[ curPosn ] != GAP:
                curBases += 1
        if direction == 1:
            return self._reference[ self._posn+1: curPosn+1 ].replace(GAP,"")
        else:
            return self._reference[ curPosn: self._posn ].replace(GAP,"")
