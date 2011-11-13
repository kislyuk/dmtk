class BaseLink:
    """Abstract base class for links."""
    
    def __init__(self): pass
    def getFirstId(self):       raise NotImplementedError
    def getSecondId(self):      raise NotImplementedError
    def getFirstStrand(self):   raise NotImplementedError
    def getSecondStrand(self):  raise NotImplementedError
    def getFirstLength(self):   raise NotImplementedError
    def getSecondLength(self):  raise NotImplementedError
    def getMinSpan(self):       raise NotImplementedError
    def getMaxSpan(self):       raise NotImplementedError
    def getTypeString(self):    raise NotImplementedError

class SequenceLink(BaseLink):
    """A link derived from a split or shredded sequence"""
    
    def __init__(self, firstId, secondId, span=0, firstLength=0, secondLength=0, firstStrand="+", secondStrand="+"):
        self.firstId = firstId
        self.secondId = secondId
        self.span = span
        self.firstLength  = firstLength
        self.secondLength = secondLength
        self.firstIsRepeat  = False
        self.secondIsRepeat = False
        self.firstStrand = firstStrand
        self.secondStrand = secondStrand

    def getFirstId(self):
        return self.firstId
    
    def getSecondId(self):
        return self.secondId

    def getFirstStrand(self):
        return self.firstStrand
    
    def getSecondStrand(self):
        return self.secondStrand

    def getFirstLength(self):
        return self.firstLength
    
    def getSecondLength(self):
        return self.secondLength

    def getMinSpan(self):
        return self.span
        
    def getMaxSpan(self):
        return self.span
    
    def __str__(self):
        return '<link type="%s" firstId="%s" secondId="%s" firstStrand="+" secondStrand="+" minSpan="%i" maxSpan="%i"/>' \
            % (self.getTypeString(), self.firstId, self.secondId, self.span, self.span)

    def getTypeString(self):
        return "SequenceLink"
        
class GenomeLink(BaseLink):
    """ A link derived from two hits to the reference genome. Used for creating exactly correct test case scaffolds."""

    def __init__(self, firstHit, secondHit, minSpan=None, maxSpan=None, tolerance=100):
        self.tolerance = tolerance
        self.firstHit  = firstHit
        self.secondHit = secondHit
        self.minSpan = minSpan
        self.maxSpan = maxSpan
        self.ensureSpans()
    
    def getFirstId(self):
        return self.firstHit.query_id
    
    def getSecondId(self):
        return self.secondHit.query_id

    def getFirstStrand(self):
        return self.firstHit.query_strand
    
    def getSecondStrand(self):
        return self.secondHit.query_strand

    def getMinSpan(self):
        return self.minSpan
        
    def getMaxSpan(self):
        return self.maxSpan
    
    def ensureSpans(self):
        if self.minSpan == None or self.maxSpan == None:
            self.setSpanFromTarget()

    def getFirstLength(self):
        return self.firstHit.query_length
    
    def getSecondLength(self):
        return self.secondHit.query_length

    def setSpanFromTarget(self):
        # TODO make sure they are fully aligned
        midpoint = (self.secondHit.target_start - self.firstHit.target_end) / 2
        self.minSpan = midpoint - self.tolerance
        self.maxSpan = midpoint + self.tolerance

    def getTypeString(self):
        return "GenomeLink"
     

    
# TODO possibly rename to HitLink?
class BambusLink(BaseLink):
    """ A link derived from two hits, either from the same sequence or two related sequences (e.g. strobes)"""

    def __init__(self, firstHit, secondHit, minSpan=None, maxSpan=None, tolerance=100):
        self.tolerance = tolerance
        self.firstHit  = firstHit
        self.secondHit = secondHit
        self.minSpan = minSpan
        self.maxSpan = maxSpan
        self.ensureSpans()
    
    def getFirstId(self):
        return self.firstHit.query_id
    
    def getSecondId(self):
        return self.secondHit.query_id

    def getMinSpan(self):
        return self.minSpan
        
    def getMaxSpan(self):
        return self.maxSpan
    
    def ensureSpans(self):
        # TODO make sure they are non-overlapping
        if self.minSpan == None or self.maxSpan == None:
            self.setSpanFromQuery()
            # self.minSpan = 0
            # self.maxSpan = self.tolerance

    def getFirstLength(self):
        return self.firstHit.query_length
    
    def getSecondLength(self):
        return self.secondHit.query_length

    def setSpanFromQuery(self):
        innerSpan = self.secondHit.query_start - self.firstHit.query_end
        if innerSpan < 0:
            # we have hits that overlap on the query. Bambus will not like this.
            adjust = abs(innerSpan) / 2 + 1

            self.firstHit.query_end    -= adjust
            if self.firstHit.query_strand == "+": self.firstHit.target_end   -= adjust
            else:                                 self.firstHit.target_start += adjust

            self.secondHit.query_start += adjust
            if self.secondHit.query_strand == "+": self.secondHit.target_start += adjust
            else:                                  self.secondHit.target_end   -= adjust

        outerSpan = self.secondHit.query_end - self.firstHit.query_start
        self.minSpan = max(outerSpan - self.tolerance, 0)
        self.maxSpan = outerSpan + self.tolerance

    def getTypeString(self):
        return "BambusLink"
     

