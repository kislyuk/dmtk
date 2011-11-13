__doc__="""Model class for a structure which stores multiple aligned sequences."""

from itertools import izip, imap

from dmtk.model.AlignedSequence import AlignedSequence
from dmtk.model.Range import Range

class MSA:
    """Storage class for modeling a collection of AlignedSequences.
    """
    def __init__(self):
        
        self._msa_data = {}
        self._msa_data['sequences'] = {}
        # we need to store an ordered list of keys so iteration is
        # well defined for this object
        self._msa_data['sequence_names'] = []
        self._msa_data['alignment_length'] = 0
        self._msa_data['has_names'] = False
        self._msa_data['has_reference'] = False
        self._msa_data['reference_name'] = None
        self._msa_data['max_seq_id'] = -1 #need this ID as sequence IDs when names are not provided.
        
    def add( self, aSeq, name=None, strand='+', isReference=False ):
        if name: 
            self._msa_data['has_names'] = True
        else:
            name = str(self._msa_data['max_seq_id'])
            self._msa_data['max_seq_id'] += 1

        if isReference:
            if self.hasReference():
                raise Exception, "This MSA already contains a reference named:%s"% self.getReferenceSeqName()
            else:
                self.setReference(True)
                self.setReferenceSeqName(name)


        self._msa_data['sequences'][name] = MsaEntry(aSeq,name,strand,isReference=isReference)
        self._msa_data['sequence_names'].append( name )
        self._msa_data['alignment_length'] = max( self._msa_data['alignment_length'], aSeq.alignedRange.end )


    def setReferenceSeqName(self,n):
        "Sets Ref Seq Name"
        self._msa_data['reference_name'] = n

    def getReferenceSeqName(self):
        "Returns Ref Seq Name"
        return self._msa_data['reference_name']

    def getReferenceMSAItem(self):
        return self._msa_data['sequences'][self.getReferenceSeqName()]

    def getMSAItemBySeqName(self,seq_name):
        return self._msa_data['sequences'][seq_name]

    def setReference(self,b):
        "Sets boolean for presence of reference sequence"
        self._msa_data['has_reference'] = b

    def hasReference(self):
        "Returns boolean for presence of reference sequence"
        return self._msa_data['has_reference']

    def hasNames(self):
        "Returns true if this MSA has names stored for each sequence"
        return self._msa_data['has_names']

    def __iter__(self):
        for entry in self._msa_data['sequences'].values():
            yield entry.sequence

    def nameIter(self):
        return iter(self._msa_data['sequence_names'])
    
    def nameSeqIter(self):
        "Returns iterator that provides name, sequence pairs"
        for entry in self.entryIter():
            yield entry.name, entry.sequence

    def entryIter(self):
        "Returns iterator that provides MsaEntries"
        return imap( lambda x: self._msa_data['sequences'][x], \
            self._msa_data['sequence_names'] )
    
    def __len__(self):
        return len(self._msa_data['sequences'])
    
    def getAlignmentLength(self):
        return self._msa_data['alignment_length']
    
    def getColumnCalls(self,i):
        calls = self.getColumnCallsExt(i,isolateRef=False)
        return calls[0]

    ###Zero indexed
    def getColumnCallsExt(self, i, isolateRef=False):
        "Returns sequence of calls at column i in the alignment."
        "If isolateRef is specified, a second array is returned with the ref bases"
        if isinstance(i,int):
            pass
        else:
            raise Exception,"You must specify an integer index"

        read_calls = []
        ref_calls = []

        for aSeq in iter(self):
            if aSeq.alignedRange.contains(i):

                if isolateRef:
                    if aSeq.isReference:
                        ref_calls.append( aSeq[ i - aSeq.alignedRange.start ] )
                    else:
                        read_calls.append( aSeq[ i - aSeq.alignedRange.start ] )
                    
                else:
                    read_calls.append( aSeq[ i - aSeq.alignedRange.start ] )


        return (read_calls,ref_calls)
    

    ##Warning.. this appears to be 0 indexed, but not inclusive of the last integer specified. 0-5 will return subsequences of 5, not 6 length...
    def getRangeCallsExt(self, range, strict=True,isolateRef=False):
        """
        Returns sequence of calls in range in the alignment. 
        If strict is True, then only reports reads that completely contain the range.
        If strict is False, then reports any reads that overlap the range at all.
        if isolateRef is True, the reference sequence information is saved to a second array, else it is all  saved to the first.
        """
        ref_calls = []
        read_calls = []

        for name, aSeq in self.nameSeqIter():
            if (strict and aSeq.alignedRange.containsRange(range)) \
                or  \
               (not strict and range.intersects(aSeq.alignedRange)):
                subseq = ""
                # print str(range)
                for i in range:
                    idx = i - aSeq.alignedRange.start  
                    if (-1 < idx and idx < len(aSeq)):
                        subseq += aSeq[ idx ] 
                    else:
                        subseq += "-"

                if isolateRef:
                    if aSeq.isReference:
                        ref_calls.append(subseq)
                    else:
                        read_calls.append(subseq)
                    
                else:
                    read_calls.append(subseq)
                
        return(read_calls,ref_calls)



    def getRangeCalls(self, range, strict=True):
        calls = self.getRangeCallsExt(range,strict)[0]
        return calls

    
    def getConsistentRanges(self, minLength=None, maxLength=None):
        """Returns list of ranges with consistent sequences, where
	   consistent means all sequences extend throughout the entire range, 
	   subject to the minLength and maxLength constraints.
        
            i.e. the following MSA would return four ranges
                starting and ending at each |
        
            |-----------|----------|----------|
                        |----------|----------|-------------|
                                   |----------|
        """
        bounds =  [ i.alignedRange.start   for i in self ]
        bounds += [ i.alignedRange.end     for i in self ]
        bounds.sort()
        
        if (maxLength!=None):
            for boundIdx in range(len(bounds)-1):
                if (bounds[boundIdx+1] - bounds[boundIdx] > maxLength):
                    bounds += range(bounds[boundIdx] + maxLength, bounds[boundIdx+1], maxLength)
            bounds.sort()
        
        if (minLength!=None):
            newBounds = [ bounds[0] ]
            for bound in bounds[1:] :
                if (bound - newBounds[-1] < minLength): continue
                newBounds.append(bound)
            bounds = newBounds
            bounds.sort()
             
        ranges = []
        for boundIdx in range(len(bounds)-1):
            ranges.append( Range( bounds[boundIdx], bounds[boundIdx+1]) )
            # print str(ranges[-1])
        return ranges
    

 
    def clone(self):
        c = MSA()
        #c._sequences = [ MsaEntry(a.clone()) for a in iter(self) ]
        #c._alignmentLength = self._alignmentLength
        
        for orig_seq in self.entryIter():
            aSeq = AlignedSequence(sequence=orig_seq.sequence,seqName=orig_seq.name, \
                isReference=orig_seq.isReference() )

            c.add(aSeq,name=orig_seq.name,strand=orig_seq.strand,\
                isReference=orig_seq.isReference() )

        c._msa_data['alignment_length'] = self.getAlignmentLength()
        c._msa_data['has_names'] = self.hasNames()

        return c
       

class MsaEntry:
    def __init__( self, sequence=None, name=None, strand='+', isReference=False ):
        self.sequence = sequence
        self.name = name
        self.strand = strand
        self.sequence.isReference = isReference  # TODO change to method call?

    def isReference(self):
        return self.sequence.isReference
