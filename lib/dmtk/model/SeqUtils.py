__doc__="""Useful methods for sequence manipulation (and avoiding BioPython dependencies)."""
import sys
import string

COMP_TABLE = string.maketrans( 'ACGTRYKMWSHVBDNXacgtrykmwshvbdnx-', \
                               'TGCAYRMKWSDBVHNXtgcayrmkwsdbvhnx-' )

def revcomp( sequence ):
    """
    Return the reverse-complemented string for an 8-bit DNA string
    This method uses string.translate and won't work with
    unicode strings.
    """
    comp = sequence.translate( COMP_TABLE )
    # this is a fast string reverse from the python cookbook
    return comp[::-1]

def indexToKmer( index, k ):
    """Converts an index (integer) into a k-mer string.
    Inverse function of kmerToIndex."""
    kmer = [ 'A' for i in xrange(k) ]
    x = index
    j = k-1
    while x>0:
        kmer[ j ] = BASE_STRING[x & 3]
        x >>= 2
        j -= 1
    return "".join(kmer)
    
def kmerToIndex( kmer ):
    """Converts a k-mer into a unique integer (against all other k-mers
    for a fixed k)
    Useful for constructing efficient lookup tables."""
    index = base2int[ ord(kmer[0]) ]
    for ch in kmer[1:]:
        iBase = base2int[ ord(ch) ]
        index <<= 2
        index += iBase
    return index
    
base2int = []
BASE_STRING = 'ACGT'
 
def initializeLookupTables():
    BASE_SEQ = 'ACGTacgt'
    INDEX_SEQ = [ 0, 1, 2, 3, 0, 1, 2, 3 ]
    global base2int, bitCount
    base2int = [ -1 for i in xrange(256) ]
    for i in xrange(256):
        ch = chr(i)
        ind = BASE_SEQ.find(ch)
        if ind>=0:
            base2int[ i ] = INDEX_SEQ[ ind ]

initializeLookupTables()
