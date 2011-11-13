from h5py import *
from numpy import *
import bisect

class CmpH5Format:
    """
    Provides a lightweight way to manage the differences between the
    different CmpH5 files. This class is like a name service for cmpH5
    files.
    """
    def __init__(self, cmpH5):
        if ('Version' in cmpH5.attrs):
            self.VERSION = cmpH5.attrs['Version']
        
        self.ALN_INFO             = 'AlnInfo'
        self.REF_INFO             = 'RefInfo'
        self.MOVIE_INFO           = 'MovieInfo'
        self.REF_GROUP            = 'RefGroup'
        self.ALN_GROUP            = 'AlnGroup'
        self.ALN_INDEX_NAME       = 'AlnIndex'
        self.ALN_INDEX            = '/'.join([self.ALN_INFO, self.ALN_INDEX_NAME])

        self.REF_GROUP_ID         = '/'.join([self.REF_GROUP, 'ID'])
        self.REF_GROUP_PATH       = '/'.join([self.REF_GROUP, 'Path'])
        self.REF_OFFSET_TABLE     = '/'.join([self.REF_GROUP, 'OffsetTable'])

        self.ALN_GROUP_ID         = '/'.join([self.ALN_GROUP, 'ID'])
        self.ALN_GROUP_PATH       = '/'.join([self.ALN_GROUP, 'Path'])

        (self.ID, self.ALN_ID, self.MOVIE_ID, self.REF_ID, self.TARGET_START, self.TARGET_END, self.RC_REF_STRAND, 
         self.HOLE_NUMBER, self.SET_NUMBER, self.STROBE_NUMBER, self.MOLECULE_ID, self.READ_START, self.READ_END, self.MAP_QV, 
         self.N_MATCHES, self.N_MISMATCHES, self.N_INSERTIONS, self.N_DELETIONS, self.OFFSET_BEGIN, self.OFFSET_END, 
         self.N_BACK, self.N_OVERLAP) = range(0, 22)

        self.extraTables = [ '/'.join([self.ALN_INFO, x]) for x in cmpH5[self.ALN_INFO].keys() if not x == self.ALN_INDEX_NAME]


def rightmostBinSearch(vec, val):
    """
    Return the rightmost position in the vector vec of val. If val is
    absent then we return the leftmost position of the value:
    min(vec[vec > val]). If val is greater than all elements in vec we
    return len(vec).
    """
    assert(len(vec) > 0)

    i = bisect.bisect_left(vec, val)

    if (len(vec) == i):
        return(i)
    
    while (i + 1 < len(vec) and vec[i + 1] == val):
        i += 1

    return(i)
    

             
def leftmostBinSearch(vec, val):
    """
    Return the leftmost position in the vector vec of val. If val is
    absent then we return the lefternmost position for the value:
    max(vec[vec < val]). The time complexity here is potentially worse
    than log(n) because of the extra step of walking backwards.
    """
    assert(len(vec) > 0)
    i = bisect.bisect_left(vec, val)

    if (i == 0):
        return(i)
    elif (i == len(vec)):
        v = vec[i-1]
        i -= 1
    else:
        v = vec[i]

    if (v > val):
        i -= 1

    while (i > 0 and vec[i-1] == vec[i]):
        i -= 1

    return(i)


def numberWithinRange(s, e, vec):
    """
    Compute the number of elements in vec (where vec is sorted), such
    that each element, e obeys the following constraint: start <= e <
    end.
    """
    lI = leftmostBinSearch(vec, s)
    rI = rightmostBinSearch(vec, e)
    return(len(filter(lambda x : s <= x < e, vec[lI:rI])))
    

def computeIndices(tStart, tEnd):
    """
    Given a sorted (tStart, tEnd) compute a two-column matrix with
    columns nBack and nOverlap.

    nBack is defined as follows: given a position j, nBack is the
    offset for the smallest i, such that tEnd[i] > tStart[j], i.e., j
    - nBack == i.

    nOverlap is defined as follows: given a position j, nOverlap is |
    tEnd[i] > tStart[j] | for all i in 1,...,j - 1.
    """
    res = zeros(2 * len(tStart), dtype = "int32").reshape(len(tStart), 2)

    for i in range(len(tStart) - 1, 0, -1):
        j = i - 1
        nBack = 0
        nOver = 0

        while (j >= 0):
            if (tEnd[j] > tStart[i]):
                nOver += 1
                nBack = (i - j)
            j -= 1

        res[i, 0] = nBack
        res[i, 1] = nOver
    return(res)


def computeIndicesDP(tStart, tEnd):
    """
    Given a sorted tStart, tEnd compute a two-column matrix with
    columns nBack and nOverlap.

    nBack is defined as follows: given a position j, nBack is the
    offset for the smallest i, such that tEnd[i] > tStart[j], i.e., j
    - nBack == i.

    nOverlap is defined as follows: given a position j, nOverlap is |
    tEnd[i] > tStart[j] | for all i in 1,...,j - 1.

    This is a dynamic programming implementation and is *substantially*
    faster than computeIndices.
    """
    res = zeros(2 * len(tStart), dtype = "int32").reshape(len(tStart), 2)
    sortedEnds = sort(tEnd)

    for i in range(1, len(tStart)):
        nEnding = numberWithinRange(tStart[i - 1], tStart[i], sortedEnds - 1)

        if (nEnding == 0):
            res[i, 0] = res[i - 1, 0] + 1            
            res[i, 1] = res[i - 1, 1] + 1
        else:
            res[i, 1] = res[i - 1, 1] - nEnding + 1

            advance = 0
            for k in range(i - 1 - res[i - 1, 0], i):
                if (tEnd[k] > tStart[i]):
                    break
                advance += 1
            res[i, 0] = res[i - 1, 0] - advance + 1

    return(res)


def computeRefIndexTable(refIDs, refIDVector):
    """
    Compute a table of offsets for refID given the unique refIDs in
    the file. This implementation expects refIDs which go from [1,
    Inf).
    """
    def table(vector, bins):
        counts = dict(zip(bins, [0]*len(bins)))
        for i in xrange(0, len(vector)):
            counts[vector[i]]+= 1
        return counts

    bc = table(refIDVector, refIDs).values()
    sm = cumsum(bc)
    offsetStart = concatenate((array([0]), sm[:-1]))

    return array(zip(refIDs, offsetStart, sm), dtype = "uint32")

def getOverlappingRanges(tStart, tEnd, nBack, nOverlap, rangeStart, rangeEnd):
    """
    Return indices overlapping the range defined by [rangeStart,
    rangeEnd]. Here tStart, tEnd, nBack, nOverlap are vectors of
    length n sorted according to tStart and tEnd. The vectors nBack
    and nOverlap are typically produced by computeIndices[DP].
    """
    assert(rangeEnd >= rangeStart and 
           len(tStart) == len(tEnd) == len(nBack) == len(nOverlap))

    lM = leftmostBinSearch(tStart, rangeStart) 
    lM = lM - nBack[lM]
    rM = rightmostBinSearch(tStart, rangeEnd + .5)

    assert(rM >= lM and rM >= 0 and lM >= 0)

    if (lM == rM):
        return(array([], dtype = "uint32"))
    else:
        idxs   = array(range(lM, rM), dtype = "uint32")
        toKeep = array([True]*len(idxs))

        for i in range(0, len(idxs)):
            if (tEnd[idxs[i]] <= rangeStart):
                toKeep[i] = False
            else:
                continue
        return(idxs[toKeep])


def projectIntoRange(tStart, tEnd, rangeStart, rangeEnd):
    """
    Here, I project reads into the range defined by [rangeStart,
    rangeEnd]. Coverage can be most efficiently calculated by first
    obtaining all reads overlapping the range using the
    getOverlappingRanges function then projecting them into the same
    or smaller range.
    """
    assert(len(tStart) == len(tEnd))

    res = array([0] * (rangeEnd - rangeStart + 1))

    for i in range(0, len(tStart)):
        s = max(rangeStart, tStart[i]) - rangeStart
        e = min(rangeEnd,   tEnd[i] - 1) - rangeStart

        if (e >= s):
            res[s:(1 + e)] = res[s:(1 + e)] + 1
    return(res)





#############################################################################
##
## cmp.h5 wrappers.
##
#############################################################################
def getReadsInRange(cmpH5, coords, justIndices = False):
    """
    Return an ndarray representing the portion of the reads which
    overlap the range specfied by coords, where coords is a
    three-tuple composed of (refSeqID, rangeStart, rangeEnd).  Here,
    cmpH5 is an hdf5 object representing a pointer to a sorted cmp.h5
    file. The sorting of a file can be obtained using the command line
    tool sortCmpH5.py. Additionally, the format argument is used to
    specify whether the file is Astro or Springfield; the default is
    Springfield.
    """
    format = CmpH5Format(cmpH5)

    alignmentIndex = cmpH5[format.ALN_INDEX]
    offsets = cmpH5[format.REF_OFFSET_TABLE].value 
    refSeq, rangeStart, rangeEnd = coords
        
    offStart, offEnd = offsets[offsets[:,0] == refSeq, 1:3].ravel()
       
    if (offEnd - offStart > 0):
        refAlignIdx = alignmentIndex[offStart:offEnd, ]
        idxs = getOverlappingRanges(refAlignIdx[:,format.TARGET_START], refAlignIdx[:,format.TARGET_END],
                                    refAlignIdx[:,format.N_BACK], refAlignIdx[:,format.N_OVERLAP], 
                                    rangeStart, rangeEnd)
    else:
        ## This looks strange, but the idea is that a rowless matrix
        ## still has columns and these are what I want to preserve --
        ## h5py objects cannot be subset by a vector of length 0,
        ## however, numpy allows this.
        refAlignIdx = alignmentIndex[1:2, ]
        idxs = array([], dtype = 'uint32')

    if justIndices:
        return(idxs + offStart)
    else:
        return(refAlignIdx[idxs,])

def getCoverageInRange(cmpH5, coords):
    """
    Return a vector of length: coords[2] - coords[1] + 1 where each
    element represents the number of reads overlapping that position
    in the cmp.h5 file.
    """
    format = CmpH5Format(cmpH5)
    reads  = getReadsInRange(cmpH5, coords)

    if (reads.shape[0] == 0):
        return array([0]*(coords[2] - coords[1] + 1))
    else:
        return(projectIntoRange(reads[:,format.TARGET_START], reads[:,format.TARGET_END],
                                coords[1], coords[2]))
    

################################################################################
## 
## Testing/Timing code
## 
################################################################################
# import numpy.random
# import time
# import unittest
# import os
# import pdb

# class testSortedCmpH5Tools(unittest.TestCase):
#     def __brute_force_lm_search(self, vec, val):
#         if (val not in vec):
#             nvec = vec[ vec < val ]
#             if (len(nvec) == 0):
#                 return(0)
#             val = max(nvec)
#         for i in range(0, len(vec)):
#             if (vec[i] == val):
#                 break
#         return(i) 

#     def __brute_force_rm_search(self, vec, val):
#         if (val not in vec):
#             nvec = vec[ vec > val ]
#             if (len(nvec) == 0):
#                 return(len(vec))
#             val = min(nvec)
#             return(bisect.bisect_left(vec, val))
#         else:
#             return(bisect.bisect_right(vec, val) - 1)
    
#     def __brute_force_number_in_range(self, s, e, vec):
#         return(len(filter(lambda x : s <= x < e, vec)))

#     def __generate_positions(self, size, coverage, lScale = 50):
#         NN = size*coverage
#         tS = random.randint(0, size, NN)
#         tE = tS + array(map(int, random.exponential(lScale, NN) + 1))
#         ar = array([tS, tE]).transpose()
#         ar = ar[lexsort((tE, tS)),]
#         return(ar)
    
#     def __compare_implementations(self, size, coverage = 1):
#         NN = size * coverage
#         ar = self.__generate_positions(size, coverage)
#         res = computeIndices(ar[:,0], ar[:,1])
#         resDP = computeIndicesDP(ar[:,0], ar[:,1])
#         self.assertTrue(sum(res[:,0:2] == resDP[:,0:2]) == NN*2)

#     def __brute_force_search(self, tStart, tEnd, nBack, nOverlap, start, end):
#         toKeep = array([False]*len(tStart))
#         res = array(range(0, len(tStart)))

#         for i in range(0, len(tStart)):
#             # four cases to deal with.
#             if (tStart[i] >= start and tStart[i] <= end):
#                 toKeep[i] = True
#             elif (tEnd[i] > start and tEnd[i] <= end):
#                 toKeep[i] = True
#             elif (tStart[i] <= start and tEnd[i] >= end):
#                 toKeep[i] = True
#             elif (tStart[i] >= start and tEnd[i] <= end):
#                 toKeep[i] = True
#             else:
#                 continue
#         return(res[toKeep])

#     def __brute_force_get_reads(self, aIdx, start, end, format):
#         if aIdx.shape[0] == 0:
#             return aIdx

#         idxs = self.__brute_force_search(aIdx[:,format.TARGET_START],
#                                          aIdx[:,format.TARGET_END],
#                                          aIdx[:,format.N_BACK],
#                                          aIdx[:,format.N_OVERLAP],
#                                          start, end)
#         return(aIdx[idxs,])
        

#     def test_getOverlappingIndices(self):
#         for i in [100, 500]:
#             for j in [.1, 1, 5, 10]:
#                 for k in range(0, 10):
#                     ar = self.__generate_positions(i, j)
#                     idx = computeIndicesDP(ar[:,0], ar[:,1])
#                     aArray = hstack((ar, idx))
#                     s = random.randint(0, i, 1)
#                     e = int(1 + random.exponential(30, 1))
#                     x = getOverlappingRanges(aArray[:,0], aArray[:,1], aArray[:,2], aArray[:,3], s, s + e)
#                     y = self.__brute_force_search(aArray[:,0], aArray[:,1], aArray[:,2], aArray[:,3], s, s + e)
#                     self.assertTrue(all(sort(x) == sort(y)))


#     def test_constructIndices(self):
#         for i in [100, 200]:
#             for j in [1, 5]:
#                 self.__compare_implementations(i, j)

        
#     def test_leftmostBinSearch(self):
#         for j in range(0, 100):
#             a = sort(random.randint(0, 100, 100))
#             v = random.randint(0, 100, 1)
#             self.assertEqual(leftmostBinSearch(a, v), 
#                              self.__brute_force_lm_search(a, v))
        
#     def test_rightmostBinSearch(self):
#         for j in range(0, 100):
#             a = sort(random.randint(0, 100, 100))
#             v = random.randint(0, 100, 1)
#             self.assertEqual(rightmostBinSearch(a, v), 
#                              self.__brute_force_rm_search(a, v))

#     def test_numberWithinRange(self):
#         for j in range(0, 100):
#             a = sort(random.randint(0, 100, 100))
#             s,e = sort(random.randint(0, 100, 2))
#             self.assertEqual(numberWithinRange(s,e,a), 
#                              self.__brute_force_number_in_range(s,e,a))

#     def test_projectIntoRange(self):
#         tStart = array([1,1,1,1,1,2,2,2,2,10,20])
#         tEnd   = array([2,3,4,5,6,3,4,5,6,15,25])
#         self.assertTrue(all(projectIntoRange(tStart, tEnd, 1, 5) == array([5, 8, 6, 4, 2])))
#         self.assertTrue(all(projectIntoRange(tStart, tEnd, 20, 25) == array([1, 1, 1, 1, 1, 0])))

#     def test_computeRefIndexes(self):
#         refIDs = [1,1,1,1,1,2,2,2,2,2,1,1,1,1,5]
#         tbl    = computeRefIndexTable(array([1,2,3,4,5]), array(refIDs))
#         utbl   = [ 1,  0,  9,  2,  9, 14, 3, 14, 14, 4, 14, 14, 5, 14, 15]
#         self.assertEqual(sum(tbl.ravel() == utbl), len(utbl))
        
#     def test_getReadsInRange(self):
#         ## here is something not-portable, what is the right idiom?
#         files = [("~/projects/software/assembly/system-test/mock-data/outputs/rpal_three_ref_1x_coverage.sorted.cmp.h5", "Springfield"),
#                  ("~/aligned_reads.cmp.h5", "Springfield")]
      
#         for f in files:
#             cmpH5   = File(os.path.expanduser(f[0]), 'r')
#             format  = CmpH5Format(cmpH5)
#             refSeqs = cmpH5[format.REF_OFFSET_TABLE].value[:,0]

#             for nothing in range(0, 100):
#                 refSeq  = refSeqs[random.randint(0, len(refSeqs), 1)]
#                 start   = random.randint(0, 1e5, 1)
#                 stop    = int(1 + random.exponential(50, 1)) + start
#                 reads   = getReadsInRange(cmpH5, (refSeq, start, stop), format = format)

#                 alnIdx  = cmpH5['AlignmentIndex'].value
#                 readsBF = self.__brute_force_get_reads(alnIdx[alnIdx[:,format.REF_ID] == refSeq,], 
#                                                        start, stop, format = format)

#                 self.assertTrue(all(reads.shape == readsBF.shape))
#                 self.assertTrue(all(sort(reads[:,format.TARGET_START]) == sort(readsBF[:,format.TARGET_START])))
         

# if __name__ == "__main__":
#     suite1 = unittest.TestLoader().loadTestsFromTestCase(testSortedCmpH5Tools)
#     unittest.TextTestRunner(verbosity=2).run(suite1)

