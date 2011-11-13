"""
Sorted cmp.h5 column iterator.

Author: Andrey Kislyuk
"""

import sys, logging
import numpy, scipy
from dmtk.io import cmph5
from dmtk.io.cmph5 import CmpH5SortingTools

class ColumnIterator:
    """ Supports iteration on alignment columns in reference coordinates.
    To identify the reference to iterate over, requires ReferenceAnnotation cmph5.refInfoIterator()
    
    The default iterator returns (ref_pos, [(AlignmentHit, i), ...]) tuples where i indexes into the columns of the corresponding AlignmentHit.
    """
    def __init__(self, cmpH5, ref_info, ref_start = 0, ref_end = None, pulse_tables=['PulseWidth', 'IPD'],
                 target_coverage_ceiling = None):
        if not cmpH5.isSorted():
            raise StandardError("Iteration not supported on unsorted cmp.h5")
        if cmpH5.version < "1.2.0":
            raise StandardError("Iteration not supported on pre-1.2.0 cmp.h5, use cmpH5Convert.py")
        
        self.ref_start, self.ref_end = ref_start, ref_end
        self.aln_hit_cache = {}
        
        self.cmp, self.ref_info, self.pulse_tables = cmpH5, ref_info, pulse_tables
        try:
            self.ref_group = self.cmp.refGroupFromFullName(ref_info.fullName)
        except KeyError:
            raise StopIteration()
        
        if len(self.ref_group) < 1: raise StopIteration()
        
        self.ref_start = int(self.ref_start)
        if self.ref_end is None or self.ref_end > int(self.ref_info.length): self.ref_end = int(self.ref_info.length)
        
        self.aln_group_aln_index_rows_cache = {}
        self.num_aln_hits = 0
        for aln_group in self.ref_group.alnGroupIterator():
            self.aln_group_aln_index_rows_cache[aln_group.id] = set(aln_group._myAlignmentIndexRows())
            self.num_aln_hits += aln_group.numAlnHits

        # TODO: dynamic stride adjustment (double/halve stride if nreads under/over target)
        stride_base = 5000
        self.stride = stride_base*self.ref_info.length/self.num_aln_hits
        if self.stride < 100: self.stride = 100
        if self.stride > 100000: self.stride = 100000
        if self.stride < 1000:
            self.stride = int(self.stride/100)*100
        else:
            self.stride = int(self.stride/1000)*1000
        
        if self.stride > self.ref_end - self.ref_start + 1: self.stride = self.ref_end - self.ref_start + 1
        
        self.target_coverage_ceiling = target_coverage_ceiling
        if target_coverage_ceiling:
            subreadlengths = []
            for aln_group in self.ref_group.alnGroupIterator():
                my_subreadlengths = [cmpH5["/AlnInfo"].asRecArray()[row]["tEnd"] - cmpH5["/AlnInfo"].asRecArray()[row]["tStart"] for row in self.aln_group_aln_index_rows_cache[aln_group.id]]
                subreadlengths.extend(my_subreadlengths)
            self.expected_subread_length = int(numpy.mean(my_subreadlengths))
            self.num_coverage_trim_events = 0
        
        logging.debug("%s started on %s:%s..%s in %s, stride %d" % (str(self), self.ref_group.name, self.ref_start, self.ref_end, self.cmp.h5File.filename, self.stride))
    
    def __iter__(self):
        num_iters = 0
        
        # TODO: dynamic stride adjustment (double/halve stride if nreads under/over target)
        for coords_start in range(self.ref_start, self.ref_end+1, self.stride):
            cur_coords = (self.ref_group.id, coords_start, coords_start+self.stride)
            reads_in_range = CmpH5SortingTools.getReadsInRange(self.cmp.h5File, cur_coords, justIndices = True)
            if len(reads_in_range) < 1: continue
            
            if self.target_coverage_ceiling:
                max_nreads = self.target_coverage_ceiling * self.stride / self.expected_subread_length
                if len(reads_in_range) > max_nreads:
                    reads_in_range = reads_in_range[::len(reads_in_range)/max_nreads][:max_nreads]
                    self.num_coverage_trim_events += 1
            
            reads_in_range = set(reads_in_range)
            
            # logging.debug("Processing "+self.ref_group.name+": "+str(coords_start)+", reads in range: "+str(len(reads_in_range)))
            
            for r in self.aln_hit_cache.keys():
                if r not in reads_in_range:
                    del self.aln_hit_cache[r]
            
            for aln_group in self.ref_group.alnGroupIterator():
                aln_group_rir = reads_in_range & self.aln_group_aln_index_rows_cache[aln_group.id]
                self.aln_hit_cache.update([(r, aln_group._alnHitFromRowNum(r, auxTables=self.pulse_tables, computeAlnPosTranslations=True)) for r in aln_group_rir if r not in self.aln_hit_cache])

            rstarts, rends = {}, {}
            cur_rows = set()
            for r, aln_hit in self.aln_hit_cache.iteritems():
                if aln_hit.target_start >= coords_start and aln_hit.target_start < coords_start+self.stride:
                    rstarts.setdefault(aln_hit.target_start, [])
                    rstarts[aln_hit.target_start].append(r)
                
                if aln_hit.target_end >= coords_start and aln_hit.target_end < coords_start+self.stride:
                    rends.setdefault(aln_hit.target_end, [])
                    rends[aln_hit.target_end].append(r)
                
                if aln_hit.target_start <= coords_start and aln_hit.target_end > coords_start:
                    cur_rows.add(r)            

            for ref_pos in range(coords_start, coords_start+self.stride):
                if ref_pos > self.ref_end: raise StopIteration()
                #logging.debug("%s: pos %d, cache size %d, range %d .. %d" % (str(self), ref_pos, len(self.aln_hit_cache),  coords_start, coords_start+self.stride))

                if ref_pos in rends:
                    cur_rows.difference_update(rends[ref_pos])
                if ref_pos in rstarts:
                    cur_rows.update(rstarts[ref_pos])
                
                if len(cur_rows) < 1: continue
                
                #aln_contexts = [(self.aln_hit_cache[r], self.aln_hit_cache[r].targetToAlnPos[ref_pos-self.aln_hit_cache[r].target_start]) for r in cur_rows]
                aln_contexts = [(self.aln_hit_cache[r], self.aln_hit_cache[r]._target2aln[ref_pos-self.aln_hit_cache[r].target_start]) for r in cur_rows]
                
                num_iters += 1
                yield ref_pos, aln_contexts
                
        #logging.debug("%s completed %d iterations on %s:%s..%s in %s" % (str(self), num_iters, self.ref_group.name, self.ref_start, self.ref_end, self.cmp.h5File.filename))
        if self.target_coverage_ceiling and self.num_coverage_trim_events > 0:
            logging.info("%s: coverage trimmed in %d windows" % (str(self), self.num_coverage_trim_events))


class RolledColumnIterator:
    ''' Supports iteration on unrolled circular templates, collating contexts per template position
    '''
    def __init__(self, cmpH5, ref_info, template_length, pulse_tables=['PulseWidth', 'IPD']):
        self.cmpH5, self.ref_info, self.template_length, self.pulse_tables = cmpH5, ref_info, int(template_length), pulse_tables

    def __iter__(self):
        col_iterators = {}
        iter_buffer = {}
        for tile_start in range(0, int(self.ref_info.length), self.template_length):
            col_iterators[tile_start] = ColumnIterator(self.cmpH5, self.ref_info, tile_start, tile_start+self.template_length, self.pulse_tables).__iter__()
            iter_buffer[tile_start] = None
        
        logging.debug("%s started on %s in %s, template length %d, %d child iterators" %
                      (str(self), self.ref_info.fullName, self.cmpH5.h5File.filename, self.template_length, len(col_iterators)))

        do_ref_col_sanity_check = True
        for template_pos in range(0, self.template_length):
            # logging.debug("%s: pos %d" % (str(self), template_pos))
            contexts_this_pos = []
            # replenish cache and remove iterators which are done
            for tile_start, iterator in col_iterators.iteritems():
                if iterator is None: continue
                if iter_buffer[tile_start] is not None: continue
                try:
                    iter_buffer[tile_start] = iterator.next()
                except StopIteration:
                    del iter_buffer[tile_start]
                    col_iterators[tile_start] = None
            
            # pop items from cache only if reference positions match
            for tile_start, i in iter_buffer.iteritems():
                ref_pos, aln_contexts = i
                if ref_pos % self.template_length == template_pos:
                    contexts_this_pos.extend(aln_contexts)
                    iter_buffer[tile_start] = None

            if do_ref_col_sanity_check and len(contexts_this_pos) > 0: # sanity check: assert same template nucleotide in all contexts
                col_nt = [aln_hit.alignedTarget[aln_pos] for aln_hit, aln_pos in contexts_this_pos]
                if len(set(col_nt)) != 1:
                    raise StandardError('Target nucleotides do not match for template position %d, check template length %d' % (template_pos, self.template_length))
                do_ref_col_sanity_check = False
            
            yield template_pos, contexts_this_pos


class CmpAlnCol(object): 
    """Represents a column of basecalls and associated data in the MSA.
    Unused, from deprecated pre-1.2 code."""
    def __init__( self, parent, reads=None, offset=0 ):
        self._parent = parent
        # list of CmpBaseEvidence
        self._bases = []
        self.refPos = 0
        self.offset = offset

    def append( self, evidence ):
        "Add a basecall (CmpBaseEvidence) to this column"
        self._bases.append( evidence )

    def __iter__( self ):
        "Iterator over CmpBaseEvidence"
        return iter(self._bases)

    def __len__( self ):
        "Number of bases in this column" 
        return len(self._bases)

    def __str__( self ):
        "for debugging"
        bases = ''.join( [ e.basecall for e in self._bases ] )
        strands = ''.join( [ e.strand for e in self._bases ] )
        return 'CAC { refPos=%d; offset=%d; bases=%s; strands=%s }' % \
            ( self.refPos, self.offset, bases, strands )

class CmpBaseEvidence(object):
    """Struct containing the evidence for an individual
    base and associated information in an alignment column.
   
    Subclasses can extend the information content.
    Unused, from deprecated pre-1.2 code.
    """
    def __init__( self, basecall, strand, read_id='' ):
        self.basecall = basecall
        # '+', '-'
        self.strand = strand
        self.read_id = read_id


# Deprecated, from pre-1.2 code.

#from collections import deque
#import itertools as it
#
#def sortedColumnIterator(cmpH5, refSeqName, refStart = 0, refEnd = None, __blockSize = 200, __noWarn = False):
#
#    format = CmpH5Format(cmpH5)
#    if (not __noWarn):
#        print "This is not production grade code yet! If you insist, and don't want to see this message, then call the function with __noWarn = True"
#
#    SMAP = {1:"-", 0:"+"}
#
#    refGroup = cmpH5[refSeqName]
#    refGroupID = refGroup.id
#
#    offsets = cmpH5['SortedRefIndexTable'].asNumPy
#    offStart, offEnd = offsets[offsets[:,0] == refGroupID, 1:3].ravel()
#    refAlignIdx = cmpH5['/AlnInfo/AlnIndex']._dataset[offStart:(offEnd + 1), ]
#    maxRefPos = max(refAlignIdx[:, format.TARGET_END])
#
#    def myPop(q):
#        if (not len(q)):
#            return(None)
#        else:
#            return(q.pop())
#    
#    def myPush(q, v):
#        q.appendleft(v)
#    
#    def getReadGroup(rgID):
#        rgPath = cmpH5["/AlnGroup/Path"][rgID]
#        return cmpH5[ rgPath ]
#
#    def _sortedColumnIterator(refStart, refEnd): 
#        aArrayCache = {}
#
#        def getAlignmentPair(idx):
#            rgID = reads[idx, format.ALN_ID]
#            if (not aArrayCache.has_key(rgID)):
#                aArrayCache[reads[idx, format.ALN_ID]] = getReadGroup(rgID)["AlnArray"]._dataset.value
#            aA = aArrayCache[rgID]
#            return(aA[currentPos[idx] + reads[idx, format.OFFSET_BEGIN]])
#
#        def isGap(idx):
#            return(getAlignmentPair(idx) in array([16, 32, 64, 128, 240]))
#
#        def readBase(idx):
#            return(alignmentPairMap[getAlignmentPair(idx)][0])
#
#        ## obtain the reads within the region.
#        idxs = getOverlappingRanges(refAlignIdx[:,format.TARGET_START], refAlignIdx[:,format.TARGET_END],
#                                    refAlignIdx[:,format.N_BACK], refAlignIdx[:,format.N_OVERLAP], 
#                                    refStart, refEnd)
#        if (not len(idxs)):
#            raise StopIteration()
#
#        reads = refAlignIdx[idxs, ]
#
#        currentPos = map(lambda x : x if x > 0 else 0, refStart - reads[:,format.TARGET_START])
#        currentPos = currentPos + reads[:,format.RC_REF_STRAND]*(reads[:,format.OFFSET_END] - reads[:,format.OFFSET_BEGIN])
#        step = array([1, -1])[reads[:, format.RC_REF_STRAND]]
#        
#        begin = 0
#        for refPos in range(refStart, refEnd + 1):
#            inList = []
#        
#            ## what reads overlap this ref position, 
#            ## given that their sorted. There might
#            ## be a substantially faster way to do this. 
#            for i in range(begin, reads.shape[0]):
#                read = reads[i,]
#
#                if (read[format.TARGET_START] > refPos):
#                    break
#                elif (read[format.TARGET_END] < refPos):
#                    begin = i 
#                    continue
#                else:
#                    inList.append((i, read))
#
#            ## make a queue for each read. 
#            queues = [ deque() for i in range(0, len(inList)) ]
#
#            for i, (idx, read) in enumerate(inList):
#                myPush(queues[i], readBase(idx))
#                currentPos[idx] += step[idx]
#
#                while (currentPos[idx] <= (read[format.TARGET_END] - read[format.TARGET_START])):
#                    if (isGap(idx)):
#                        myPush(queues[i], readBase(idx))
#                        currentPos[idx] += step[idx]
#                    else:
#                        break
#
#            cOff = 0
#            while (True):
#                alnCol = CmpAlnCol(refGroup)
#                alnCol.refPos = refPos
#                alnCol.offset = cOff
#
#                for i in xrange(0, len(queues)): #, read in zip(queues, inList):
#                    x = myPop(queues[i])
#                    if (x):
#                        alnCol.append(CmpBaseEvidence(x, SMAP[inList[i][1][format.RC_REF_STRAND]], inList[i][1][format.ID]))
#                    else:
#                        alnCol.append(CmpBaseEvidence('-', SMAP[inList[i][1][format.RC_REF_STRAND]], inList[i][1][format.ID]))
#                yield alnCol
#
#                cOff += 1
#                if (not any(map(lambda x : len(x) > 0, queues))):
#                    break
#
#    ## 
#    ## This is the main driver loop for the iterator. 
#    ##
#    currentStart = refStart
#    currentEnd   = currentStart + __blockSize
#    while (True):
#        itr = _sortedColumnIterator(currentStart, currentEnd)
#        for col in itr:
#            yield col
#        
#        if (currentStart > maxRefPos): 
#            raise StopIteration()
#        else:
#            currentStart = currentEnd + 1
#            currentEnd = currentStart + __blockSize
#
#
#def unsortedColumnIterator( cmpH5, refSeqName, refPosition=0 ):
#    """Returns iterator over CmpAlnCol using slow logic which handles
#    the older unsorted cmp.h5."""
#    globalIndex = cmpH5["/AlnInfo/AlnIndex"]._dataset
#    tStartCol = cmpH5._colName2Index["tStart"]
#    tEndCol = cmpH5._colName2Index["tEnd"]
#    refSeqIdCol = cmpH5._colName2Index["RefGroupID"]
#    readGroupIdCol = cmpH5._colName2Index["AlnGroupID"]
#
#    sortI = numpy.argsort(globalIndex[...,tStartCol])
#    sortStartT = numpy.sort(globalIndex[...,tStartCol])
#    endT = (globalIndex[...,tEndCol][x] for x in sortI)
#
#    refSeqGroup = cmpH5[ refSeqName ]
#    refGroups = cmpH5["/RefGroup"].asRecArray( )
#    refInfos = cmpH5["/RefInfo"].asRecArray( )
#    refInfoID = refGroups[ refGroups["Path"] == refSeqName ][0]["RefInfoID"]
#    refLength = refInfos[ refInfos["ID"] == refInfoID ][0]["Length"]
#    
#    if refLength < 0:
#        raise Exception( "Column iteration requires reference length in ref group (RefGroup=%s)" % refSeqName )
#
#    theRefSeqId = cmpH5["/RefGroup"].getID( Path=refSeqName )
#    refSelector = ( theRefSeqId==globalIndex[...,refSeqIdCol][x] for x in sortI )
#    sortedIndex = [x[:3] for x in it.izip(sortI, sortStartT, endT, refSelector ) if x[3] ]
#
#    s1 = bisect.bisect_left(sortStartT, refPosition)
#    sortedIndex2 = sortedIndex[:s1]
#    sortedIndex2.sort(key=lambda x:x[2])                  
#    keys = [x[2] for x in sortedIndex2]
#    s2 = bisect.bisect_left(keys, refPosition+1)  
#    
#    rg_path_table = cmpH5["/AlnGroup/Path"] #StrDataset( cmpH5._group["ReadGroupPath"] )
#    def get_read_group( row ):
#        readGroupPath = rg_path_table[ int(row[ readGroupIdCol ]) ]
#        return cmpH5[ readGroupPath ]
#
#    def add_u2a_map( hit ):
#        # target_map is mapping of unaligned to aligned coords
#        hit.target_map = numpy.zeros( len(hit.alignedTarget), 'uint16' )
#        i, nGaps = 0, 0
#        for ch in hit.alignedTarget:
#            if ch=='-':
#                nGaps += 1
#            else:
#                hit.target_map[i] = nGaps + i
#                i += 1
#
#    def read_ids_to_hits( ids ):
#        for idx, row in it.imap( lambda x: ( x[0], globalIndex[x[0],...] ), ids ):
#            readGroup = get_read_group( row )
#            hit = readGroup._alnHitFromRowNum( idx )
#            add_u2a_map( hit )
#            yield hit
#    #
#    # reads which intersect the current reference position
#    #
#    readIds = sortedIndex2[s2:]  #all reads end after the currentRefPos
#    # all reads starts at the currentRefPos                
#    while s1 < len(sortedIndex) and sortedIndex[s1][1] <= refPosition:  
#        readIds.append(sortedIndex[s1])
#        s1+=1
#
#    currentReads = []                                     
#    for hit in read_ids_to_hits( readIds ):
#        currentReads.append( hit )
#
#    currentRefPos = refPosition
#    currentInsertOffset = 0
#    
#    while currentRefPos < refLength:
#        alignColumn = CmpAlnCol( refSeqGroup )
#        alignColumn.refPos = currentRefPos
#        alignColumn.offset = currentInsertOffset
#        moreInsertions = False
#
#        def bases_for_ref( offset ):
#            endOffset = len(hit.alignedTarget)
#            for i, ch in enumerate(hit.alignedTarget[(offset+1):]):
#                if ch!='-':
#                    endOffset = offset + i + 1
#                    break
#            return hit.alignedQuery[ offset:endOffset ]
#
#        for hit in currentReads:
#            unalignedOffset = currentRefPos - hit.target_start
#            offset = hit.target_map[ unalignedOffset ]
#            readBases = bases_for_ref( offset )
#            strand = hit.target_strand
#            if currentInsertOffset < len(readBases):
#                moreInsertions = True
#                bc = readBases[currentInsertOffset]
#                alignColumn.append( CmpBaseEvidence( bc, strand, hit.query_id ) )
#            else:
#                alignColumn.append( CmpBaseEvidence( '-', strand, hit.query_id ) )
#
#        if moreInsertions:
#            yield alignColumn
#            currentInsertOffset += 1
#        else:
#            if len(currentReads)==0:
#                yield alignColumn
#            currentRefPos += 1
#            currentInsertOffset = 0
#            currentReads = [ h for h in currentReads if h.target_end > currentRefPos ] #faster way to remove reads?
#
#            readIds = []  #all reads end after the currentRefPos
#            while s1 < len(sortedIndex) and sortedIndex[s1][1] <= currentRefPos:
#                readIds.append(sortedIndex[s1])
#                s1+=1
#
#            for hit in read_ids_to_hits( readIds ):
#                currentReads.append( hit )
