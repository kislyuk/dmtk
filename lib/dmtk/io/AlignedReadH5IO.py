import sys
import h5py
import numpy
import bisect

import unittest
from dmtk.model.AlignmentHit import AlignmentHit

__version__="$Revision: #40 $ $Change: 50848 $"

# This constant is/will be defined in the spec as a signed integer value for no data.
INT_NODATA = -99999

cMap = dict(zip('NACGTacgt-','NTGCAtgca-'))

base2hexMap = {"A":1, "C":2, "G":4, "T":8, "a":1, "c":2, "g":4, "t":8, "-":0, "N":15}
rBase2hexMap = {"T":1, "G":2, "C":4, "A":8, "t":1, "g":2, "c":4, "a":8, "-":0, "N":15}
alignmentPairMap = {}
rAlignmentPairMap = {}

for b1 in ['A','C','G','T','-', 'N']:
    for b2 in ['A','C','G','T','-', 'N']:
        alignmentPairMap[ (base2hexMap[b1]<< 4) | base2hexMap[b2] ] = (b1,b2) 

for b1 in ['A','C','G','T','-', 'N']:
    for b2 in ['A','C','G','T','-', 'N']:
        rAlignmentPairMap[ (rBase2hexMap[b1]<< 4) | rBase2hexMap[b2] ] = (b1,b2) 

Basemap = numpy.array(['-', 'A', 'C', '-', 'G', '-', '-', '-', 'T', 
                       '-', '-', '-', '-', '-', '-', '-', '-', 'N'])
RBasemap = numpy.array(['-', 'T', 'G', '-', 'C', '-', '-', '-', 'A',
                        '-', '-', '-', '-', '-', '-', '-', '-', 'N'])

def revCompSeq(s):
    return "".join([cMap[c] for c in s[::-1]])

class AlignedReadH5F(object):

    def __init__(self, filename, mode):
        self._h5FileHandle = h5py.File(filename, mode)
        self._assignAlignmentIndexCol()
        self._refSeqName2Id = {}
        self._readGroupPath2Id = {}
        self.mode = mode
        self._vlType = h5py.new_vlen(str)

        if mode in ['a','w']:
            self._createIndexTables()
        if mode != 'w':
            self._globalIndex = self._h5FileHandle["AlignmentIndex"]
            self._refSeqName = self._h5FileHandle["RefSeqName"]
            self._refSeqID = self._h5FileHandle["RefSeqID"]
            self._readGroupPath = self._h5FileHandle["ReadGroupPath"]
            self._readGroupPathID = self._h5FileHandle["ReadGroupPathID"]
    
        self._updateRefSeqDict()
        self._updateReadGroupDict()

    def _assignAlignmentIndexCol(self):
        self._indexColNames = ("AlignmentId", "ReadGroupId", 
                               "refSeqId", "tStart", "tEnd", "AlignedStrand", 
                               "uniqReadId", "subReadId", "rStart", "rEnd",
                               "nM", "nMM", "nIns", "nDel")
        self._indexColNum = len(self._indexColNames)
            
    def _createIndexTables(self):
        initialIndexSize = 100
        if "AlignmentIndex" not in self._h5FileHandle:
            self._globalIndex = self._h5FileHandle.create_dataset("AlignmentIndex", (initialIndexSize, self._indexColNum), dtype=h5py.h5t.NATIVE_INT32 , maxshape=(None,self._indexColNum))
            self._globalIndex.attrs.create("lastRow", data=0, dtype="u4")
            for i in range(self._indexColNum):
                self._globalIndex.attrs.create("ColName%02d" % i, data=self._indexColNames[i])
        
        initialRefSeqNameTableSize = 32
        if "RefSeqName" not in self._h5FileHandle:
            #use two simple column tables to avoid compound type (Well, I do like compound type but it might cause some problems for implementing it in .NET.), JC
            self._refSeqName = self._h5FileHandle.create_dataset("RefSeqName", (initialRefSeqNameTableSize,), self._vlType, maxshape=(None,)) #Reference name can not be more than 64 byte C string
            self._refSeqID = self._h5FileHandle.create_dataset("RefSeqID", (initialRefSeqNameTableSize,), "u4", maxshape=(None,)) #Reference seq. ID is 4 byte un-sign integer 
            self._refSeqID.attrs.create("lastRow", data=0, dtype="u4")
        
        initialReadGrouoPathSize = 128
        if "ReadGroupPath" not in self._h5FileHandle:
            #use two simple column tables to avoid compound type (Well, I do like compound type but it might cause some problems for implementing it in .NET.), JC
            self._readGroupPath = self._h5FileHandle.create_dataset("ReadGroupPath", (initialReadGrouoPathSize,), self._vlType, maxshape=(None,)) #ReadGroupPath can not be more than 64 byte C string
            self._readGroupPathID = self._h5FileHandle.create_dataset("ReadGroupPathID", (initialReadGrouoPathSize,), "u4", maxshape=(None,)) #ReadGroupPath ID is 4 byte un-sign integer 
            self._readGroupPathID.attrs.create("lastRow", data=0, dtype="u4")
        

    def _updateRefSeqDict(self):
        lastRow = self._refSeqID.attrs["lastRow"]
        if lastRow != 0:
            ids = self._refSeqID[:lastRow]
            names = self._refSeqName[:lastRow]    
            self._refSeqName2Id = dict(zip(names, ids))
        else:
            self._refSeqName2Id = {}


    def _updateReadGroupDict(self):
        lastRow = self._readGroupPathID.attrs["lastRow"]
        if lastRow != 0:
            ids = self._readGroupPathID[:lastRow]
            names = self._readGroupPath[:lastRow]
            self._readGroupPath2Id = dict(zip(names, ids))
        else:
            self._readGroupPath2Id = {}

    def setRootAttr(self, key, value):
        self._h5FileHandle.attrs.create(key, data=value)  #should make this type-safe

    def addNewReference(self, refSeqGroupName, annotationString = None):
        if refSeqGroupName not in self._h5FileHandle: # only create a new group if the reference group does not exist
            refSeqGroup = self._h5FileHandle.create_group(refSeqGroupName)
            if annotationString:
                refSeqGroup.attrs["annotationString"] = annotationString
            #update the Key-Value pair for the reference tables
            lastRow = self._refSeqID.attrs["lastRow"]
            
            #resize the table if necessary
            refSeqTableSize = self._refSeqID.shape[0]
            if lastRow + 1 > refSeqTableSize: 
                  self._refSeqID.resize( (lastRow+32, ) )
                  self._refSeqName.resize( (lastRow+32, ) )
            
            self._refSeqID[lastRow] = lastRow + 1
            self._refSeqName[lastRow] = refSeqGroupName
            
            self._refSeqID.attrs.modify("lastRow", lastRow + 1)
            self._refSeqName2Id[refSeqGroupName] = lastRow + 1
            
        #self._updateRefSeqDict()

    def addReadGroup(self, refSeqName, readGroupName):  #for Astro the readGroupName is the "expid-runid" as default
        refSeqGroup = self.getHDF5ObjectByPath(refSeqName)
        readGroup = refSeqGroup.create_group(readGroupName)

        alignmentArray = readGroup.create_dataset("AlignmentArray", (10*500,), dtype=h5py.h5t.NATIVE_UINT8, maxshape=(None,))
        alignmentArray.attrs.create("lastRow", data=0, dtype="u4")

        alignmentIndex = readGroup.create_dataset("AlignmentIndex", (100,self._indexColNum), dtype=h5py.h5t.NATIVE_INT32, maxshape=(None,self._indexColNum))
        alignmentIndex.attrs.create("lastRow", data=0, dtype="u4")

        groupPath = refSeqGroup.name + "/" + readGroupName
        #update the Key-Value pair for the reference tables
        lastRow = self._readGroupPathID.attrs["lastRow"]
            
        #resize the table if necessary
        tableSize = self._readGroupPathID.shape[0]
        if lastRow + 1 > tableSize: 
            self._readGroupPathID.resize( (lastRow+128, ) )
            self._readGroupPath.resize( (lastRow+128, ) )

        self._readGroupPathID[lastRow] = lastRow + 1
        self._readGroupPath[lastRow] = groupPath
        self._readGroupPath2Id[groupPath] = lastRow+1
        self._readGroupPathID.attrs.modify("lastRow", lastRow + 1)
        
        #self._updateReadGroupDict()
        
    def close(self):
        self._h5FileHandle.close()
        
    def getHDF5ObjectByPath(self, path):
        return self._h5FileHandle[path]
          
    def writeAnAlignment(self, alignment):
        raise "Not implemented Yet"
        pass
    
    def writeAlignments(self, efSeqName, readGroupName, alignmentDataList):
        raise "Not implemented Yet"
        pass

class AlignedAstroReadH5F(AlignedReadH5F):
    
    def __init__(self, filename, mode):
        AlignedReadH5F.__init__(self, filename, mode)
        if mode in ['a','w']:
            initialMovieTableSize = 1024
            if "MovieName" not in self._h5FileHandle:
                #use two simple column tables to avoid compound type (Well, I do like compound type but it might cause some problems for implementing it in .NET.), JC
                self._movieName = self._h5FileHandle.create_dataset("MovieName", (initialMovieTableSize,), self._vlType, maxshape=(None,)) #Movie name can not be more than 64 byte C string
                self._movieID = self._h5FileHandle.create_dataset("MovieID", (initialMovieTableSize,), "u4", maxshape=(None,)) #Movie seq. ID is 4 byte un-sign integer 
                self._movieID.attrs.create("lastRow", data=0, dtype="u4")
        if mode != 'w':
            self._movieName = self._h5FileHandle["MovieName"]
            self._movieID = self._h5FileHandle["MovieID"]
        
        self._updateMovieDict()

    def _parseReadId(self, id):

        # set defaults
        movieDate = "m000000"; movieTime = "000000"; 
        machine = "2nd"; block = "b00"; 
        segment = "0"; panel = "p0"; subreadId = 0;
        zmwX = "x0"; zmwY = "y0"; runName = "1000000-1000";

        # parse the various flavors of ids
        id = id.split("|")[0] # remove any annotations
        fields = id.split("_")
        if len(fields) == 8:
            zmwX, zmwY, runName, movieDate, movieTime, machine, \
                panel, block = fields
            subreadId = 0
            if len(block.split(".")) > 1:
                ( block, segment ) = block.split(".") # _block.subread
                subreadId = segment[1:]
        # ids that come directly from pls.h5 seem to be missing runcode
        elif len(fields) == 7:
            zmwX, zmwY, movieDate, movieTime, machine, \
                panel, block = fields
            subreadId = 0
            if len(block.split(".")) > 1:
                ( block, segment ) = block.split(".") # _block.subread
                subreadId = segment[1:]
        elif len(fields) == 9:
            zmwX, zmwY, runName, movieDate, movieTime, machine, \
                panel, block, segment = fields
            subreadId = int(segment[1:])
        elif len(fields) == 5:
            zmwX, zmwY, runName, panel, subreadId = fields
            subreadId = int(subreadId[1:])
        elif len(fields) == 3 and "-" in fields[2]:
            zmwX, zmwY, runName = fields
        elif len(fields) == 1:
            zmwX = "x%i" % int(fields[0]) # TODO is this really the best way to handle naive ids?
        else:
            # the line after this is broken.  At least print out some error state
            raise SystemExit, "Can't parse read id %s" % id

        return zmwX, zmwY, runName, movieDate, movieTime, machine, \
                    panel, block, segment, subreadId

    def constructAlignData(self, hit):

        # have to grab out the cr=X,Y if it exists so that
        # our read_start and read_end are accurate.
        # parses out |cr=12,32
        cr = (0, 0)
        for annotation in hit.query_id.split("|"):
            if annotation.split("=")[0] == "cr":
                cr = [ int(x) for x in annotation.split("=")[1].split(",") ]

        zmwX, zmwY, runName, movieDate, movieTime, machine, \
            panel, block, segment, subreadId = self._parseReadId(hit.query_id)

        al = abs( hit.query_end - hit.query_start )
        nCorrect = al - hit.nIns - hit.nDel - hit.nMismatch

        zmwX = zmwX[1:]
        zmwY = zmwY[1:]
        movieName = "_".join( (movieDate, movieTime, machine, panel, block) )
        
        panel = panel[1:]
        strand = hit.query_strand
        targetStrand = hit.target_strand

        alignedRead = hit.alignedQuery
        alignedTarget = hit.alignedTarget
        if strand == "-":
            alignedRead = revCompSeq(alignedRead)
            alignedTarget = revCompSeq(alignedTarget)
        if strand != targetStrand:
            relativeStrand = 1
        else:
            relativeStrand = 0

        tStart = hit.target_start
        tEnd = int(hit.target_end)
        if tStart > tEnd:
            tStart, tEnd = tEnd, tStart

        rStart = int(hit.query_start) + cr[0]
        rEnd = int(hit.query_end) + cr[0]

        if rStart > rEnd:
            rStart, rEnd = rEnd, rStart

        expId, runId = runName.split("-")
       
        zscore=int(hit.zScore * 10000) if hit.zScore != None else INT_NODATA
        nMM = hit.nMismatch
        nM = nCorrect
        nIns = hit.nIns
        nDel = hit.nDel

        alignData = {"RefSeq": hit.target_id, "tStart":tStart, "tEnd":tEnd, 
                     "AlignedStrand":relativeStrand, "MovieName":movieName, 
                     "ExpId":expId, "RunId":runId, "Panel":panel, 
                     "x":zmwX, "y":zmwY, "SubreadId":subreadId, 
                     "rStart":rStart, "rEnd":rEnd, "Z":zscore, 
                     "nM":nM, "nMM":nMM, "nIns":nIns, "nDel":nDel, 
                     "alignedRead":alignedRead, "alignedTarget":alignedTarget, 
                     "RunName":runName}

        return alignData

    def writeHit(self, hit):

        self.addNewReference(hit.target_id, annotationString="")
        alignData = self.constructAlignData(hit)
        self.writeAnAlignment(alignData)

    def writeHitsSameReference(self, hits):
        self.addNewReference(hits[0].target_id, annotationString="")
        alignmentDataList = []
        for hit in hits:
            alignmentDataList.append( self.constructAlignData(hit) )
        refSeq = alignmentDataList[0]["RefSeq"]
        readGroup = alignmentDataList[0]["ExpId"] + \
            "-" + alignmentDataList[0]["RunId"]
        movieName = alignmentDataList[0]["MovieName"]
        self.writeAlignments( refSeq, readGroup, movieName, alignmentDataList)

    def writeAnAlignment(self, alignmentData):
        refSeq = alignmentData["RefSeq"]
        readGroup = alignmentData["ExpId"] + "-" + alignmentData["RunId"]
        movieName = alignmentData["MovieName"]
        self.writeAlignments( refSeq, readGroup, movieName, [alignmentData])
 
    def _updateMovieDict(self):
        lastRow = self._movieID.attrs["lastRow"]
        if lastRow != 0:
            ids = self._movieID[:lastRow]
            names = self._movieName[:lastRow]    
            self._movie2Id = dict(zip(names, ids))
        else:
            self._movie2Id = {}

    def addANewMovie(self, movieName):
        lastRow = self._movieID.attrs["lastRow"]
            
        #resize the table if necessary
        movieTableSize = self._movieID.shape[0]
        if lastRow + 1 > movieTableSize: 
            self._movieID.resize( (lastRow+32, ) )
            self._movieName.resize( (lastRow+32, ) )
            
        self._movieID[lastRow] = lastRow + 1
        self._movieName[lastRow] = movieName
        self._movie2Id[movieName] = lastRow + 1
            
        self._movieID.attrs.modify("lastRow", lastRow + 1)
        
            

    def _assignAlignmentIndexCol(self):
        self._indexColNames = ("AlignmentId", "ReadGroupId", "MovieId", 
                               "RefSeqId", "tStart", "tEnd", "AlignedStrand", 
                               "ExpId", "RunId", "Panel", "x", "y", "SubreadId", "rStart", "rEnd",
                               "Z", "nM", "nMM", "nIns", "nDel", "offset_begin", "offset_end")  #For Astros only now Aug. 12, 2009, JC, the ReadGroupId should be replaced with H5R reference.....
        
        self._indexColNum = len(self._indexColNames)
        self._colName2Index = dict(zip(self._indexColNames,range(self._indexColNum)))
   
    def writeAnAlignment(self, alignmentData):
        refSeq = alignmentData["RefSeq"]
        readGroup = alignmentData["ExpId"] + "-" + alignmentData["RunId"]
        movieName = alignmentData["MovieName"]
        self.writeAlignments( refSeq, readGroup, movieName, [alignmentData])
    
    def writeAlignments(self, refSeqName, readGroupName, movieName, alignmentDataList):

        try:
            refSeqGroup = self.getHDF5ObjectByPath(refSeqName)
        except KeyError:
            self.addNewReference(refSeqName)
            refSeqGroup = self.getHDF5ObjectByPath(refSeqName)
        refSeqId = self._refSeqName2Id[refSeqName]
        refSeqIdCol = self._colName2Index["RefSeqId"]
        
        try:
            readGroup = refSeqGroup[readGroupName]
        except KeyError:
            self.addReadGroup(refSeqGroup.name, readGroupName)
            readGroup = refSeqGroup[readGroupName]
            
        if movieName not in self._movie2Id:
            self.addANewMovie(movieName)
        movieId = self._movie2Id[movieName]
        movieIdCol = self._colName2Index["MovieId"]
        
        readGroupId = self._readGroupPath2Id[refSeqGroup.name + "/" + readGroupName ]
        readGroupIdCol = self._colName2Index["ReadGroupId"]
        
        alignmentArray = readGroup["AlignmentArray"]
        alignmentIndex = readGroup["AlignmentIndex"]
        globalIndex = self.getHDF5ObjectByPath("/AlignmentIndex")
        
        alignmentArrayLastRow = alignmentArray.attrs["lastRow"]
        alignmentIndexLastRow = alignmentIndex.attrs["lastRow"]
        globalIndexLastRow = globalIndex.attrs["lastRow"]
        
        alignmentIdCol = self._colName2Index["AlignmentId"]
        
        offset_begin = alignmentArrayLastRow
        offset_end = alignmentArrayLastRow
        offsetBeginColIndex = self._colName2Index["offset_begin"]
        offsetEndColIndex = self._colName2Index["offset_end"]
        
        alignmentIndexData = []
        alignmentData = []
        alignmentId = globalIndexLastRow + 1
        
        for adl in alignmentDataList:
            
            indexData =  numpy.array([ int(adl.get(self._indexColNames[col], numpy.nan)) for col in range(self._indexColNum) ], dtype=numpy.int32)
            
            indexData[movieIdCol] = movieId
            indexData[readGroupIdCol] = readGroupId
            indexData[refSeqIdCol] = refSeqId
            
            alignedRead = adl["alignedRead"]
            alignedTarget = adl["alignedTarget"]
            alignLength = len(alignedRead)  
            
            offset_end = offset_begin + alignLength
            indexData[offsetBeginColIndex] = offset_begin
            indexData[offsetEndColIndex] = offset_end 
            
            indexData[alignmentIdCol] = alignmentId
            alignmentId += 1
            
            alignmentIndexData.append( indexData )
            alignmentData += [ ( ( base2hexMap[alignedRead[x]] << 4) | base2hexMap[alignedTarget[x]] ) for x in xrange(alignLength) ] 
            alignmentData += [0] # zero padded at the end of each of the alignment array
            
            offset_begin = offset_end + 1 # zero padded at the end of each of the alignment array
            
        newIndex = numpy.array(alignmentIndexData)
        newIndexSize = len(newIndex)
        newAlignmentData = numpy.array(alignmentData, dtype=numpy.uint8)
        newAlignmentSize = len(newAlignmentData)
        
        if alignmentIndex.shape[0] < alignmentIndexLastRow + newIndexSize:
            alignmentIndex.resize( (alignmentIndexLastRow + newIndexSize, self._indexColNum ) )
        alignmentIndex[alignmentIndexLastRow:(alignmentIndexLastRow + newIndexSize),...] =  newIndex   
        alignmentIndex.attrs.modify('lastRow', alignmentIndexLastRow + newIndexSize)

        if globalIndex.shape[0] < globalIndexLastRow + newIndexSize:
            globalIndex.resize( (globalIndexLastRow + newIndexSize, self._indexColNum ) )
        globalIndex[globalIndexLastRow:(globalIndexLastRow + newIndexSize),...] =  newIndex   
        globalIndex.attrs.modify('lastRow', globalIndexLastRow + newIndexSize)

        if alignmentArray.shape[0] < alignmentArrayLastRow +  newAlignmentSize:
            alignmentArray.resize( (alignmentArrayLastRow +  newAlignmentSize,) )
        alignmentArray[alignmentArrayLastRow:(alignmentArrayLastRow +  newAlignmentSize)] = newAlignmentData;
        alignmentArray.attrs.modify('lastRow', alignmentArrayLastRow +  newAlignmentSize)

class AlignedAstroReadH5FDAO(object):  #this should be subclassed from a more general interface if necessary
    def __init__(self, h5f):  #h5f should be a instance of AlignedAstroReadH5F, we can force to do a type checking here if necessary
        self.h5f = h5f 
        self._globalIndex = self.h5f.getHDF5ObjectByPath('/AlignmentIndex').value

    def getRefSeqs(self):
        return [ str(seq) for seq in self.h5f.getHDF5ObjectByPath('/RefSeqName').value if seq ]
        

    def _getAlignmentData(self, refSeqGroup, rowIndex):
        
        globalIndex = self._globalIndex 
        
        AlignmentId, ReadGroupId, MovieId,\
        RefSeqId, tStart, tEnd, AlignedStrand,\
        ExpId, RunId, Panel, x, y, SubreadId, rStart, rEnd,\
        Z, nM, nMM, nIns, nDel, offset_begin, offset_end = globalIndex[rowIndex, ...]
        
        runName = "%7d-%04d" % (ExpId, RunId)
        ReadGroup = refSeqGroup[runName]

        alignmentArrayDS = ReadGroup['AlignmentArray']
        alignmentArray = alignmentArrayDS[offset_begin:offset_end] 

        binRBases = (alignmentArray & 0xf0) >> 4; 
        binTBases = (alignmentArray & 0x0f) ;
        
        if AlignedStrand == 0:
            rSeq = "".join(Basemap[binRBases])
            tSeq = "".join(Basemap[binTBases])
        else:
            rSeq = "".join(RBasemap[binRBases[::-1]])
            tSeq = "".join(RBasemap[binTBases[::-1]])

        tPosMap = {}
        offset = 0
        for t in xrange(tStart, tEnd):
            endOffset = offset+1
            while endOffset < len(tSeq) and tSeq[endOffset] == "-":
                endOffset += 1
            #print t, offset, endOffset, rSeq[offset:endOffset], tSeq[offset:endOffset]
            tPosMap[t] = (offset, endOffset)
            offset = endOffset
        readId = "x%d_y%d_%s_p%d_s%d" % (x, y, runName, Panel, SubreadId)
        #readId = "%d" % AlignmentId
        return [ tStart, tEnd, tPosMap, AlignedStrand, rSeq, tSeq, readId ]
  
    def _getPerBaseInfo( self, readGroup ):
        """Returns a dictionary containing all of the per-base 
        info that we will want to return with our alignment hit."""

        if 'AlignmentArray' not in readGroup:
            return None

        alignmentArrayDS = readGroup['AlignmentArray']
        dataSize = len(alignmentArrayDS)
        
        # fetch all to memory for speeding up, it 
        # requires explicitly slicing coordinate to copy the data 
        alignmentArray = alignmentArrayDS[0:dataSize] 
        
        ### these are done in numpy, fast,.,
        binRBases = (alignmentArray & 0xf0) >> 4; 
        binTBases = (alignmentArray & 0x0f) ;
        rSeqAll = "".join(Basemap[binRBases])
        tSeqAll = "".join(Basemap[binTBases])

        return { "tSeq":tSeqAll, "rSeq":rSeqAll }

    def _hitFromPerBaseInfo( self, perBaseInfo, offsetBegin, offsetEnd, hitToModify=None ):
        """Given a dictionary of the per-base metrics, 
        writes the info into an AlignmentHit object
        (or subclass) and returns it.
        
        If hitToModify is specified, that hit is modified 
        and returned, instead of generating a new hit."""
 
        hit = hitToModify if hitToModify else AlignmentHit()
        
        hit.alignedQuery = perBaseInfo["rSeq"][offsetBegin:offsetEnd]
        hit.alignedTarget = perBaseInfo["tSeq"][offsetBegin:offsetEnd]
        hit.hasAlignment = True

        return hit

    def getAlignedReadIterator(self, refSeqName):
        refSeqName = refSeqName.strip("/")
        refSeqGroup = self.h5f.getHDF5ObjectByPath('/'+refSeqName)
        theRefSeqId = self.h5f._refSeqName2Id[refSeqName]
        mId2movieName = self.h5f._movie2Id
        movieName2mId = dict(zip(mId2movieName.values(), mId2movieName.keys()))
       

        for readGroupName, readGroup in refSeqGroup.iteritems():
            lastRow = readGroup['AlignmentIndex'].attrs['lastRow']
            indexArray = readGroup['AlignmentIndex'][0:lastRow]
            perBaseInfo = self._getPerBaseInfo( readGroup )

            for rowIndex in xrange(0, lastRow):
                rowData = indexArray[rowIndex]
                AlignmentId, ReadGroupId, MovieId,\
                RefSeqId, tStart, tEnd, AlignedStrand,\
                ExpId, RunId, Panel, x, y, SubreadId, rStart, rEnd,\
                Z, nM, nMM, nIns, nDel, offset_begin, offset_end = rowData

                runName = "%7d-%04d" % (ExpId, RunId)
                
                movieName = movieName2mId[MovieId]
                readId = "x%d_y%d_%s_%s_s%d" % (x, y, runName, movieName, SubreadId)
                hit = self._hitFromPerBaseInfo( perBaseInfo, offset_begin, offset_end )
                                                 
                hit.setMatchInfo(nM, nMM, nIns, nDel, tEnd-tStart) #hit.alignedLength = tEnd- tStart? or len(tSeq)?, JC

                #hit.setMatchInfo(nM, nMM, nIns, nDel, tEnd-tStart, rSeq, tSeq) #hit.alignedLength = tEnd- tStart? or len(tSeq)?, JC
                hit.query_id = readId
                hit.query_start = rStart
                hit.query_end = rEnd
                hit.query_strand = "+" 
                hit.query_length = rEnd - rStart

                hit.target_id = refSeqName
                hit.target_start = tStart
                hit.target_end = tEnd
                hit.target_strand = "+" if AlignedStrand == 0 else "-"
                hit.target_length = tEnd - tStart

                hit.zScore = 1.0*Z/10000 if Z != INT_NODATA else None
                yield hit 



    def getColumnIterator(self, refSeqName, refPosition):
        return self._getColumnIterator(refSeqName, refPosition, exportRefOffset=False)
     
    def getRefseqOffsetColumnIterator(self, refSeqName, refPosition):
        return self._getColumnIterator(refSeqName, refPosition, exportRefOffset=True)
     
    def _getColumnIterator(self, refSeqName, refPosition, exportRefOffset=False):
        globalIndex = self._globalIndex 
        tStartCol = self.h5f._colName2Index["tStart"]
        tEndCol = self.h5f._colName2Index["tEnd"]
        refSeqIdCol = self.h5f._colName2Index["RefSeqId"]
        
        sortI = numpy.argsort(globalIndex[...,tStartCol])
        sortStartT = numpy.sort(globalIndex[...,tStartCol])
        endT = [globalIndex[...,tEndCol][x] for x in sortI]
        refSeqId = [globalIndex[...,refSeqIdCol][x] for x in sortI]

        refSeqName = refSeqName.strip("/")
        refSeqGroup = self.h5f.getHDF5ObjectByPath('/'+refSeqName)
        theRefSeqId = self.h5f._refSeqName2Id[refSeqName]
        sortedIndex = [x for x in zip(sortI, sortStartT, endT, refSeqId ) if x[3] == theRefSeqId]
       
        s1 = bisect.bisect_left(sortStartT, refPosition)
        sortedIndex2 = sortedIndex[:s1]
        sortedIndex2.sort(key=lambda x:x[2])
        keys = [x[2] for x in sortedIndex2]
        s2 = bisect.bisect_left(keys, refPosition+1)
        
        readIds = sortedIndex2[s2:]  #all reads end after the currentRefPos
        #all reads starts at the currentRefPos
        while s1 < len(sortedIndex) and sortedIndex[s1][1] <= refPosition: 
            readIds.append(sortedIndex[s1])
            s1+=1
            
        currentReads = []
        for rId in readIds:  
            currentReads.append( self._getAlignmentData( refSeqGroup, rId[0] ))
        
        currentRefPos = refPosition
        currentInsertOffset = 0
        while 1:
            alignColumn = []
            moreInsertions = 0
            for r in currentReads:
                tStart, tEnd, tPosMap, AlignedStrand, rSeq, tSeq, readId = r
                offset = tPosMap[currentRefPos][0]
                endOffset = tPosMap[currentRefPos][1]
                readBases = rSeq[offset:endOffset]
                #print len(readBases), currentInsertOffset, endOffset, offset
                if currentInsertOffset < len(readBases):
                    moreInsertions = 1
                    alignColumn.append( (readBases[currentInsertOffset], AlignedStrand, readId))
                else:
                    alignColumn.append( ('-',  r[3], readId))
                    
            if moreInsertions == 1:
                if exportRefOffset:
                    yield alignColumn, currentRefPos, currentInsertOffset
                else:
                    yield alignColumn
                currentInsertOffset += 1
            else:
                currentRefPos += 1
                currentInsertOffset = 0
                currentReads = [ r for r in currentReads if r[1] > currentRefPos ] #faster way to remove reads?
                
                readIds = []  #all reads end after the currentRefPos
                while s1 < len(sortedIndex) and sortedIndex[s1][1] <= currentRefPos: 
                    readIds.append(sortedIndex[s1])
                    s1+=1

                for rId in readIds:  
                    currentReads.append( self._getAlignmentData( refSeqGroup, rId[0] ))

                if len(currentReads) == 0 and s1 >= len(sortedIndex):
                    break

class testAlignedReadH5F(unittest.TestCase):
    def setUp(self):
        self.testDirectory = "../../test/data/"

    def test_addNewReference(self):
        self.h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_1.h5","w")
        for i in range(10):
            self.h5f.addNewReference("testRefSeq%03d" % i)
        self.h5f._h5FileHandle.close()
 
    def test_addNewReadGroup(self):
        self.h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_2.h5","w")
        self.h5f.addNewReference("testRefSeq1", annotationString="Hello World")
        g = self.h5f.getHDF5ObjectByPath("/testRefSeq1")
        for i in range(10):
            self.h5f.addReadGroup(g.name, "readGroup%03d" % i)
            rg = self.h5f.getHDF5ObjectByPath(g.name+"/readGroup%03d" % i)
        self.h5f._h5FileHandle.close()

    def test_writeAlignmentReads(self):
        self.h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_3.h5","w")
        self.h5f.addNewReference("testRefSeq1", annotationString="Hello World")
        g = self.h5f.getHDF5ObjectByPath("/testRefSeq1")
        self.h5f.addReadGroup(g.name, "readGroup01")
        
        l = """x29_y54,2000049-0006,"m090718_111806_Uni_p2_b20.f1",+,0,526,-,27997,27426,4.191,25,70,107,324,GTGAGGGGC-AGGCAA-ATTAGAAA-AATTT-GGAGATCAGGGTTTGAGTTTTTTAACTATAGCT-GATTTC-ACTCATTTTCTAATTTTTTATAA-CAAACTAAAAAT--ATAAAACTCC-GAATC-T-A--TTG-A-CCTTAAAAC-A-CTCGG-AGTA--TCCACT-TTAAC---CTTTCTTGTTTTTTTCGGTATAACA--TCCG-TTG-TGGTC-TTATAACAAT--CCTTTT-G-ACC----GCTAACTTTGTTT-TAAAGT-TTAA-ACATAAGAACTTGCACCCTAATGTTTCCTCGGA--T-TTCGCTCTTTT--T-GGTAA--GATTATC----TATCTT--CGTC-CCA-CAAAACGAGGAGGGAGAAAAAAATAAGAAAAAATAAGCTTAACCACAAAGAT-CGAACATTCCAGGGAACATTAGTGTGCAAGTCTATCTACATTATAAAAAAGTCACGGGACAA-AGGTTATTTTATAAGACCAC-GAATT-GTTAACATCAGCGGCGCGCGCGGCAA---ACTGGTTGCTAATAGCGTTGCCGGACCAA---AGACAACTTAACT-CGA--ATGATGTCCGGT,GGGTGGGGCGAAGCAACATT-GAAACAATTTTGGAGAT-AGG-TTTGAGTTTT--AACTATAGCTTGATTTCCAATCATTTTCTAATTTTTTATAATCAAACTAAAAATTAATAAAACTCGAGAATCCTGACTTTGGAGCCTTTAAACTAGCTCGGCAGTACTTTCACTATTAACACCCTTTCTTGTTTTT---GGTTTAACAATTCCGGTTGGTGGTCCTTATAACAATAACCTTTTTGTACCCCCCGCTGACTTTTTTAATAAAGTGTTAAGAAATAAGA-CTTGCACACTAATGTTTCCTGGGAACTGTTCGATCTTTTGGTAGCTAACTTTTTATCCTTTTATTTTTACTTCTTCATTATAACCCTTA-TCAAACTTACTTTTTAAAAAACGCGCATTGCTCCAAA-ATCCGTA-ATGCTGTAGATCTCCAG-CAGC--GTGTAGATA-GTTCGA----AGACCC--GACAATTTCTCAATTTACCCGACTTCAAAATTCGTTTTAATTCGCTTTCAATACGGCCACCCACTG--TGCAAAAAGCGAT-CGAAAAGAAGCGAGGGAAAGTAGTTGAGAATATGATGTCCAGT"""
        l = l.split(",")
        zmwID,runName,movieName,\
        strand,start,end,\
        targetStrand,targetStart,targetEnd,\
        Z,nInsert,nDelete,nMismatch,nCorrect,\
        alignedRead,alignedTarget = l 

        if "." in movieName:
            movieName, subreadId = movieName.strip('"').split(".")
            subreadId = subreadId[1:]
        else:
            subreadId = 0
            
        zmwX,zmwY = zmwID.split("_")
        zmwX=zmwY[1:]
        zmwY=zmwY[1:]
        panel = movieName.split("_")[3][1:]

        if strand == "-":
            alignedRead = revCompSeq(alignedRead)
            alignedTarget = revCompSeq(alignedTarget)
        if strand != targetStrand:
            relativeStrand = 1
        else:
            relativeStrand = 0
        expId = runName.split("-")[0]
        runId = runName.split("-")[1]
        zscore=int(float(Z)*10000)
        nMM = nMismatch
        nM = nCorrect
        nIns = nInsert
        nDel = nDelete
        alignData = {"RefSeq":"testRefSeq1", 
                     "tStart":targetStart, "tEnd":targetEnd, "AlignedStrand":relativeStrand, 
                     "ExpId":expId, "RunId":runId, "Panel":panel, "x":zmwX, "y":zmwY, "SubreadId":subreadId, 
                     "rStart":start, "rEnd":end,
                     "Z":zscore, "nM":nM, "nMM":nMM, "nIns":nIns, "nDel":nDel, 
                     "alignedRead":alignedRead, "alignedTarget":alignedTarget}
        
        #self.h5f.writeAnAlignment(alignData)
        
        self.h5f.writeAlignments( "testRefSeq1", runName, movieName, [alignData]*100)
                
        self.h5f._h5FileHandle.close()
        
    def emitMSA(self, reads, tpos):
        reads = [ r for r in reads if r[1] > tpos]
        coverage = len(reads)
        baseCount = {}
        insertBaseCount = {}
        if len(reads) == 0: return reads
        for r in reads:
            rSeq = r[4]
            tSeq = r[5]
            tPosMap = r[2]
            offset = tPosMap[tpos][0]
            endOffset = tPosMap[tpos][1]
            readBases = rSeq[offset:endOffset]
            alignedBase = readBases[0]
            baseCount[alignedBase] = baseCount.get(alignedBase, 0) + 1
            for c in readBases[1:]:
                insertBaseCount[c] = insertBaseCount.get(c,0) + 1
            #print tpos, rSeq[offset:endOffset], tSeq[offset:endOffset]
        tbase = tSeq[offset]
        for b in ["A","C","G","T","-"]:
            print tpos, tbase, b, baseCount.get(b, 0), coverage
        for b in ["A","C","G","T"]:
            print tpos, "-", b, insertBaseCount.get(b, 0), coverage
        return reads

    def test_writeCompareCSVReads(self):
        
        h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_4.h5","w")
        h5f.addNewReference("testRefSeq1", annotationString="Hello World")
        g = h5f.getHDF5ObjectByPath("/testRefSeq1")
        f = open(self.testDirectory+"AlignedReadH5IO_test.csv")
        allAlignedReadData = {}
        for l in f:
            l = l.split(",")
            if l[0] == "zmwID": continue
            zmwID,runName,movieName,\
            strand,start,end,\
            targetStrand,targetStart,targetEnd,\
            Z,nInsert,nDelete,nMismatch,nCorrect,\
            alignedRead,alignedTarget = l 
    
            if "." in movieName:
                movieName, subreadId = movieName.strip('"').split(".")
                subreadId = subreadId[1:]
            else:
                subreadId = 0
                
            zmwX,zmwY = zmwID.split("_")
            zmwX=zmwX[1:]
            zmwY=zmwY[1:]
            panel = movieName.split("_")[3][1:]
    
            if strand == "-":
                alignedRead = revCompSeq(alignedRead)
                alignedTarget = revCompSeq(alignedTarget)
            if strand != targetStrand:
                relativeStrand = 1
            else:
                relativeStrand = 0
            
            tStart = int(targetStart)
            tEnd = int(targetEnd)
            if tStart > tEnd:
                tStart, tEnd = tEnd, tStart
    
            rStart = int(start)
            rEnd = int(end)
            if rStart > rEnd:
                rStart, rEnd = rEnd, rStart

                
            expId = runName.split("-")[0]
            runId = runName.split("-")[1]
            zscore=int(float(Z)*10000) 
            nMM = nMismatch
            nM = nCorrect
            nIns = nInsert
            nDel = nDelete
            alignData = {"RefSeq":"testRefSeq1", 
                         "tStart":tStart, "tEnd":tEnd, "AlignedStrand":relativeStrand, 
                         "ExpId":expId, "RunId":runId, "Panel":panel, "x":zmwX, "y":zmwY, "SubreadId":subreadId, 
                         "rStart":rStart, "rEnd":rEnd,
                         "Z":zscore, "nM":nM, "nMM":nMM, "nIns":nIns, "nDel":nDel, 
                         "alignedRead":alignedRead, "alignedTarget":alignedTarget}
            allAlignedReadData.setdefault( (runName,movieName),[])
            allAlignedReadData[(runName,movieName)].append(alignData)
        #self.h5f.writeAnAlignment(alignData)
        
        for rn in allAlignedReadData:
            h5f.writeAlignments( "testRefSeq1", rn[0], rn[1], allAlignedReadData[rn])
        
        globalIndex = h5f.getHDF5ObjectByPath('/AlignmentIndex').value
        sortI = numpy.argsort(globalIndex[...,4])
        currentReads = []
        currentTPos = -1
        refSeqGroup = h5f.getHDF5ObjectByPath('/testRefSeq1')
        theRefSeqId = h5f._refSeqName2Id['testRefSeq1']
        
        for i in sortI:
            
            AlignmentId, ReadGroupId, MovieId,\
            RefSeqId, tStart, tEnd, AlignedStrand,\
            ExpId, RunId, Panel, x, y, SubreadId, rStart, rEnd,\
            Z, nM, nMM, nIns, nDel, offset_begin, offset_end = globalIndex[i,...]

            if RefSeqId != theRefSeqId: continue
            #if alignedStrand == 0: continue
            if currentTPos != -1 and tStart != currentTPos:
                for tpos in range(currentTPos, tStart):
                    currentReads = self.emitMSA(currentReads, tpos)
            currentTPos = tStart
            runName = "%7d-%04d" % (ExpId, RunId)
            ReadGroup = refSeqGroup[runName]
            dataSet = ReadGroup['AlignmentArray']
            if AlignedStrand == 0:
                alignData = [alignmentPairMap[c] for c in dataSet[offset_begin:offset_end]]
            else:
                alignData = [rAlignmentPairMap[c] for c in dataSet[offset_begin:offset_end][::-1]]
    
            rSeq,tSeq = zip(*alignData)
            rSeq, tSeq = "".join(rSeq), "".join(tSeq)
            tPosMap = {}
            offset = 0
            for t in xrange(tStart, tEnd):
                endOffset = offset+1
                while endOffset < len(tSeq) and tSeq[endOffset] == "-":
                    endOffset += 1
                #print t, offset, endOffset, rSeq[offset:endOffset], tSeq[offset:endOffset]
                tPosMap[t] = (offset, endOffset)
                offset = endOffset
    
            currentReads.append( [tStart, tEnd, tPosMap, AlignedStrand, rSeq, tSeq ] )
                              
        h5f._h5FileHandle.close()
        

class TestH5FDAO(unittest.TestCase):
    def setUp(self):
        self.testDirectory = "../../test/data/"

    def test_read(self):
        h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_4.h5","r")
        dao = AlignedAstroReadH5FDAO(h5f)
        colIterator = dao.getColumnIterator('testRefSeq1',0)
        efaStrings = {}
        i = 0
        while 1:
            try:
                col = colIterator.next()
            except StopIteration:
                break
            if i > 1000: break #for testing purpose only output 1000 column
            for elm in col:
                #print elm
                efaStrings.setdefault(elm[2], (i,[]) )
                efaStrings[elm[2]][1].append(elm[0])
            i += 1
            
        for tag in efaStrings:
            print ">"+tag+"_%d" % efaStrings[tag][0]
            print "".join(efaStrings[tag][1])

    def test_pairwiseAlignmentIterator(self):
        h5f = AlignedAstroReadH5F(self.testDirectory+"ARTest_4.h5","r")
        dao = AlignedAstroReadH5FDAO(h5f)
        arIter = dao.getAlignedReadIterator('testRefSeq1')
        for alignment in arIter:
            print alignment
            print alignment.alignedQuery
            print alignment.alignedTarget
    
if __name__ == "__main__":
    suite1 = unittest.TestLoader().loadTestsFromTestCase(testAlignedReadH5F)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(TestH5FDAO)
    unittest.TextTestRunner(verbosity=2).run(suite1)
    unittest.TextTestRunner(verbosity=2).run(suite2)
    # # a.test_writeCompareCSVReads()
    # unittest.TextTestRunner(verbosity=2).run(suite)

