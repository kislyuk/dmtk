"""
a prototype of the alignment storage model in hdf5
"""

from tables import *
import numpy  
import bisect

class SeqIndexTable(IsDescription):
    seqId = StringCol(128, pos=1)
    seqRangeStart = Int32Col(pos=2) 
    seqRangeEnd =  Int32Col(pos=3)
    seqRangeFirst = Int32Col(pos=4)
    seqRangeLast = Int32Col(pos=5)
    seqReadTableName = StringCol(128, pos=5)

class ReadIndexTable(IsDescription):
    hitId = Int32Col(pos=1)
    readName = StringCol(64, pos=2)
    readLength = Int32Col(pos=3)
    alignScore = Int16Col(pos=4)
    queryStart = Int32Col(pos=5)
    queryEnd = Int32Col(pos=6)
    targetStart = Int32Col(pos=7)
    targetEnd = Int32Col(pos=8)
    nIns = Int32Col(pos=9)
    nDel = Int32Col(pos=10)
    nMM = Int32Col(pos=11)
    #tableStartRow = Int32Col(pos=12)
    #tableEndRow = Int32Col(pos=13)
    
class SeqReadTable(IsDescription):
    hitId = Int32Col(pos=1)
    readPos = Int32Col(pos=2)
    readBase = StringCol(1,pos=3)
    seqPos = Int32Col(pos=4)
    seqBase = StringCol(1,pos=5)
    alignAnnotation = StringCol(1,pos=6)
    alignOrientation = StringCol(1,pos=7)
    
    
class SeqReadPulseFeatureTable(IsDescription):
    #TODO: is there a better way to do this? Class inheritance is broken for IsDescription class 
    hitId = Int32Col(pos=1)
    readPos = Int32Col(pos=2)
    readBase = StringCol(1,pos=3)
    seqPos = Int32Col(pos=4)
    seqBase = StringCol(1,pos=5)
    alignAnnotation = StringCol(1,pos=6)
    alignOrientation = StringCol(1,pos=7)
    
    t1 = Float32Col(pos=8)
    dt = Float32Col(pos=9)
    ts = Float32Col(pos=10)
    PNR = Float32Col(pos=11)
    SNR = Float32Col(pos=12)
    minChi2dof = Float32Col(pos=13)
    minDeltaChi2dof = Float32Col(pos=14)
    pkmax = Float32Col(pos=15)
    Zpkmax = Float32Col(pos=16)
    xcen = Float32Col(pos=17)
    xcenDYE = Float32Col(pos=18)



class AlignedSeqRead(object):

    alignMetaFields = ('readName', 'readLength', 'alignScore',
                       'queryStart', 'queryEnd', 'targetStart', 'targetEnd',
                       'nIns', 'nDel', 'nMM')
    
    alignDataFields = ('readPos', 'readBase', 'seqPos', 'seqBase', 'alignAnnotation', 'alignOrientation')
    seqReadTableType = SeqReadTable
    
    def __init__(self, alignMeta, alignData):
        
        for fieldName in self.alignMetaFields:
            self.__dict__[fieldName] = alignMeta[fieldName]
            
        if self.targetStart > self.targetEnd:
            self.alignOrientation = "-"
        else:
            self.alignOrientation = "+"
        self.alignData = {}
        dataLength = -1
        
        for fieldName in self.alignDataFields:
            self.alignData[fieldName] = alignData[fieldName]
            if dataLength == -1:
                dataLength = len(alignData[fieldName])
            elif dataLength != len(alignData[fieldName]):
                raise DataError('All data array in alignData should have the same length.')
        self.dataLength = dataLength
    
    def getMetaDataRow(self, hitId = 0):
        row = [hitId]
        for fieldName in self.alignMetaFields:
            row.append(self.__dict__[fieldName])
        #row = row + [-1, -1] #for table starting row and table ending row
        return row
    
    def getSeqDataTable(self, hitId = 0, fields = None):
        if fields == None:
            fields = self.alignDataFields
        allData = []
        for i in xrange(self.dataLength):
            data = [hitId]
            for f in fields:
                data.append(self.alignData[f][i])
            allData.append(tuple(data))
        return allData

    def getSeqDataTableInRange(self, hitId = 0, seqRange = (0,0) ,fields = None):
        if fields == None:
            fields = self.alignDataFields
        allData = []
        ss, se = seqRange
        for i in xrange(self.dataLength):
            data = [hitId]
            seqPos = self.alignData['seqPos'][i]
            if ss <= seqPos and seqPos < se: #TODO: check this carefully
                for f in fields:
                    data.append(self.alignData[f][i])
                allData.append(tuple(data))
        return allData

    
class AlignedSeqReadWithPulseFeature(AlignedSeqRead):
    alignDataFields = ('readPos', 'readBase', 'seqPos', 'seqBase', 'alignAnnotation', 'alignOrientation',
                       't1', 'dt', 'ts',
                       'PNR', 'SNR', 'minChi2dof', 'minDeltaChi2dof',
                       'pkmax', 'Zpkmax', 'xcen', 'xcenDYE')
    seqReadTableType = SeqReadPulseFeatureTable
    
    
class AlignedReadStorage(object):
    
    def __init__(self, hdf5filename, mode ='a'):
        self.hdf5filename = hdf5filename
        if mode in ['a','w','r']:
            self.hdf5handle = openFile(self.hdf5filename, mode = mode, title = "PBI Alignment Store Prototype")
        else:
            raise IOError("Ivalid mode to open a HDF5 file: mode %s" % mode)

    def createTableForSeq(self, seqName, seqRanges=None, seqReadTableType=AlignedSeqReadWithPulseFeature.seqReadTableType):
        seqGroup = self.hdf5handle.createGroup("/", seqName, "all reads about the sequence " + seqName)
        self.hdf5handle.createTable("/"+seqName, seqName+"_readIndex", ReadIndexTable, "Read Index")
        self.hdf5handle.createTable("/"+seqName, seqName+"_seqIndex", SeqIndexTable, "Sequence Index")
        
        for r in xrange(len(seqRanges)-1):
            r = tuple(seqRanges[r:r+2])
            tableName = seqName+"_%08d_%08d" % r
            self.hdf5handle.createTable("/"+seqName, tableName, 
                                        seqReadTableType, "Seq Read Table %d to %d " % r,
                                        expectedrows=50000)
            t = self.hdf5handle.getNode("/"+seqName, tableName)
            t.attrs.seqSorted = 'False'
            t.attrs.readSorted = 'False'
            seqIndexTable = self.hdf5handle.getNode("/"+seqName, seqName+"_seqIndex")
            seqIndexTable.append([(seqName, r[0], r[1], -1, -1, tableName)])
        seqIndexTable.flush()

    def getSeqIndexTable(self, seqName, seqStart):
        seqIndexTable = self.hdf5handle.getNode("/"+seqName, seqName+"_seqIndex");

        targetTableName = [ (x['seqReadTableName'],
                             x['seqRangeStart'],
                             x['seqRangeEnd']) for x in seqIndexTable.where( 
                            '(seqRangeStart <= %d) & (seqRangeEnd > %d)' % (seqStart, seqStart) ) ]

        if len(targetTableName) > 1:
            #FIXME: do the error handling properly
            raise "Multiple tables for a Read: %s, %d" %( targetTableName, seqStart )
        elif len(targetTableName) == 1:
            return targetTableName[0]
        else:
            return None

    def getSeqIndexTables(self, seqName, seqStart, seqEnd):
        seqIndexTable = self.hdf5handle.getNode("/"+seqName, seqName+"_seqIndex");
        #targetTableName = [ x['seqReadTableName'] for x in seqIndexTable.iterrows() 
        #                    if x['seqRangeStart'] <= seqStart and x['seqRangeEnd'] > seqStart ]
        
        
        if seqStart > seqEnd:
            seqStart, seqEnd = seqEnd+1, seqStart+1  #TODO: check this carefully
        targetTableName = [ (x['seqReadTableName'],
                             x['seqRangeStart'],
                             x['seqRangeEnd']) for x in seqIndexTable.where( 
                            '(seqRangeStart < %d) & (seqRangeEnd > %d)' % (seqEnd, seqStart) ) ]
        return targetTableName
        
    def writeARead(self, seqName, alignedSeqRead):
        #TODO: determine which table to write
        table = self.hdf5handle.getNode("/"+seqName, seqName+"_readIndex");
        hitId = table.nrows + 1
        table.append([alignedSeqRead.getMetaDataRow(hitId)])
        table.flush()

        seqStart = alignedSeqRead.targetStart
        seqEnd = alignedSeqRead.targetEnd

        if seqStart > seqEnd:
            seqStart, seqEnd = seqEnd, seqStart

        targetTableNames = self.getSeqIndexTables(seqName, seqStart, seqEnd)

        for tableName, ss, se in targetTableNames:
            table = self.hdf5handle.getNode("/"+seqName, tableName)
            #print hitId, tableName, ss, se, seqStart, seqEnd
            table.append(alignedSeqRead.getSeqDataTableInRange(hitId, (ss, se) ) )
            table.flush()
            table.attrs.seqSorted = 'False'
            table.attrs.readSorted = 'False'
        
    def sortSeqReadTableBySeqPos(self, seqName, seqPos):
        targetTableName = self.getSeqIndexTable(seqName, seqPos)
        table = self.hdf5handle.getNode("/"+seqName+"/"+targetTableName);
        d = table[:]
        d.sort(order=['seqPos','hitId'])
        table[:] = d
        table.flush()
    
    def _sortTableBySeqPos(self, table):
        d = table[:]
        if len(d) != 0 and table.attrs.readSorted == 'False':
            d.sort(order=['seqPos','hitId'])
            table[:] = d
            table.flush()
            table.attrs.readSorted = 'False'
            table.attrs.seqSorted = 'True'    
        table.flush()
        
    def _sortTableByRead(self, table):
        d = table[:]
        if len(d) != 0 and table.attrs.readSorted == 'False':
            d.sort(order=['hitId','readPos'])
            table[:] = d
            table.flush()
            table.attrs.readSorted = 'True'
            table.attrs.seqSorted = 'False'
        table.flush()

    def sortAllSeqReadTableBySeqPos(self, seqName):
        seqReadTableIndex = self.hdf5handle.getNode("/"+seqName, seqName+"_seqIndex");
        for tableName in [x['seqReadTableName'] for x in seqReadTableIndex.iterrows()]:
            table = self.hdf5handle.getNode("/"+seqName, tableName);
            self._sortTableBySeqPos(table)

    def sortAllSeqReadTableByRead(self, seqName):
        seqReadTableIndex = self.hdf5handle.getNode("/"+seqName, seqName+"_seqIndex");
        for tableName in [x['seqReadTableName'] for x in seqReadTableIndex.iterrows()]:
            table = self.hdf5handle.getNode("/"+seqName, tableName);
            self._sortTableByRead(table)

    def getAlignmentIterator(self, seqName, startPos, endPos=-1, selectedHitIds = None):
        #TODO: need to work so the iterator can go across table boundary

        seqPos = startPos
        comMap = dict(zip('ACGT-','TGCA-'))
        tableInfo = self.getSeqIndexTable(seqName, seqPos)
        if tableInfo == None:
            raise StopIteration
        
        tableName, seqRangeStart, seqRangeEnd = tableInfo
        table = self.hdf5handle.getNode("/"+seqName, tableName);
        if table.attrs.seqSorted == 'False':
            self._sortTableBySeqPos(table)

        seqPosArray = table.cols.seqPos[:]

        while 1:
            #left, right = seqPosArray.searchsorted(seqPos, 'left'), seqPosArray.searchsorted(seqPos, 'right')
            
            left, right = bisect.bisect_left(seqPosArray, seqPos),  bisect.bisect_right(seqPosArray, seqPos)
            if left == right:
                if left == 0 and len(seqPosArray) != 0:
                    seqPos = min(seqPosArray)
                    continue
            
                elif right == len(seqPosArray) or len(seqPosArray) == 0:
                    seqPos = seqRangeEnd
                    tableInfo = self.getSeqIndexTable(seqName, seqPos)
                    if tableInfo == None:
                        raise StopIteration
    
                    tableName, seqRangeStart, seqRangeEnd = tableInfo
                    table = self.hdf5handle.getNode("/"+seqName, tableName);
                    
                    if table.attrs.seqSorted == 'False':
                        self._sortTableBySeqPos(table)
                    seqPosArray = table.cols.seqPos[:]
                    continue
                else:
                    seqPos += 1
                    continue
                
            if endPos > 0 and seqPos > endPos:
                raise StopIteration

            rIdToBases = {}
            for d in table[left:right]:
                hitId, readPos, readBase, seqPos, seqBase, annotation, orientation = d
                if selectedHitIds != None and hitId not in selectedHitIds: continue
                if orientation == '-':
                    readBase = comMap[readBase]
                    seqBase = comMap[seqBase]
                if 'ref' not in rIdToBases and seqBase != '-':
                    rIdToBases['ref'] = [seqBase, '+']
                #if 'ref' not in rIdToBases:
                #    rIdToBases['ref'] = [[seqBase], ['+']]

                rIdToBases.setdefault(hitId, [[],[]])
                rIdToBases[ hitId ][0].append(readBase)
                rIdToBases[ hitId ][1].append(orientation)

                #rIdToBases[ hitId ] = (readBase, orientation) 
            maxInsLen = max( [len(x[0]) for x in rIdToBases.values() ]) 
            for rId, data in rIdToBases.items():
                bases, orientation = data
                l = len(bases)
                rIdToBases[rId] = (seqPos, '-' * (maxInsLen - l) + "".join(bases), orientation)
            yield rIdToBases
            seqPos += 1
            
    def getMSA(self, seqName, startPos, endPos=-1, hitIds = None):
        rMap = dict(zip('ATGCatgc-','TACGtagc-'))
        MSASeqs = {}
        aIterator = self.getAlignmentIterator(seqName, startPos, endPos, selectedHitIds = hitIds)
        for x in aIterator:
            for id, alignedSeq in x.items():
                MSASeqs.setdefault(id, []).append(alignedSeq)
        
        hitId2readTag = {}
        readIndex = self.hdf5handle.getNode("/"+seqName, seqName+"_readIndex");
        
        for row in readIndex:
            hitId2readTag[row['hitId']] = row['readName']
        
        seqPos2alignPos = []
        alignPos = 0
        for seqs in MSASeqs['ref']:
            seqPos = seqs[0]
            alignPos += len(seqs[1])
            seqPos2alignPos.append( (seqPos, alignPos) )
        seqPos2alignPos = dict(seqPos2alignPos)
        
        rtnList =[]
        for id,seqs in MSASeqs.items():
            if id != 'ref':
                name = ">%s|start=%d|rId=%d" % (hitId2readTag[id],
                                        seqPos2alignPos[seqs[0][0]] - len( seqs[0][1])+1, id )
            else:
                name = ">%s" % seqName
            #print  "".join([x[1] for x in seqs])
            msaSeq = "-"*(seqPos2alignPos[seqs[0][0]] - len( seqs[0][1])) + "".join([x[1] for x in seqs])
            rtnList.append( (name, msaSeq) )
        return rtnList
    
    def pluralityCall(self, seqName, start = 0, hitIds = None):
        rMap = dict(zip('ATGCatgc-','TACGtagc-'))
        table = self.hdf5handle.getNode("/"+seqName, seqName+"_readIndex");
        hitId = table.cols.hitId[:]
        readName = table.cols.readName[:]
        id2name = dict(zip(hitId,readName))
        
        aIterator = self.getAlignmentIterator(seqName, start, selectedHitIds = hitIds)
        seq = []
        refSeq  =[]
        countM = []
        alignedPos = []

        for x in aIterator:
            alignedBases = zip(*[ y[1][1][:] for y in x.items() if y[0] != 'ref' ])
            #print x
            #alignedDirection =zip(* [ y[1][2] for y in x.items() if y[0] != 'ref' ])
            for i  in xrange(len(alignedBases)):
                bases = alignedBases[i]
                alignedPos.append(x['ref'][0])
                #direction = alignedDirection[i]
                bCount = {}
                for b in bases:
                    bCount[b] = bCount.get(b,0) + 1
                ca = [ z[::-1] for z in bCount.items() ]
                ca.sort()  #TODO: check if there is tie....
                ca = ca[::-1]
                seq.append(ca[0][1])
                countM.append(bCount)
                refSeq.append( x['ref'][1][i] )
                #print [(id2name[rid],x[rid]) for rid in x.keys() if rid != 'ref']
                #print x['ref'][0], x['ref'][1][i], ca[0][1], "".join(bases)
                #print x['ref'][0], x['ref'][1][i], ca[0][1], "".join(alignedDirection)
        tmp=[x for x in zip(seq, refSeq, countM, alignedPos) if x[0]!='-' or x[1] != '-' ]
        return ["".join([x[0] for x in tmp]) , "".join([x[1] for x in tmp]), [x[2] for x in tmp], [x[3] for x in tmp] ]
                    
    def close(self):
        self.hdf5handle.close()
    
def test():
    import psyco
    #psyco.bind(AlignedReadStorage)
    #psyco.bind(AlignedSeqRead)
    #ReadTableType = AlignedSeqReadWithPulseFeature
    ReadTableType = AlignedSeqRead
    
    aStore = AlignedReadStorage("../data/test.h5", 'w')
    aStore.createTableForSeq("testSeq", (0,2000),  seqReadTableType=ReadTableType.seqReadTableType)
    aStore.close()
    
    aStore = AlignedReadStorage("../data/test.h5", 'a')  
    dataFile = open("../data/1590017-0010_Jackalope_p1_ana.csv")
    dataFile.readline()
    data = dataFile.read().split("\n")
    dataFile.close()
    
    indexFileName = "../data/1590017-0010_Jackalope_p1_ana_index.csv"
    try:
        indexFile = open(indexFileName)
    except IOError:
        index = {}
        return
        
    index = {}
    for l in indexFile:
        l = l.strip()
            
        if l[0] == "#":
            continue
        else:
            (idnum, tag, x, y, readLength, score,
            startTime, stopTime,
            queryStart, queryEnd,
            targetStart, targetEnd,
            MMCount, delCount, insCount, 
            startLine, dataLength) = l.split(',')
            index[int(idnum)] =  (tag,int(x), int(y), int(readLength), int(score),
                                  float(startTime), float(stopTime),
                                  int(queryStart), int(queryEnd),
                                  int(targetStart), int(targetEnd),
                                  int(MMCount), int(delCount), int(insCount), int(startLine), int(dataLength))
            readName = tag+"_%s_%s" % (x,y)
            readMeta = {'readName':readName,
                        'readLength':int(readLength),
                        'alignScore':int(score),
                        'queryStart':int(queryStart),
                        'queryEnd':int(queryEnd),
                        'targetStart':int(targetStart),
                        'targetEnd':int(targetEnd),
                        'nIns':int(insCount),
                        'nDel':int(delCount),
                        'nMM':int(MMCount)}
            readData = {}
            orientation = "+" if int(targetStart) < int(targetEnd) else "-"
            for l in data[int(startLine)-1:(int(startLine)+int(dataLength)-1)]:
                (alignmentId, qindex, qbase, tindex, tbase, alignAnnotation, contMatch, 
                 t1, dt, ts, PNR, SNR, minChi2dof, minDeltaChi2dof, pkmax, Zpkmax, xcen, xcenDYE) = l.split(',')
                
                readData.setdefault('readPos',[]).append(int(qindex))
                readData.setdefault('readBase',[]).append(qbase)
                readData.setdefault('seqPos',[]).append(int(tindex))
                readData.setdefault('seqBase',[]).append(tbase)
                readData.setdefault('alignAnnotation',[]).append(alignAnnotation)
                readData.setdefault('alignOrientation',[]).append(orientation)
                if ReadTableType == AlignedSeqReadWithPulseFeature:
                    readData.setdefault('t1',[]).append(float(t1))
                    readData.setdefault('dt',[]).append(float(dt))
                    readData.setdefault('ts',[]).append(float(ts))
                    readData.setdefault('PNR',[]).append(float(PNR))
                    readData.setdefault('SNR',[]).append(float(SNR))
                    readData.setdefault('minChi2dof',[]).append(float(minChi2dof))
                    readData.setdefault('minDeltaChi2dof',[]).append(float(minDeltaChi2dof))
                    readData.setdefault('pkmax',[]).append(float(pkmax))
                    readData.setdefault('Zpkmax',[]).append(float(Zpkmax))
                    readData.setdefault('xcen',[]).append(float(xcen))
                    readData.setdefault('xcenDYE',[]).append(float(xcenDYE))
            read = ReadTableType(readMeta, readData)
            
            aStore.writeARead('testSeq', read)
    
    #aStore.sortAllSeqReadTableBySeqPos('testSeq')
    #aStore.sortAllSeqReadTableByRead('testSeq')
    indexFile.close()
    aStore.getMSA("testSeq", 500, 700)
    
    aStore.close()
    

if __name__ == '__main__':
    test()
