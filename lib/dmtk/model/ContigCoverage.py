__doc__="""Encapsulates coverage information for per-contig coverage plots"""
import os
import numpy as np


MAX_BIN = "maxBin"
MEANS = "means"
MEAN_DEPTH_COV = "meanDepthOfCov"
AVE_REGION_SIZE = "aveRegionSize"
PERCENT_MISSING_BASES = "percentMissingBases"

class ContigCoverage:
    
    def __init__(self, seqId, name=None):
        """Encapsulates sequence info relevant to one chart"""
        self.seqid = seqId
        if name == None:
            name = seqId
        
        self.name = name
        self.xData = []
        
        self._numRecords = 0;
#        self._cumulativeNormMeans = 0

        self._totalCoverage = 0
        self._refEnd = 0
        self._refStart = None
        
        self._cumulativeRegionSizes = 0
        self._numBases = 0
        self._missingBases = 0
        self._windowSize = 0
        
        self.yDataMean = []
        self.yDataStdevPlus = []
        self.yDataStdevMinus = []
        self.fileName = "coveragePlot_%s%s"%(self.seqid,".png")

        
    def _toDictionary(self, attributeString):
        """convert an attribute line to a dictionary"""
        return dict(item.split("=") for item in attributeString.split(";"))
    
    def addData(self, gff3Record):
        """Append x,y data from this record to the contig graph"""
        
        self._numRecords = self._numRecords + 1
        
        if self._refStart == None:
            self._refStart = gff3Record.start
        
        self.xData.append(gff3Record.start);
        dict = self._toDictionary(gff3Record._getAttributeString())
        stats = dict['cov2'].split(",")
        mean = float(stats [0])
        stddev = float(stats [1])
        
        self.yDataMean.append( mean )
        self.yDataStdevPlus.append( mean+stddev )
        
        #what is better - recalculating each time data is added, or storing 2 more arrays?
        regSize = (gff3Record.end - gff3Record.start) + 1

        self._totalCoverage +=  mean * regSize
        self._refEnd = gff3Record.end
        
        self._cumulativeRegionSizes = self._cumulativeRegionSizes + regSize
        
        #the second value of gaps pair is missing bases for region
        self._missingBases = self._missingBases + int(dict['gaps'].split(",")[1]) 
        
        #assumption: regions are continuous
        if self._numBases < gff3Record.end:
            self._numBases = gff3Record.end
                
        #clip at zero
        lowerBound = mean-stddev
        if lowerBound < 0:
            lowerBound = 0
        self.yDataStdevMinus.append( lowerBound  )
        
        
        
    def meanCoveragePerBase(self):
        """Get the normalized coverage per base"""
#        return float (self._cumulativeNormMeans / self._numRecords)
        return self._totalCoverage / float (self._refEnd-self._refStart + 1)

    
    def aveRegionSize(self):
        """Get the average chunk size of this contig"""
        if self._numRecords == 0:
            return 0
        return float (self._cumulativeRegionSizes / self._numRecords)
    
    def missingBases(self):
        """Get number missing bases"""
        return self._missingBases
    
    def numBases(self):
        """Get number bases"""
        return self._numBases
        
        
        
def getReferenceCoverageStats(contigList):
    """Get a dictionary of coverage stats for the list of contigs"""
    statsDict = {}

    maxbin = 0
    means = []
    cumulativeAveRegionSize = 0
    totalNumBases = 0
    totalMissingBases = 0
    
    contigCoverages = []
    numContigs = len(contigList)
    for cc in contigList:
      
        contigCoverages.append( cc.meanCoveragePerBase() * cc.numBases() )
        cumulativeAveRegionSize += cc.aveRegionSize()
        
        if len( cc.yDataMean ) > 0:
            means.extend( cc.yDataMean )
            localmax = max( cc.yDataMean )
            if localmax > maxbin:
                maxbin = localmax
            
        totalNumBases += cc.numBases()
        totalMissingBases += cc.missingBases()
            
    statsDict[MAX_BIN] = maxbin
    statsDict[MEANS] = means

    aveRegionSize = float (cumulativeAveRegionSize / numContigs)    
    statsDict[MEAN_DEPTH_COV] = sum( contigCoverages ) / totalNumBases
    
    statsDict[AVE_REGION_SIZE] = int(aveRegionSize)
    statsDict[PERCENT_MISSING_BASES] = ( float(totalMissingBases)/float(totalNumBases) ) * 100
    return statsDict

  