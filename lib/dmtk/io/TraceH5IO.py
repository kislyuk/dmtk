## provide read only access to fetch data from trc.h5


import h5py
import numpy

class TrcH5File(object):
    def __init__(self, fn):
        self._h5handle = h5py.File(fn, 'r')
        f = self._h5handle # just a local alias
        self.basemap = dict(zip( f["/ScanData/DyeSet"].attrs["BaseMap"][0], (0,1,2,3) ))
        self.channle2base = dict(zip( (0,1,2,3), f["/ScanData/DyeSet"].attrs["BaseMap"][0]))
        try:
            self.frameRate = f["/ScanData/AcqParams"].attrs["FrameRate"][0]
        except TypeError:
            self.frameRate = f["/ScanData/AcqParams"].attrs["FrameRate"]
        print self.frameRate
        self._decodeMap = f["/TraceData/Codec/Decode"]
        self._traceData = f["/TraceData/Traces"]
        self._holeXY = f ["/TraceData/HoleXY"]
        self._holeNumber = f ["/TraceData/HoleNumber"][:]
        self._xy2holeNumber = dict( zip(  [(c[0], c[1]) for c in self._holeXY], self._holeNumber ) )
        self._holeNumber2index = dict( zip( self._holeNumber , range(len(self._holeNumber ) ) ) )

    def getTraceDataForHoleAtXY(self, x, y, startTime = None, endTime = None):
        hn = self._xy2holeNumber[ (x,y) ] 
        return self.getTraceDataForHole(hn, startTime, endTime)

    def getTraceDataForHole(self, holeNum, startTime = None, endTime = None):  
        """ if startTime and endTime are "None", fetch whole trace """
        print holeNum, self._xy2holeNumber[ (25,66) ]
        hIdx = self._holeNumber2index[holeNum]
        trData4Hole = self._traceData[hIdx][:]
        basemap = self.basemap
        decodeMap = self._decodeMap[:]
        (nC, nF) = trData4Hole.shape
        if startTime == None:
            startFrame = 0
        else:
            startFrame = int(startTime * self.frameRate)
        if endTime == None:
            endFrame =nF
        else:
            endFrame = int( endTime * self.frameRate )
            if endFrame >= nF:
                endFrame = nF -1 
               
        rtn = [ decodeMap.take(trData4Hole[i][startFrame:endFrame]) for i in xrange(len(basemap)) ]
        #rtn = [ decodeMap.take(trData4Hole[startFrame:endFrame][i]) for i in xrange(len(basemap)) ]
        return rtn, startFrame, endFrame, self.frameRate

def simpltTest():
    #trh5 = TrcH5File("./m100112_193152_Twe_p2_b15.trc.h5")
    #trh5 = TrcH5File("./m091125_165109_Cog_p2_b30.trc.h5")
    trh5 = TrcH5File("./m100810_153248_Orb_p1_b20.trc.h5")
    print trh5.getTraceDataForHole(2294, startTime = 1570, endTime = 1680)
    #print trh5.getTraceDataForHole(30)
    #print trh5.getTraceDataForHoleAtXY(30,28)

if __name__ == "__main__":
    simpltTest()
