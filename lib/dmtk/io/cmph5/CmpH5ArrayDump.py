from re import findall
import math
import numpy as n

Z_A = -77.27
Z_C = 0.08540
Z_S = 0.001210
FIXED_LENGTH = 500
BASES = ['t','g','a','c']
QV_THRESHOLDS = [6,8,10]

class ArrayReporter(object):
    """Generates numpy record arrays from .cmp.h5"""

    def __init__(self, the_cmph5):
        self._cmph5 = the_cmph5

    def get_perReadAlnInfo(self):
        columns = [('zmwid', int),
                   ('runcode', n.dtype('S80',1)),
                   ('movie', n.dtype('S100',1)),
                   ('subreadId', n.int32),
                   ('readStrand', n.dtype('S1',1)),
                   ('readStart', n.int32),
                   ('readEnd', n.int32),
                   ('targetStrand', n.dtype('S1',1)),
                   ('targetStart', n.int32),
                   ('targetEnd', n.int32),
                   ('zScore', n.int32),
                   ('nInsert', n.int32),
                   ('nDelete', n.int32),
                   ('nMismatch', n.int32),
                   ('nCorrect', n.int32),
                   ('template', n.dtype('S80',1))]
        data = n.recarray((self._cmph5.numAlnHits,),dtype=columns)

        aIdx = self._cmph5["/AlnInfo"].asRecArray()
        
        data['zmwid'] = aIdx['HoleNumber']
        tempRunMovie = n.array([('%07d-%04d' % (self._cmph5.getExp(x), self._cmph5.getRun(x)),self._cmph5['MovieName'][x]) for x in aIdx['MovieId']])
        data['runcode'] = tempRunMovie[:,0]
        data['movie'] = tempRunMovie[:,1]
        data['subreadId'] = aIdx['SubreadId']
        data['readStrand'] = '+'
        data['targetStrand'] = n.array(['+' if x == 0 else '-' for x in aIdx['RCRefStrand']])
        data['readStart'] = aIdx['rStart']
        data['readEnd'] = aIdx['rEnd']
        data['targetStart'] = aIdx['tStart']
        data['targetEnd'] = aIdx['tEnd']
        data['zScore'] = self._cmph5['ZScore']
        data['nInsert'] = aIdx['nIns']
        data['nDelete'] = aIdx['nDel']
        data['nMismatch'] = aIdx['nMM']
        data['nCorrect'] = data['readEnd']-data['readStart']-data['nInsert']-data['nDelete']-data['nMismatch']
        refDict = dict([(self._cmph5['/RefGroup'].getID(Path=x[1]),x[0]) for x in self._cmph5._fullName2RefGroup.iteritems()])
        data['template'] = n.array(['"%s"' % (refDict[x]) for x in aIdx['RefSeqId']])

        return data
         
    def get_perMoviePMInfo(self):
        """Get per movie Pulse Metric info (per channel) for
        each movie in the .cmp.h5 file"""
        m_idx = self._cmph5.indexCols.index("MovieId")
        movies = dict([(i,[]) for i in self._cmph5["MovieName"].values()])
        moviemap = dict(zip(list(self._cmph5['MovieName']),list(self._cmph5['MovieID'])))

        columns = [('runcode', n.dtype('S80',1)),
                   ('movie', n.dtype('S80',1)),
                   ('channel', n.int32),
                   ('pulseWidth', float),
                   ('pkmid', float),
                   ('pulseRate', float)]

        data = n.core.recarray((len(movies.keys())*4,),dtype=columns)
        data.fill('0')

        # Check if the .cmp.h5 has pulse info loaded
        try:
            almtest = self._cmph5.alnHitIterator().next()
            if not almtest.hasPulseInfo:
                return data
        except StopIteration:
            return data

        for refseq in self._cmph5.refGroupIterator():
            for rdgrp in refseq.readGroupIterator():
                rdgrp_movies = rdgrp.alnHitsByMovie()

                for movieid,allalnhits in rdgrp_movies.iteritems():
                    for alnhit in allalnhits:
                        row = self._get_alnhitPM(alnhit)
                        if row != []:
                            movies[movieid].append(row)

        for movidx,movieStats in enumerate(movies.iteritems()):
            alnhits = n.array(movieStats[1])
            for chid,idx in enumerate(range(12)[::3]):
                mstridx = (movidx*4)+chid
                movieid = moviemap[movieStats[0]]
                data[mstridx]['runcode'] = '%d-%04d' % (self._cmph5["Experiment"][movieid],self._cmph5["Run"][movieid])
                data[mstridx]['movie'] = movieStats[0]
                data[mstridx]['channel'] = chid
                data[mstridx]['pulseWidth'] = n.median(alnhits[:,idx])   if alnhits != [] else 0.0
                data[mstridx]['pkmid']      = n.median(alnhits[:,idx+1]) if alnhits != [] else 0.0 
                data[mstridx]['pulseRate']  = n.median(alnhits[:,idx+2]) if alnhits != [] else 0.0

        return data

    def get_perMovieAlnInfo(self):
        """Get per movie alignment statisitcs for each movie in
        the .cmp.h5 file"""
        m_idx = self._cmph5._colName2Index["MovieId"]
        qs_idx = self._cmph5._colName2Index["rStart"]
        qe_idx = self._cmph5._colName2Index["rEnd"]
        nI_idx = self._cmph5._colName2Index["nIns"]
        nD_idx = self._cmph5._colName2Index["nDel"]
        nMM_idx = self._cmph5._colName2Index["nMM"]
        ref_idx = self._cmph5._colName2Index["RefSeqId"]
        rdg_idx = self._cmph5._colName2Index["ReadGroupId"]
        ob_idx = self._cmph5._colName2Index["offset_begin"]
        oe_idx = self._cmph5._colName2Index["offset_end"]      

        refdict = {}
        for idx,ref in enumerate(self._cmph5.refInfoIterator()):
            refdict[idx+1] = (ref.name, ref.fullName)
        
        almIdx = self._cmph5['AlignmentIndex']
        zscore = self._cmph5['ZScore']
        movNumb = len(self._cmph5['MovieName'].values())
        columns = [('movie', n.dtype('S80',1)),
                   ('runcode', n.dtype('S80',1)),
                   ('nReads', n.int32),
                   ('medianZ', float),
                   ('medianAcc', float),
                   ('medianLength', float),
                   ('nReadsAboveZ3', float),
                   ('medianZaboveZ3', float),
                   ('AZ50', float),
                   ('Z95', float),
                   ('AZ95', float),
                   ('medianAccAboveZ3', float),
                   ('medianLengthAboveZ3', float),
                   ('nQVs>=6', n.int32),
                   ('nQVs>=8', n.int32),
                   ('nQVs>=10', n.int32),
                   ('refseq', n.dtype('S80',1))]
        data = n.core.recarray((movNumb*len(refdict),),dtype=columns)
        data.fill('0')

        for movidx in xrange(movNumb):
            for refidx in refdict.keys():               
                t_movidx = movidx*len(refdict) + refidx - 1
                data[t_movidx]['movie'] = self._cmph5['MovieName'][movidx+1]
                data[t_movidx]['runcode'] = '%d-%04d' % (self._cmph5["Experiment"][movidx+1],self._cmph5["Run"][movidx+1])
                data[t_movidx]['refseq'] = '"%s"' % refdict[refidx][1] 
                sliceidx = n.nonzero((almIdx[:,m_idx] == movidx+1) & (almIdx[:,ref_idx] == refidx))[0]
                if len(sliceidx) > 0:
                    almrows = almIdx[sliceidx,:]                              
                    if len(sliceidx) == 1:
                        almrows = almrows.reshape((1,len(almrows))) 

                    zrows = zscore.asNumPy[sliceidx]
                    
                    data[t_movidx]['nReads'] = almrows.shape[0]
                    data[t_movidx]['medianZ'] = n.median(zrows)

                    almlength = almrows[:,qe_idx] - almrows[:,qs_idx]
                    nError = n.array(almrows[:,nI_idx] + almrows[:,nD_idx] + almrows[:,nMM_idx],dtype=float)
                    data[t_movidx]['medianAcc'] = n.median(1 - (nError/almlength))*100
                    data[t_movidx]['medianLength'] = n.median(almlength)

                    z3idx = n.nonzero(zrows>3)[0]
                    almrows = almrows[z3idx,:]
                    zrows = zrows[z3idx,:]
                    data[t_movidx]['nReadsAboveZ3'] = almrows.shape[0]
                    data[t_movidx]['medianZaboveZ3'] = n.median(zrows)
                    data[t_movidx]['AZ50'] = 100.0 * self._accFromZ(data[t_movidx]['medianZaboveZ3'])
                    data[t_movidx]['Z95'] = self._quantile(list(zrows),0.95)        
                    data[t_movidx]['AZ95'] = 100.0 * self._accFromZ(data[t_movidx]['Z95'])

                    almlength = almrows[:,qe_idx] - almrows[:,qs_idx]
                    nError = n.array(almrows[:,nI_idx] + almrows[:,nD_idx] + almrows[:,nMM_idx],dtype=float)
                    data[t_movidx]['medianAccAboveZ3'] = n.median(1 - (nError/almlength))*100
                    data[t_movidx]['medianLengthAboveZ3'] = n.median(almlength)

                    # Calculate QVs
                    for rgpId in n.unique(almrows[:,rdg_idx]):
                        pairs = zip(almrows[almrows[:,rdg_idx] == rgpId,ob_idx],almrows[almrows[:,rdg_idx] == rgpId,oe_idx])
                        QVs = self._cmph5.readGroupById(rgpId)['QualityValue'].asNumPy / 10
                        
                        for qvt in QV_THRESHOLDS:
                            qvs = sum([len(x[x>=qvt]) for x in [QVs[m[0]:m[1]] for m in pairs]])
                            data[t_movidx]['nQVs>=%s' % str(qvt)] = data[t_movidx]['nQVs>=%s' % str(qvt)] + qvs
                            
        return data
        

    ######################
    ## HELPER FUNCTIONS ##
    ######################
    def _get_alnhitInfo(self,alnhit,columns):
        """Return a record array containing alignment info for
        a single read"""
        l = n.core.records.fromarrays(n.zeros(shape=len(columns)),dtype=columns)
        l['zmwid'] = alnhit.alignmentIndexData["HoleNumber"]
        l['runcode'] = alnhit.experiment + "-" + alnhit.run
        l['movie'] = alnhit.query_id.split('/')[0]
        l['subreadId'] = alnhit.alignmentIndexData['SubreadId']
        l['readStrand'] = alnhit.query_strand
        l['readStart'] = alnhit.query_start
        l['readEnd'] = alnhit.query_end
        l['targetStrand'] = alnhit.target_strand
        l['targetStart'] = alnhit.target_start
        l['targetEnd'] = alnhit.target_end
        l['zScore'] = alnhit.zScore
        l['nInsert'] = alnhit.nIns
        l['nDelete'] = alnhit.nDel
        l['nMismatch'] = alnhit.nMismatch
        l['nCorrect'] = alnhit.query_end-alnhit.query_start-alnhit.nIns-alnhit.nDel-alnhit.nMismatch
        l['template'] = '"%s"' % alnhit.target_id

        return l
      
    def _get_alnhitPM(self,alnhit):
        """Return a list containing per read, per channel pulse
        metrics"""
        l = [] # [CH_#_medPW, CH_#_pkmid, rate]
        qseq = alnhit.alignedQuery.lower()
        pd = alnhit.pulseInfo

        rate = (alnhit.target_end-alnhit.target_start)/(pd['StartTime'][-1] -pd['StartTime'][0])
        
        for chbase in BASES:
            idx = n.nonzero(n.array(list(alnhit.alignedQuery)) == chbase)[0]
            if len(idx) == 0:
                return [] # Not sure about throwing away reads like this
            else:
                l.append(n.median(pd['PulseWidth'][idx]))
                l.append(n.median(pd['pkmid'][idx]))
                l.append(rate)
            
        return l
    
    def _quantile(self,v, quantile):
        """Calculate quantile given list of values"""
        if len(v)==0:
            return 0.0
        n = len(v)
        nq = int(round(quantile * float(n)))
        if nq > n-1: nq = n-1
        v.sort()
        return v[nq]

    def _accFromZ(self,z):
        """Calculate accuracy from Z score"""
        a = Z_S * z * math.sqrt( Z_A*Z_A + 1.0 )
        y = math.exp( a - Z_C - Z_A / ( FIXED_LENGTH+20.0 ) )
        acc = y / (1.0+y)
        return acc



    
