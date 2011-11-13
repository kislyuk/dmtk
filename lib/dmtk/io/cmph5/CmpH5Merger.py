#!/usr/bin/env python
import numpy as n
import os
import sys
import logging
import subprocess
import h5py
from re import findall,sub
import datetime

from dmtk.io.cmph5 import factory
from dmtk.io.VersionUtils import getSmrtAnalysisVersionInfo

# NOTES:
# (1) Consensus is not getting coppied over

class CmpH5Merger():
    """
    Class for merging cmp.h5 of similar specs into one file with or without
    preserving Pulse Feature information. Sorting is NOT preserved, should be
    faster to merge then sort than merge-sort (check!).
    """
    def __init__(self, filenames, outfile, forceMerge=False):
        self._outfile = outfile
        self._seedFN = filenames[0]
        self._FNs = filenames[1:]
        self._forceMerge = forceMerge
        self._setup()
    
    def run(self):
        """
        Merge cmp.h5 files in _FNs onto seeding cmp.h5 file _seedFN
        """              
        for fin in self._FNs:
            cmph5_in = h5py.File(fin,'r')
                        
            if self._validateCmpH5(cmph5_in, self._forceMerge):
                self.extendMovieInfo(cmph5_in)
                self.extendRefGroup(cmph5_in)
                self.extendAlnGroup(cmph5_in)
                self.extendAlnInfo(cmph5_in)
            else:
                logging.debug('[%s] is not valid for merging!' % cmph5_in.filename)
                # raise Exception('Invalid cmp.h5 for merging!')

        # Reset sorting (might not be necessary!)
        if self.cmph5_out['/AlnInfo/AlnIndex'].shape[0] != 0 and 'Index' in self.cmph5_out.attrs:
            self.cmph5_out['/AlnInfo/AlnIndex'][:,self._IDX['nBackRead']] = 0
            self.cmph5_out['/AlnInfo/AlnIndex'][:,self._IDX['nReadOverlap']] = 0
            self.cmph5_out.attrs.__delitem__('Index')
        
        self.cmph5_out.close()

        # Reset MoleculeID and set Log entry
        t_cmph5 = factory.create(self._outfile,'a')
        t_cmph5.resetMoleculeIDs()
        t_cmph5.log('CmpH5Merger.py',
                    getSmrtAnalysisVersionInfo()[0],
                    str(datetime.datetime.now().isoformat()),
                    ' '.join(sys.argv),
                    'Merging')
        t_cmph5.close()
        
    #################
    # Merge methods #
    #################
    
    def extendAlnInfo(self, cmph5_in):
        """
        Merge cmph5_in's /AlnInfo with the seeds.
        """       
        logging.info('Extending alninfo from [%s]' % cmph5_in.filename.split('/')[-1])

        aIdx_in = cmph5_in['/AlnInfo/AlnIndex'].value
        aIdx = self.cmph5_out['/AlnInfo/AlnIndex']
        aIdx_in[:,self._IDX['AlnID']] = aIdx_in[:,self._IDX['AlnID']] + aIdx.shape[0]
        aIdx_in[:,self._IDX['MovieID']] = n.array([self.F_MovieID[x] for x in aIdx_in[:,self._IDX['MovieID']]], dtype='uint32')
        aIdx_in[:,self._IDX['AlnGroupID']] = n.array([self.F_AlnGroupPath[x] for x in aIdx_in[:,self._IDX['AlnGroupID']]], dtype='uint32')
        aIdx_in[:,self._IDX['RefGroupID']] = n.array([self.F_RefGroupPath[x] for x in aIdx_in[:,self._IDX['RefGroupID']]], dtype='uint32')

        self._extendDset(aIdx, aIdx_in)
        
        for subgrp in [key for key in self.cmph5_out['/AlnInfo'].keys() if key != 'AlnIndex']:
            dout = self.cmph5_out['/AlnInfo'][subgrp]
            newVal = cmph5_in['/AlnInfo'][subgrp].value
            self._extendDset(dout, newVal)

        self.cmph5_out['/AlnInfo'].attrs.modify('nRow', self.cmph5_out['/AlnInfo/AlnIndex'].shape[0])

    def extendAlnGroup(self,cmph5_in):
        """
        Merge cmph5_in's /AlnGroup with the seeds.
        """
        logging.info('Extending alngroups from [%s]' % cmph5_in.filename.split('/')[-1])
        
        cache_AGP = cmph5_in['/AlnGroup/Path'].value.tolist()
        lastRGRPID = n.max(self.cmph5_out['/AlnGroup/ID'].value.tolist())
        rgpMap = dict(zip(cache_AGP,cmph5_in['/AlnGroup/ID'].value.tolist()))
        allRGRPs = [str(x) for x in cache_AGP]
        for rgp in allRGRPs:
            if not self.C_RGPToCopy.get(rgp):
                self.C_RGPToCopy[rgp] = rgp

        alngrp = [[None]*len(self.C_RGPToCopy) for i in xrange(3)] # newAlnGrpID, oldAlnGrpID, newAlnGrpPath
        cg = 0
        for subgrp in self.C_RGPToCopy.iteritems():
            if subgrp[0].endswith("NULL_GROUP") and subgrp[1] in self.cmph5_out:
                continue
            self.cmph5_out.copy(cmph5_in[subgrp[0]], subgrp[1]) # Most expensive operation!
            lastRGRPID += 1
            alngrp[0][cg] = lastRGRPID
            alngrp[1][cg] = rgpMap[subgrp[0]]
            alngrp[2][cg] = subgrp[1]
            cg += 1

        dout = self.cmph5_out['/AlnGroup/Path']
        newVal = n.array(alngrp[2][:cg], dtype=object)
        self._extendDset(dout, newVal)

        dout = self.cmph5_out['/AlnGroup/ID']
        newVal = n.array(alngrp[0][:cg], dtype='uint32')
        self._extendDset(dout, newVal)

        # Get mappings to fix AlignmentIndex
        for i in zip(alngrp[1][:cg],alngrp[0][:cg]):
            self.F_AlnGroupPath[i[0]] = i[1]

        self.cmph5_out['/AlnGroup'].attrs.modify('nRow', self.cmph5_out['/AlnGroup/ID'].shape[0])

    def extendRefGroup(self, cmph5_in):
        """
        Merge cmph5_in's /RefGroup with the seeds.
        """
        logging.info('Extending refgroups from [%s]' % cmph5_in.filename.split('/')[-1])

        t_refDict = self._getRefDict(cmph5_in)
        lastRefID = n.max(self.cmph5_out['/RefGroup/ID'].value)
        lastRefName = n.sort(self.cmph5_out['/RefGroup/Path'].value)[-1]       
        self.F_AlnGroupPath = {}
        self.C_RGPToCopy = {} # input rgrp: output rgrp
        changeMap = {} # inputCmph5 RefGroupPath: outputCmph5 RefGroupPath
        if self.refDict != t_refDict:
            for ref in t_refDict:
                exists = False
                for backref in self.refDict:
                    if self.refDict[backref] == t_refDict[ref]:
                        changeMap[ref] = backref
                        exists = True
                if not exists:
                    newRefName = ref 
                    if ref in self.refDict: 
                        newRefName = '/ref%06d' % (int(findall('ref(\\d+)',lastRefName)[0])+1)
                        lastRefName = newRefName
                    lastRefID += 1
                    
                    changeMap[ref] = newRefName
                    self.outIDDict[newRefName] = lastRefID
        else:
            t = self.cmph5_out['/RefGroup/Path'].value.tolist()
            changeMap = dict(zip(t,t))

        inIDDict = dict(zip(cmph5_in['/RefGroup/Path'].value,cmph5_in['/RefGroup/ID'].value))
        self.F_RefGroupPath = dict([(inIDDict[x[0]],self.outIDDict[x[1]]) for x in changeMap.items()])

        cache_RIFN = self.cmph5_out['/RefInfo/FullName'].value
        cache_RIID = self.cmph5_out['/RefInfo/ID'].value
        t_lastRefInfoID = n.max(self.cmph5_out['/RefGroup/RefInfoID'])
        
        refg = [[None]*len(changeMap) for i in xrange(3)] # ID, Path, RefInfoID
        refi = [[None]*len(changeMap) for i in xrange(4)] # ID, FullName, Length, MD5
        ci = 0
        cg = 0
        for oldRefPath in changeMap:
            newRefPath = changeMap[oldRefPath]
            if not self.refDict.get(newRefPath):                
                # RefGroupID
                self.refDict[newRefPath] = t_refDict[oldRefPath]
                refg[0][cg] = self.outIDDict[newRefPath]

                # RefGroupPath
                refg[1][cg] = newRefPath

                # RefInfo
                if self.refDict[newRefPath][0] in cache_RIFN:
                    idx = n.nonzero(cache_RIFN == self.refDict[newRefPath][0])[0]
                    lastRefInfoID = cache_RIID[idx][0]
                    refg[2][cg] = lastRefInfoID
                else:
                    t_lastRefInfoID = t_lastRefInfoID + 1
                    refg[2][cg] = t_lastRefInfoID

                    refi[0][ci] = t_lastRefInfoID
                    refi[1][ci] = self.refDict[newRefPath][0]
                    refi[2][ci] = self.refDict[newRefPath][1]
                    refi[3][ci] = self.refDict[newRefPath][2]
                    ci += 1
                cg += 1                

                # Add new RefGroup
                self.cmph5_out.create_group(newRefPath)
                for subgrp in cmph5_in[oldRefPath].values():
                    if not findall('Consensus',subgrp.name):
                        self.C_RGPToCopy[subgrp.name] = self._fixRefId(newRefPath[1:], subgrp.name)

            elif self.refDict.get(newRefPath) and newRefPath != oldRefPath:
                for subgrp in cmph5_in[oldRefPath].values():
                    if not findall('Consensus',subgrp.name):
                        self.C_RGPToCopy[subgrp.name] = self._fixRefId(newRefPath[1:], subgrp.name)

        if ci:
            for i,dset in enumerate([('ID','uint32'),('FullName',object),('Length','uint32'),('MD5',object)]):
                dout = self.cmph5_out['/RefInfo'][dset[0]]
                newVal = n.array(refi[i][:ci], dtype=dset[1])
                self._extendDset(dout, newVal)
        
        if cg:
            for i,dset in enumerate([('ID','uint32'),('Path',object),('RefInfoID','uint32')]):
                dout = self.cmph5_out['/RefGroup'][dset[0]]
                newVal = n.array(refg[i][:cg], dtype=dset[1])
                self._extendDset(dout, newVal)

        self.cmph5_out['/RefGroup'].attrs.modify('nRow', self.cmph5_out['/RefGroup/ID'].shape[0])
        self.cmph5_out['/RefInfo'].attrs.modify('nRow', self.cmph5_out['/RefInfo/ID'].shape[0])

    def extendMovieInfo(self, cmph5_in):
        """
        Merge cmph5_in's /MovieInfo with the seeds.
        Notes: Using .value and NOT .value.tolist() to avoid strange
        events when comparing string with numpy.object_ 
        """
        logging.info('Merging %d movies from [%s] into the existing %d' % 
                     (cmph5_in['/MovieInfo/Name'].shape[0],                      
                      cmph5_in.filename.split('/')[-1],
                      self.cmph5_out['/MovieInfo/Name'].shape[0]))

        # Extend MovieInfo ONLY with new movies
        cache_MIMin = cmph5_in['/MovieInfo/Name'].value
        cache_MIMout = self.cmph5_out['/MovieInfo/Name'].value
        newMovies = [m for m in cache_MIMin if m not in cache_MIMout]
        tin_movieDict = dict(zip(cache_MIMin, cmph5_in['/MovieInfo/ID'].value))
        tout_movieDict = dict(zip(cache_MIMout, self.cmph5_out['/MovieInfo/ID'].value))
        if newMovies:
            newidx = [cmph5_in['/MovieInfo/Name'].value == m for m in newMovies]
            newidx = reduce(lambda x,y: x|y,newidx)
            oldidx = n.logical_not(newidx)
            self.F_MovieID = [(tin_movieDict[m],tout_movieDict[m]) for m in cache_MIMin[oldidx]]
        
            for dsetName in ['Exp','Name','Run']:
                if dsetName in self.cmph5_out['/MovieInfo'].keys():
                    dout = self.cmph5_out['/MovieInfo'][dsetName]
                    din = cmph5_in['/MovieInfo'][dsetName]
                    self._extendDset(dout, din.value[newidx])

            # MovieID
            offsetID = n.max(self.cmph5_out['/MovieInfo/ID'].value) + 1
            dout = self.cmph5_out['/MovieInfo/ID']
            din = cmph5_in['/MovieInfo/ID'].value[newidx]
            newVal = n.array(range(offsetID, offsetID + din.shape[0]), dtype='uint32')
            self._extendDset(dout, newVal)

            self.F_MovieID.extend(zip(din,newVal))
            self.F_MovieID = dict(self.F_MovieID)

            # Fix MasterDataset nRow
            self.cmph5_out['/MovieInfo'].attrs.modify('nRow', self.cmph5_out['/MovieInfo/ID'].shape[0])
        else:
            self.F_MovieID = dict([(tin_movieDict[m],tout_movieDict[m]) for m in cache_MIMin])
        
    #################
    # Setup methods #
    #################

    def _setup(self):
        """
        Setup CmpH5Merger using the first input cmp.h5 as the yardstick
        for comparing the rest and making sure they meet the same standards
        """
        badSeed = True
        while badSeed:            
            self.cmph5_out = h5py.File(os.path.abspath(self._seedFN), 'r')
            self._valDict = self._getValDict(self.cmph5_out, self._forceMerge)

            # If it is empty and we have at least one more cmpH5, make that the seed and toss this one.
            if not self._valDict['AlnInfoDSets'] and len(self._FNs) > 0:
                logging.info('Rejecting [%s] as seed file (0-size AlignmentIndex or AlnGroups don\'t share the same datasets' % self._seedFN)
                self._seedFN = self._FNs[0]
                self._FNs = self._FNs[1:]
                self.cmph5_out.close()
            else:
                badSeed = False
                self.cmph5_out.close()
                logging.info('Copying seed cmp.h5 from [%s]' % self._seedFN)
                subprocess.call('cp %s %s' % (os.path.abspath(self._seedFN), self._outfile), shell=True)
                self.cmph5_out = h5py.File(self._outfile, 'r+')                            

        # Check if the seed is empty
        if not self._valDict['AlnInfoDSets'] and not len(self._FNs):
            logging.info('The seed [%s] has no alignments and there are no more files left to merge' % self._seedFN)

        else:
            # Get AlignmentIndex column indices
            self._IDX = factory.create(self._seedFN,'r')._colName2Index

            # Private methods for processing
            self._fixRefId = lambda newSt,oldSt: sub('ref\d+', newSt, oldSt)

            # Sanitize the seed
            self._sanitizeSeed()

            # Create reference dictionary
            self.refDict = self._getRefDict(self.cmph5_out)
            self.outIDDict = dict(zip(self.cmph5_out['/RefGroup/Path'].value,self.cmph5_out['/RefGroup/ID'].value))

            # Remove consensus 
            for rgrp in self.cmph5_out['/RefGroup/Path'].value.tolist():
                if 'Consensus' in self.cmph5_out[rgrp].keys():
                    del self.cmph5_out[rgrp+'/Consensus']

            # Remove sorting 
            if 'SortedRefIndex' in self.cmph5_out['/RefGroup'].keys():
                del self.cmph5_out['/RefGroup/SortedRefIndex']

    def _sanitizeSeed(self):
        """
        Reset all IDs found in /AlnGroup/ID and /MovieInfo/ID, to series
        of integers starting at 1 and ending at the number of rows of
        each dataset.
        """
        aIdx = self.cmph5_out['/AlnInfo/AlnIndex']
        
        # Fix MovieID and ReadGroupPathID
        dsetMap = {'/MovieInfo/ID':'MovieID','/AlnGroup/ID':'AlnGroupID'}
        for dsetName in dsetMap:
            idFIX = {}
            toFix = False
            for i,t_id in enumerate(self.cmph5_out[dsetName].value):
                idFIX[t_id] = i+1
                if t_id != i+1:
                    toFix = True
                    self.cmph5_out[dsetName][i] = i+1

            if idFIX and toFix:
                aIdx[:,self._IDX[dsetMap[dsetName]]] = n.array([idFIX[x] for x in aIdx[:,self._IDX[dsetMap[dsetName]]]], dtype='uint32')

    def _getValDict(self, cmph5, ignorePM=False):
        """
        Generate dictionary used for validating .cmp.h5 files on their
        own as well as against each other.
        """
        if cmph5['/AlnInfo/AlnIndex'].shape[0] == 0:
            return {'AlnInfoDSets':'','PulseMetrics':''}
        valDict = {}
        
        if not ignorePM:
            valDict['PulseMetrics'] = []
            a = [sorted(cmph5[grp].keys()) for grp in cmph5['/AlnGroup/Path'].value.tolist() if 'AlnArray' in cmph5[grp].keys()]
            valDict['PulseMetrics'] = a[0]
            if filter(lambda x: x != a[0], a):
                return {'AlnInfoDSets':'','PulseMetrics':''}

        valDict['AlnInfoDSets'] = sorted(cmph5['/AlnInfo'].keys())
        
        return valDict

    ##################
    # Helper Methods #
    ##################

    def _validateCmpH5(self, cmph5, ignorePM=False):
        """
        Validate cmph5 against the seed.
        """
        t_valDict = self._getValDict(cmph5, ignorePM)
        if not t_valDict['AlnInfoDSets']:
            logging.info('Rejecting [%s] since it has a 0-sized AlnIndex' % cmph5.filename)
        elif sorted(t_valDict['AlnInfoDSets']) != sorted(self._valDict['AlnInfoDSets']):
            t_mdsets = ','.join([d for d in self._valDict['AlnInfoDSets'] if d not in t_valDict['AlnInfoDSets']])
            logging.info('Rejecting [%s] since it is missing AlnInfo DataSets: %s' % (cmph5.filename, t_mdsets))
        if not ignorePM:
            if not t_valDict['PulseMetrics']:
                logging.info('Rejecting [%s] since AlnGroup datasets are not shared by all AlnGroups' % cmph5.filename)
            elif sorted(t_valDict['PulseMetrics']) != sorted(self._valDict['PulseMetrics']):
                t_mdsets = ','.join([d for d in self._valDict['PulseMetrics'] if d not in t_valDict['PulseMetrics']])
                logging.info('Rejecting [%s] since it is missing AlnGroup DataSets: %s' % (cmph5.filename, t_mdsets))

        return self._getValDict(cmph5, ignorePM) == self._valDict

    def _getRefDict(self,cmph5):
        """
        Encapsulate all data from /RefGroup/ into a dictionary. 
        Notes: Caching into lists necessary for speed and to avoid 
        having numpy.object_ isntances instead of strings.
        """
        rInfo = dict(zip(cmph5['/RefInfo/ID'].value.tolist(),
                         zip(cmph5['/RefInfo/FullName'].value.tolist(),
                             cmph5['/RefInfo/Length'].value.tolist(),
                             cmph5['/RefInfo/MD5'].value.tolist())))
        cache_RGP = cmph5['/RefGroup/Path'].value.tolist()
        cache_RGRID = cmph5['/RefGroup/RefInfoID'].value.tolist()
        return dict([(cache_RGP[i],rInfo[rID])for i,rID in enumerate(cache_RGRID)])

    def _extendDset(self, dout, newVal):
        """
        Extend the dout HDF5 Dataset with newVal. Set the 'lastRow'
        attribute where appropriate.
        """
        newVLen = newVal.shape[0]
        oldVlen = dout.shape[0]

        if len(dout.shape) == 1:
            newShape = (dout.shape[0] + newVLen,)
            dout.resize(newShape)
            dout[oldVlen:(oldVlen + newVLen)] = newVal
        else:
            newShape = (dout.shape[0] + newVLen, dout.shape[1])
            dout.resize(newShape)
            dout[oldVlen:(oldVlen + newVLen),:] = newVal                

        if 'lastRow' in dout.attrs.keys():
            dout.attrs.modify('lastRow', oldVlen+newVLen)
       
if __name__ == '__main__':
    import glob
    from re import findall
    import subprocess
    
    logLevel = logging.DEBUG
    logFormat = '%(asctime)s [%(levelname)s] %(message)s'
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)

    FNs = []
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_1/?.cmp.h5'))
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_2/ecoli?.cmp.h5'))
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_3/?.cmp.h5'))
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_4/?.cmp.h5'))
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_5/?.cmp.h5'))
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_6/*.cmp.h5'))
    FNs.append(sorted(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_7/control_reads*.cmp.h5'))[::-1])
    FNs.append(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/test_8/?.cmp.h5'))
    
    VERs = sorted(glob.glob('/home/UNIXHOME/diliopoulos/Projects/Software/regTests/CmpH5Merger/equivTest/mergedOld_*'))
    blist = [8]

    # BYPASS:
    # FNs.append(sorted(glob.glob('/scratch/mergeTest/test*.h5')))
    # FNs.append(sorted(glob.glob('/scratch/mergeTest/aligned*.h5')))
    # FNs.append(sorted(glob.glob('/home/UNIXHOME/diliopoulos/Bucket/aligned*.h5')))
    # from IPython.Shell import IPShellEmbed; IPShellEmbed(argv=[])()

    for i,fn in enumerate(FNs):
        print 'Running Test %d of %d...' % (i+1, len(FNs))
        outFN = './mergedNew_test_%d.cmp.h5' % (i+1)
        c = CmpH5Merger(fn, outFN)
        c.run()
        
        runCode1 = 0
        if (i+1) not in blist:
            cmd = 'cmpH5Compare.py ./%s %s' % (outFN,VERs[i])
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = p.communicate()
            retCode1 = p.returncode
            logging.info('cmpH5Compare STDOUT: %s [%s]' % (out.strip(), 'PASS' if not retCode1 else 'FAIL'))

        cmd = 'cmpH5Sort.py ./%s' % outFN
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        retCode2 = p.returncode
        logging.info('cmpH5Sort STDOUT: %s [%s]' % (out.strip(), 'PASS' if not retCode2 else 'FAIL'))

        print 'Finished test %d of %d [%s]' % (i+1, len(FNs), 'PASS' if retCode1 == retCode2 == 0 else 'FAIL')
        subprocess.call('rm ./%s' % outFN, shell=True)

