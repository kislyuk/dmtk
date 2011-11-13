#!/usr/bin/env python
import numpy as n
import os
import logging
import h5py

from dmtk.io.cmph5 import CmpH5Factory

###########
# Notes
# [1] For triming datasets in place, can't use visit() or visititems() since
#     HDF5.DataSets are protected and can't be resized when called via these
#     2 methods. Can't find a more elegant workaround than _trimDataSet
# [2] self._IDX is hardcoded due to segfault when loading CmpH5Factory

class CmpH5Splitter():

    def __init__(self, filename, outDir='', fullRefName=False):
        self._fullRefName = fullRefName
        self._seedFN = os.path.abspath(filename)
        self._outDir = os.getcwd() if not outDir else outDir
        self.cmph5_in = h5py.File(self._seedFN, 'r')
        self.refgrpDict = dict(zip(self.cmph5_in['/RefGroup/Path'],zip(self.cmph5_in['/RefGroup/ID'],self.cmph5_in['/RefGroup/RefInfoID'])))
        self.refinfoDict = dict(zip(self.cmph5_in['/RefInfo/ID'], self.cmph5_in['/RefInfo/FullName']))
        self.aIdx = self.cmph5_in['/AlnInfo/AlnIndex'].value
        
        self._IDX = {'AlnGroupID': 1, 'MovieID': 2, 'RefGroupID': 3}
        # self._IDX = CmpH5Factory.factory.create(self._seedFN,'r')._colName2Index
        
    def run(self):
        logging.info('Splitting [%s] into %d files...' % (os.path.basename(self._seedFN),len(self.refgrpDict)))
        outFNs = []
        for refgrp in self.refgrpDict:
            if self._fullRefName:
                outFN = self._outDir+'/%s.h5' % self.refinfoDict[self.refgrpDict[refgrp][1]]
            else:
                outFN = self._outDir+'/%s.cmp.h5' % refgrp[1:]
            outFNs.append(outFN)

            logging.info('Generating [%s]...' % os.path.basename(outFN))
            cmph5_out = h5py.File(outFN, 'w')

            # Copy over all essential HDF5 Groups
            for dset in ['/AlnGroup','/AlnInfo','/MovieInfo','/RefGroup','/RefInfo', '/FileLog']:
                cmph5_out.copy(self.cmph5_in[dset], dset)

            # Copy over root attributes
            map(lambda x: cmph5_out.attrs.create(x[0],x[1]), self.cmph5_in.attrs.items())

            # Trim RefGroup and RefInfo
            self._trimGroup(cmph5_out['/RefGroup'], cmph5_out['/RefGroup/Path'].value == refgrp)
            self._trimGroup(cmph5_out['/RefInfo'], cmph5_out['/RefInfo/ID'].value == self.refgrpDict[refgrp][1])

            # Trim MovieInfo and AlnGroup
            self._trimComplexGroup(cmph5_out['/MovieInfo'], 'MovieID', self.refgrpDict[refgrp][0])
            self._trimComplexGroup(cmph5_out['/AlnGroup'], 'AlnGroupID', self.refgrpDict[refgrp][0])

            # Trim AlnInfo
            self._trimGroup(cmph5_out['/AlnInfo'], self.aIdx[:,self._IDX['RefGroupID']] == self.refgrpDict[refgrp][0])

            # Copy RefGroup
            cmph5_out.copy(self.cmph5_in[refgrp], refgrp)
            cmph5_out.close()

        return outFNs

    ##################
    # Helper Methods #
    ##################
    def _trimComplexGroup(self, rgrp, idxName, refgrpID):
        rgrpName = rgrp.name
        ids = n.unique(self.aIdx[self.aIdx[:,self._IDX['RefGroupID']] == refgrpID, self._IDX[idxName]])
        idx = [self.cmph5_in[rgrpName]['ID'].value == t_id for t_id in self.cmph5_in[rgrpName]['ID'] if t_id in ids]
        idx = reduce(lambda x,y: x|y,idx)
        self._trimGroup(rgrp, idx)
    
    def _trimGroup(self, rgrp, idx):
        map(lambda x: self._trimDataSet(rgrp[x],idx), rgrp.keys())
        if 'nRow' in rgrp.attrs.keys():
            rgrp.attrs.modify('nRow', rgrp[rgrp.attrs['MasterDataset']].shape[0])
            
    def _trimDataSet(self, dout, idx):
        if len(dout.shape) == 1:
            newval = dout[idx]
        else:
            newval = dout[idx,:]
        dout.resize(newval.shape)
        dout.write_direct(newval)
        if 'lastRow' in dout.attrs.keys():
            dout.attrs.modify('lastRow', dout.shape[0])    
        
if __name__ == '__main__':   
    FN = '/home/UNIXHOME/jbullard/projects/software/R/common/122615/data/aligned_reads.cmp.h5'
    print '\n'.join(CmpH5Splitter(FN, fullRefName=True).run())
