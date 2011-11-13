#!/usr/bin/env python

import os
import sys
import tempfile
import shutil
import datetime
import h5py as H5
from numpy import *

from optparse import OptionParser
from dmtk.io.cmph5.CmpH5SortingTools import *

from dmtk.io.cmph5 import CmpH5Factory

__VERSION__ = ".64"

def __pathExists(h5, path):
    try:
        h5[path]
        return True
    except Exception, E:
        return False

def __breakInBlocks(size, nProc):
    bSize = [size / nProc] * nProc
    bSize[0] = (bSize[0] + (size - sum(bSize)))
    starts = concatenate(([0], cumsum(bSize)[:-1]))
    ends = cumsum(bSize)
    return(zip(starts, ends))

def __repackDataArrays(cH5, format, log):
    """
    Flatten read groups according to an indexed cmp.h5 file.
    """
    SORTED = "Sorted"

    alnGroups = [x for x in cH5[format.ALN_GROUP_PATH]]
    pulseDatasets = [cH5[x].keys() for x in alnGroups]
    uPulseDatasets = reduce(lambda x,y: set.union(set(x), set(y)), pulseDatasets)
    if (not all(map(lambda x : set(x) == uPulseDatasets, pulseDatasets))):
        log.error("All alignment groups need to have the same datasets.")
        raise Exception("Can only repack cmp.h5 files with consistent datasets across alignment groups.")

    readGroupPaths = dict(zip(cH5[format.ALN_GROUP_ID], [ x for x in cH5[format.ALN_GROUP_PATH]]))
    refGroupPaths = dict(zip(cH5[format.REF_GROUP_ID], [ x for x in cH5[format.REF_GROUP_PATH]]))
    uPDAndType = dict(zip(uPulseDatasets, [ cH5[readGroupPaths.values()[0]][z].dtype for z in uPulseDatasets ]))

    def getDataset(read, ds):
        return(cH5[readGroupPaths[read[format.ALN_ID]]][ds])

    def getRefGroup(gID):
        return(cH5[refGroupPaths[gID]])

    offsets = cH5[format.REF_OFFSET_TABLE].value
    sAI = cH5[format.ALN_INDEX]
    orderedRefPaths = [""] * offsets.shape[0]

    for row in xrange(0, offsets.shape[0]):
        log.msg("Processing reference group: %d of %d" % (row + 1, offsets.shape[0]))
        orderedRefPaths[row] = "/".join([getRefGroup(offsets[row, 0]).name, SORTED])

        fRow  = int(offsets[row, 1])
        lRow  = int(offsets[row, 2])
        
        ## Don't really have to do anything if there are no references
        ## which aligned.
        if (lRow == fRow):
            continue 

        ## Make a new Group.
        newGroup = getRefGroup(offsets[row, 0]).create_group(SORTED)
        log.msg("Created new read group: %s" % SORTED)

        ## Go through each read and write it into the new vector.
        reads = sAI[fRow:lRow, ]
        totalSizes = reads[:, format.OFFSET_END] - reads[:, format.OFFSET_BEGIN]
        for pulseDataset in uPulseDatasets: 
            log.msg("Processing dataset: %s" % pulseDataset)
            newDS = array([0]*sum(1 + totalSizes), dtype = uPDAndType[pulseDataset])
            currentStart = 0
            for readIdx in xrange(0, reads.shape[0]):
                read = reads[readIdx, ]
                gStart, gEnd = currentStart, currentStart + totalSizes[readIdx]
                newDS[gStart:gEnd] = getDataset(read, pulseDataset)[read[format.OFFSET_BEGIN]:read[format.OFFSET_END]]
                currentStart = gEnd + 1
            newGroup.create_dataset(pulseDataset, data = newDS, dtype = uPDAndType[pulseDataset], maxshape = None)
           
        ## After we've moved all of the data we can move the offsets.
        currentStart = 0
        for i in xrange(0, reads.shape[0]):
             reads[i, format.OFFSET_BEGIN] = currentStart
             reads[i, format.OFFSET_END] = currentStart + totalSizes[i]
             reads[i, format.ALN_ID] = row 
             currentStart = reads[i, format.OFFSET_END] + 1
        sAI[fRow:lRow,] = reads

    
    ## Now remake the AlnGroup Dataset.
    log.msg("Writing new AlnGroupPath values.")
    del(cH5[format.ALN_GROUP_PATH])
    del(cH5[format.ALN_GROUP_ID])
    cH5.create_dataset(format.ALN_GROUP_PATH, data = orderedRefPaths, 
                       dtype = H5.new_vlen(str), maxshape = None)
    cH5.create_dataset(format.ALN_GROUP_ID, data = range(0, offsets.shape[0]),
                       dtype = "int32", maxshape = None)
    for rg in readGroupPaths.values():
        del(cH5[rg])
    

def sortCmpH5(inFile, outFile, deep, jobs, log):
    """ 
    This routine takes a cmp.h5 file and sorts the AlignmentIndex
    table adding two additional columns for fast access. In addition,
    a new top-level attribute is added to the indicate that the file
    has been sorted, as well as a table to indicate the blocks of the
    alignment index associated with each reference group.
    """
    success = False;
    
    if (outFile):
        log.msg("Copying: " + inFile + " to " + outFile)
        shutil.copyfile(inFile, outFile)
        inFile = outFile

    try:
        cH5 = H5.File(inFile, 'a')
        format = CmpH5Format(cH5)
        log.msg("Read cmp.h5 with version %s" % format.VERSION)

        aI = cH5[format.ALN_INDEX]
        originalAttrs = aI.attrs.items()

        ## empty is a special case. In general, h5py handles
        ## zero-length slices poorly and therefore I don't want to
        ## make them. Therefore, I maintain the 'empty' variable to
        ## indicate that. This makes some code less pleasing, e.g.,
        ## computing the reference index data structure.
        if (aI.shape[0] == 0):
            log.warn("Warning: %s empty!" % inFile)
            success = True;
            return True; 
        
        # sort the AlignmentIndex
        aord = lexsort([aI[:,format.TARGET_END], aI[:,format.TARGET_START], 
                        aI[:,format.REF_ID]])

        assert(len(aord) == aI.shape[0])
        
        sAI = aI.value[aord,:]
        del(aI)
        log.msg("Sorted AlignmentIndex.")

        # construct reference offset datastructure.
        refSeqIDs = cH5[format.REF_GROUP_ID]
        offsets = computeRefIndexTable(refSeqIDs.value, sAI[:,format.REF_ID])
        log.msg("Constructed offset datastructure.")
        
        # fill overlap and back columns.
        for row in range(0, offsets.shape[0]):
            fRow = int(offsets[row, 1])
            lRow = int(offsets[row, 2])
            if (lRow - fRow <= 0):
                continue
            sAI[fRow:lRow, (format.N_BACK, format.N_OVERLAP)] = \
                computeIndicesDP(sAI[fRow:lRow, format.TARGET_START],
                                 sAI[fRow:lRow, format.TARGET_END])
        log.msg("Constructed indices.")

        # modify the cmp.h5 file.
        # We want to keep the chunking info on the dataset.
        del(cH5[format.ALN_INDEX])
        cH5.create_dataset(format.ALN_INDEX, data = sAI, dtype = h5t.NATIVE_UINT32,
                           maxshape = (None, None))
        
        ## If the file is already sorted there's no harm in resorting.
        if (__pathExists(cH5, format.REF_OFFSET_TABLE)):
            log.msg(format.REF_OFFSET_TABLE + " already exists, deleting.")
            del(cH5[format.REF_OFFSET_TABLE])

        ## create the offset datastructure in the file.
        cH5.create_dataset(format.REF_OFFSET_TABLE, data = offsets, 
                           dtype = h5t.NATIVE_UINT32, maxshape = (None, None))

        ## add the index attribute.
        cH5['/'].attrs.create("Index", ['REF_ID', 'TARGET_START', 'TARGET_END'])

        ## fixup attributes.
        for oA in originalAttrs:
            cH5[format.ALN_INDEX].attrs.create(oA[0], oA[1])

        ## deep repacking.
        if (deep):
            log.msg("Repacking alignment arrays.")
            __repackDataArrays(cH5, format, log)
            
        ## memory free.
        del sAI
        
        ## manage any extra datasets.
        for extraTable in format.extraTables:
            if (__pathExists(cH5, extraTable)):
                log.msg("Sorting table: %s" % extraTable)

                eTable = cH5[extraTable].value
                if (len(eTable.shape) == 1):
                    eTable = eTable[aord]
                else:
                    eTable = eTable[aord,:]

                ## save attributes, if any for re-writing below.
                originalAttrs = cH5[extraTable].attrs.items()

                del(cH5[extraTable])
                cH5.create_dataset(extraTable, data = eTable, 
                                   maxshape = tuple([None for x in eTable.shape]))
                for oA in originalAttrs:
                    cH5[extraTable].attrs.create(oA[0], oA[1])

        ## if you make it this far, set the flag.
        success = True

    except Exception, E:
        log.error(E)
        if (os.path.exists(outFile)):
            pass
        
    finally: 
        try:
            cH5.close()
        except:
            pass
        finally:
            return(success)


class Loggy:
    def __init__(self, level):
        self.level = level
    def write(self, msg, level):
        if (self.level >= level): sys.stderr.write(str(msg) + "\n")
    def error(self, msg): self.write(msg, 0)
    def warn(self, msg): self.write(msg, 1)
    def msg(self, msg): self.write(msg, 2)

def main():
    usage  = \
    """ %prog [options] input-file [output-file]
    
    Sort cmp.h5 files. If output-file is unspecified the input-file is
    overwritten. If there are a number of reference groups then the
    indexing processing can occur in parallel.  

    version: """ + __VERSION__
 
    parser = OptionParser(usage)
    parser.add_option("-s", "--silent", dest = "silent", action = "store_false", \
                          default = False, help = "print nothing.")
    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true", \
                          default = False, help = "print debugging information")
    parser.add_option("-d", "--deep", dest = "deep", action = "store_true", default = False, \
                          help = "whether a deep sorting should be conducted, i.e. sort the AlignmentArrays")
    parser.add_option("-j", "--jobs", dest = "jobs", default = 1, \
                          help = "Number of child processes to launch. This only speeds up processing if there are multiple references groups. Not yet Implemented.")
    parser.add_option("--tmpDir", dest = "tmpdir", default = "/tmp", \
                          help = "Temporary directory to use when sorting in-place.")
                      
                      
    (options, args) = parser.parse_args()

    if (not len(args)):
        parser.print_help()
        exit(1)
    
    infile = args[0]
    
    ## we do this in a temporary file because it is safer.
    if (len(args) < 2):
        ofile   = tempfile.NamedTemporaryFile(dir=options.tmpdir)
        outfile = ofile.name
    else:
        outfile = args[1]

    log = Loggy(2 if options.verbose else 1 if not options.silent else 0)
    success = sortCmpH5(infile, outfile, deep = options.deep, jobs = options.jobs, log = log)

    if (not success):
        log.error("Error during sorting. Exiting! Original file %s should still be intact." % infile)
        exit(1)
    else:
        ## add to the file log.
        cmpH5 = CmpH5Factory.factory.create(outfile, 'a')
        cmpH5.log("cmpH5Sort.py", __VERSION__, str(datetime.datetime.now()), ' '.join(sys.argv), "Sorting")
        cmpH5.close()

        if (len(args) < 2):
            shutil.copyfile(outfile, infile)
            ofile.close()
        exit(0)

    

if __name__ == "__main__":
     main()
    
