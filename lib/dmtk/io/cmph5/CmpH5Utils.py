"""A collection of classes and functions which operate on CmpH5 files.
Stored here as a library so they can be imported."""

import sys, os
import h5py

def repositoryPathToURI( path ):
    """Given a path to a reference repository entry (local or global), returns
    a URI describing that entry."""
    return "pbi://secondary/references/%s" % os.path.basename( path )

class PBCmpH5MetaDataError( Exception ):
    pass

class PBCmpH5MetaData( ):
    """Loads meta data into an existing PBCmpH5 file."""

    def __init__( self, fofnFN ):
        self.exps = {}
        self.runs = {}
        self.primaryPipeline = ""
        self.reportsFolder = ""
        self._load( fofnFN )

    def _getPipeline( self, plsPath ):
        """Given the path to a pls/bas file, returns the primary pipline that generated it."""
        plsH5 = h5py.File( plsPath, "r" )
        if  "PulseData" not in plsH5 or \
            "BaseCalls" not in plsH5["/PulseData"] or \
            "ChangeListID" not in plsH5["/PulseData/BaseCalls/"].attrs:
            raise PBCmpH5MetaDataError, \
                "Unable to locate /PulseData/BaseCalls/ChangeListID in pls file (%s)." % plsPath
        pipeline = plsH5["/PulseData/BaseCalls/"].attrs["ChangeListID"][ 0 ]
        plsH5.close()
        return str(pipeline)

    def _load( self, fofnFN ):
        reportsFolders = set( )
        primaryPipelines = set( )
        for line in open( fofnFN, 'r' ):
            pathElts = line.strip().split( os.sep )
            try:
                exp, run, report, plsFN = pathElts[-4:]
            except ValueError:
                raise PBCmpH5MetaDataError, "Unable to parse movie path from input.fofn (%s) required for PBCmpH5." % fofnFN
            # for pls/bas files stored in non-canonical locations
            # try to handle and move on
            try:
                exp = int(exp)
                run = int(run)
            except ValueError:
                exp = 0
                run = 0
            movie = plsFN.split(".")[0]
            self.exps[ movie ] = exp
            self.runs[ movie ] = run
            reportsFolders.add( report )
            primaryPipelines.add( self._getPipeline( line.strip() ) )

        self.reportsFolder = reportsFolders.pop()
        self.primaryPipeline = primaryPipelines.pop()
