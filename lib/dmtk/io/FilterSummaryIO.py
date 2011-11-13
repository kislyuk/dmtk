__doc__="""Reads a FilterSummary, returns data as an array  """
import sys
import os
import numpy as np

# In order to support variation in columns present, we need to specify their types
# Format is { Name: ( dType, f(str->dType) ) }
VALID_COLUMNS = {   "Movie"         : ( "|S64", str ),
                    "ReadId"        : ( "|S64", str ),
                    "#Bases"        : ( int,    int ),
                    "Readlength"    : ( int,    int ),
                    "ReadScore"     : ( float,  float ),
                    "Productivity"  : ( int,    int ),
                    "SequencingZMW" : ( int,    int ),
                    "PassedFilter"  : ( int,    int ),
                    "Sandwiches"    : ( int,    int ),
                    "Whitelisted"   : ( int,    int ),
                    "SNR"           : ( float, float ),
                    "ArtifactScore" : ( "|S64", str )
                }

class FilterSummary:
    
    
    
    def __init__(self, path ):
        """Path to filterSummary csv file"""
        self._path = path
        self._titles = None
        self._data = None
        self._load( )
       

    def _load(self):
        """1-time loading of the data into an array"""
        self._titles = None
        self._data = []
        
        for line in open( self._path, "r" ):
            elts = line.strip().split(",")
            if self._titles == None:
                self._titles = elts
            else:
                elts = [ VALID_COLUMNS[ name ][ 1 ]( e ) for ( name, e ) in zip( self._titles, elts ) ]
                self._data.append( tuple(elts) )

        
    
    def getDataAsNPArray(self, seqZmwsOnly=True):
        dTypes = [ ( name, VALID_COLUMNS[ name ][ 0 ] ) for name in self._titles ]
        return np.array( self._data, dtype=dTypes )

    
    
    def hasData(self):
        """returns true if data file exists and has records"""
        return len(self._data) > 0
    
    def numReads(self):
        """Get the total number of reads, filtered and unfiltered"""
        npa = self.getDataAsNPArray()
        return len( npa[ npa["SequencingZMW"] > 0 ] )
    
    def numFilteredReads(self):
        """Get the number of filtered reads"""
        npa = self.getDataAsNPArray()
        return len( npa[ npa["PassedFilter"] > 0 ] )

            
    
        
        
