#!/usr/bin/env python
"""Defines the Aligned Read Pulse Metric Interface. Wraps an AlignedReadH5F object
with an interface for adding and accessing pulse metrics from the aligned read data."""

import h5py
import os, sys
import numpy
import glob
from dmtk.io.AlignedReadH5IO import AlignedAstroReadH5F, AlignedReadH5F, AlignedAstroReadH5FDAO
from dmtk.io.MovieHDF5IO import Movie, ZMW, read2MoviePath, PLSH5Error, MoviePulseData
from dmtk.model.AlignmentHit import AlignmentHit

# This should be defined and imported from AlignedReadH5IO
( ALIGNMENT_ID, READ_GROUP_ID, MOVIE_ID, REF_SEQ_ID, T_START, T_END, ALIGNED_STRAND, \
    EXP_ID, RUN_ID, PANEL, X, Y, SUB_READ_ID, R_START, R_END, Z, N_M, N_MM, N_INS,
    N_DEL, OFFSET_BEGIN, OFFSET_END ) = range( 22 )

PULSE_FIELDS = [ "PulseWidth", "pkmid", "StartTime", "IPD", "Light", "ClassifierQV", "QualityValue", "DeletionQV" ]
BASEMAP = numpy.array(['-','A','C','-','G','-','-','-','T', '-','-','-','-','-','-','-','N'])

class AlignedPulsesError( Exception ):
    pass

class IPulseFeaturesH5F:
    """This interface simply describes the functionality that will be common
    to the Pulse Feature classes derrived from various AlignedRead classes.

    eg. AlignedAstroPulseH5F will implement this interface, and inherit from AlignedAstroReadH5F."""

    def loadMovies( self, listOfMoviePaths ):
        """Loads the specified .pls.h5 files, extracting relevant pulse metrics."""
        raise NotImplementedError

    def autoLoadMovies( self ):
        """Looks at internal read names to find the relevant movie files to load."""
        raise NotImplementedError


class AlignmentRow:

    def __init__( self, alignedHDF5F, alignData ):
        self._alignment = alignedHDF5F
        self._alignData = alignData

    def __getitem__( self, key ):
        return self._alignData[ key ]

    @property
    def readName( self ):
        return "_".join( [ "x%d" % self[X],
                           "y%d" % self[Y],
                           "%07d-%04d" % ( self[EXP_ID], self[RUN_ID] ),
                           str(self.movieName) ] ) + ".f%d" % self[SUB_READ_ID]

    @property
    def movieName( self ):
        id2Movie = dict( zip( self._alignment._movie2Id.values(), self._alignment._movie2Id.keys() ) )
        return id2Movie[ self[ MOVIE_ID ] ]
    
    @property
    def readSequence( self ):
        return self.alignedReadSequence.replace("-","")

    @property 
    def alignedReadSequence( self ):
        return self.alignedSequence[0]

    @property 
    def alignedSequence( self ):
        rawData = self.readGroup["AlignmentArray"]
        # explicit slice loads it into memory
        localAlnAry = rawData[ self[OFFSET_BEGIN] : self[OFFSET_END] ]
        binRBases = (localAlnAry & 0xf0) >> 4;
        binTBases = (localAlnAry & 0x0f);
        rSeq = "".join(BASEMAP[binRBases])
        tSeq = "".join(BASEMAP[binTBases])
        return ( rSeq, tSeq )
        
    @property
    def readGroup( self ):
        runName = "%07d-%04d" % ( self[ EXP_ID ], self[ RUN_ID ] )
        return self.refSeqGroup[ runName ]
        
    @property
    def refSeqGroup( self ):
        return self._alignment.getHDF5ObjectByPath("/"+self.refSeqName)
        
    @property
    def refSeqName( self ):
        refSeqIdToName = dict( zip( self._alignment._refSeqName2Id.values(), self._alignment._refSeqName2Id.keys() ) )
        return refSeqIdToName[ self[ REF_SEQ_ID ] ].strip("/")

class AlignedAstroPulseH5F( AlignedAstroReadH5F ):
    """Defines the Aligned Astro Read Pulse Metrics H5 File. Inherits from a standard
    AlignedAstroReadH5F class, adding read/write functionality for pulse metrics in
    the aligned read data."""

    def __init__( self, filename, mode ):
        """Takes in the path to a .cmp.h5 file, and a mode [rwa]."""
        AlignedAstroReadH5F.__init__( self, filename, mode )
        self._globalIndex = self._h5FileHandle["/AlignmentIndex"]
        self._alnIdByMovieId = {}

    def __iter__( self ):
        """Iterates over the individual aligned reads in the underlying
        AlignedReadH5F object."""

        return self.getAlignedReadIterator( )

    def checkPulseData( self, row, pulseData, baseMap ):
        """Given an AlignmentRow object and corresponding pulseData,
        checks to make sure that the pulse metrics match the read sequence."""

        alignedReadSeq = row.alignedReadSequence

        # Make sure the pulses correctly align with the bases, otherwise imagine the confusing downstream
        # results we'd get.
        try:
            pulseSeq = "".join( [ baseMap[pulseData["Channel"][i]] for i in range( row[R_START],row[R_END] ) ] )
        except IndexError:
            sys.stderr.write("Read: %s\n" % alignedReadSeq.replace("-",""))
            sys.stderr.write("Pulses: %s\n" % "".join( [ baseMap[pulse] for pulse in pulseData["Channel"] ]))
            sys.stderr.write("Basemap: %s\n" % str(baseMap))
            sys.stderr.write("Start, End: %d, %d\n" % ( row[R_START], row[R_END] ))
            sys.stderr.write("Len: %d\n" % len(pulseData["Channel"]))
            raise
        if pulseSeq != alignedReadSeq.replace("-",""):
            raise AlignedPulsesError, \
                "Pulses and Read do not match for %s: \nPulses:\n%s\nAligned:\n%s" % ( row.readName, pulseSeq, alignedReadSeq.replace("-","") )


    def loadPulseData( self, rows, movieData, readGroup ):
        """Given a set of AlignmentRow objects, loads in the relevant pulse data
        for all rows which are derrived from the specified readGroup and movie.
        Creates tables as necessary."""

        baseMap = movieData.channel2base
        
        for field in PULSE_FIELDS:
            
            # Fail silently if we don't have things like deletionQV.
            # The resulting cmp.h5 will simply not have this information.
            if not movieData.hasMetric( field ):
                #sys.stderr.write("Unable to locate %s metric... skipping.\n" % field )
                continue
        
            # Create a table if necessary
            try:
                pulseTable = readGroup[ field ]
            except KeyError:
                pulseTable = readGroup.create_dataset(  field,
                                                        readGroup["AlignmentArray"].shape,
                                                        dtype='f4',
                                                        maxshape=(None,) )
                aaLastRow = readGroup["AlignmentArray"].attrs["lastRow"]
                pulseTable.attrs.create( "lastRow", data=aaLastRow, dtype="u4" )
           
            ptInMemory = pulseTable[:]

            for row in rows:
                if row.readGroup.name != readGroup.name:
                    continue
                
                pulseData = movieData.getPulseDataFieldsFromXY( [ field, "Channel"], row[ X ], row[ Y ] )
                self.checkPulseData( row, pulseData, baseMap )

                # Populate the correct rows with these pulse data.
                rIndex = row[ R_START ]
                oBegin, oEnd = row[OFFSET_BEGIN], row[OFFSET_END]
                alnRow = readGroup["AlignmentArray"][oBegin:oEnd]
                dataSize = oEnd-oBegin
                toInsert = numpy.empty( dataSize, dtype=numpy.float )
                toInsert.fill( numpy.nan )
                toInsert[ BASEMAP[ (alnRow & 0xf0) >> 4 ] != "-" ] = pulseData[ field ][ rIndex: rIndex+dataSize ]
                
                ptInMemory[ oBegin:oEnd ] = toInsert
                
            pulseTable[:] = ptInMemory

    def loadMovies( self, moviePaths ):
        """Given a list of paths to valid .pls.h5 movies, loads in all pulse
        metrics from the movies that correspond to reads in these alignments,
        writing them to the underlying .cmp.h5 file."""

        lastMovieId = None
        lastMovie = None

        movieName2Path = dict( [ ( os.path.basename( mPath ).split(".")[0], mPath ) 
                                                for mPath in moviePaths ] )

        # Sort the global index by Movie ID for fast Reading
        lastRow = self._globalIndex.attrs["lastRow"]
        sortedIndices = numpy.argsort( self._globalIndex[:lastRow, MOVIE_ID] )

        # Now we want to run them in batch for fast writing.
        rowBatch = []
        readGroup = None
        for index in sortedIndices:
            
            #print >>sys.stderr, "index=%d" % index
            row = AlignmentRow( self, self._globalIndex[ index ] )
            #print >>sys.stderr, "globalIndex=%s" % str(self._globalIndex[index])

            movieId = row[ MOVIE_ID ]
            #print >>sys.stderr, "movieId=%d" % movieId
            if movieId != lastMovieId:
            
                if len(rowBatch) > 0:
                    #print >>sys.stderr, "loading pulse data %s, %s, %s" % (str(rowBatch), lastMovie, readGroup )
                    self.loadPulseData( rowBatch, lastMovie, readGroup )
                rowBatch = []
                readGroup = None

                lastMovie = MoviePulseData( movieName2Path[ row.movieName ] )
                lastMovieId = movieId
     
            if readGroup and row.readGroup.name != readGroup.name:
                raise AlignedPulsesError, \
                    "More than one read group found within a movie (%s,%s); this breaks the pulse loading code." % \
                                                        ( row.readGroup.name, readGroup.name )
            readGroup = row.readGroup

            rowBatch.append( row )

        # load last movie
        if len(rowBatch) > 0:
            #print >>sys.stderr, "loading pulse data %s, %s, %s" % (str(rowBatch), lastMovie, readGroup )
            self.loadPulseData( rowBatch, lastMovie, readGroup )

    def autoLoadMovies( self, reportsGlob="*DragonRabbit*", root="/mnt/data*/vol*/"):
        """Uses the read data and an assumed filesystem structure (/mnt/data*/vol*/...)
        to load in pulse metric data for all of the reads in this alignment."""
        
        moviePaths = self.getMoviePaths( reportsGlob, root )
        self.loadMovies( moviePaths )

    def getMoviePaths( self, reportsGlob="SideshowBob", root="/mnt/data*/vol*/"):
        """Uses the given inputs, as well as read data, to search for
        and return paths to all of the .pls.h5 files that went into
        making these reads."""

        moviePaths = set()
        movieGlobs = set()
        MOVIE_GLOB = root + "/%.7d/%.4d/*%s*/*p%d*.pls.h5"

        for row in self._globalIndex:
            mGlob = MOVIE_GLOB % ( int(row[EXP_ID]), int(row[RUN_ID]), reportsGlob, int(row[PANEL]) )
            movieGlobs.add( mGlob )

        for mGlob in movieGlobs:
            # The logic here is as follows:
            # Grab all matches to the permissive glob (*ReportsFolderName*).
            # If any of the matches end with the exact string, use only those, else use all.
            matchPaths = glob.glob( mGlob )
            if any( os.path.split( path )[-2].endswith( reportsGlob ) for path in matchPaths ):
                for path in matchPaths:
                    if os.path.split( path )[-2].endswith( reportsGlob ):
                        moviePaths.add( path )
            else:
                for path in matchPaths:
                    moviePaths.add( path )
               
        return moviePaths 

class AlignmentHitWithPulses( AlignmentHit ):
    """Wraps an AlignmentHit, adding pulse metrics."""
    
    def __init__( self ):
        """Set up some public attributes where we'll store the 
        pulse metrics for now."""

        AlignmentHit.__init__( self )
        self.pulses = {}   

    def __str__( self ):
        """Returns the string representation of the base class,
        plus a rather terse summary of the existing pulse metrics."""

        baseStr = AlignmentHit.__str__( self ) 
        for metric in PULSE_FIELDS:
            if metric in self.pulses:
                baseStr += "\n" + metric + " Count: " + str(len( [ 1 for pm in self.pulses[ metric ] if pm != numpy.nan ] ))
        return baseStr+"\n"

    def __getitem__( self, key ):
        """Given a pulse metric name, returns the array
        of values associated with that metric (by alignment coordinate).
        
        Otherwise it passes the key along to the base class."""
        
        if key in self.pulses:
            return self.pulses[ key ]
        else:
            return AlignmentHit.__getitem__( key )

class AlignedAstroPulseH5FDAO( AlignedAstroReadH5FDAO ):
    """Extends the functionality of the AlignedAstroReadDAO to the AlignedAstroPulse."""

    def __init__( self, h5f ):
        """Takes in an AlignedAstroPulseH5F"""
        AlignedAstroReadH5FDAO.__init__( self, h5f )

    def _getPerBaseInfo( self, readGroup ):
        """Returns the per base info of the base class, plus
        the pulse metrics."""

        pbi = AlignedAstroReadH5FDAO._getPerBaseInfo( self, readGroup )
        for metric in PULSE_FIELDS:
            try:
                pbi[ metric ] = readGroup[ metric ][:] # load into memory
            except KeyError:
                pass # If we can't find the pulse info, we just don't use it.

        return pbi

    def _hitFromPerBaseInfo( self, perBaseInfo, offsetBegin, offsetEnd, hitToModify=None ):
        """Generates a new AlignmentHitWithPulses, unless hitToModify
        is specified. Adds pulse metrics, and returns."""

        hit = hitToModify if hitToModify else AlignmentHitWithPulses()
        hit = AlignedAstroReadH5FDAO._hitFromPerBaseInfo( self, perBaseInfo, offsetBegin, offsetEnd, hit )
       
        for metric in PULSE_FIELDS:
            if metric in perBaseInfo:
                hit.pulses[ metric ] = perBaseInfo[ metric ][ offsetBegin:offsetEnd ]

        return hit
        

    def _getColumnIterator( self, refSeqName, refPosition, exportRefOffset=False ):
        """Returns a generator of alignment columns."""

        raise NotImplementedError, "Column Iterator not yet implemented for AlignedAstroPulseH5FDAO"

    
def run( ):
   
    align = AlignedAstroPulseH5F( sys.argv[1], "a" )
    align.autoLoadMovies( reportsGlob=sys.argv[2] )

    #for alignment in AlignedAstroPulseH5FDAO( align ).getAlignedReadIterator( sys.argv[3] ):
    #    print str(alignment)


if __name__ == "__main__":

    import cProfile
    import sys
    
    #run()
    cProfile.run( 'run()', 'profiling.txt' )
