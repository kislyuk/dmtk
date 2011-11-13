#!/usr/bin/env python
import sys
import os
import tempfile

from string import Template

class RScript( Template ):
    """Simple class for converting a template R file with parameters
    escaped by @ into a temporary R script.
    
    Example usage:
    script = RScript( 'template.R' )
    script.construct( {'infile':'some.csv', 'outfile':'output.png'} )
    temporaryFileName = script.save()
    script.delete()
    """
    delimiter = '@'

    def __init__( self, file=None ):
        self.fileName = None
        self.output = ''
        if file:
            self.load( file )

    def load( self, fileName ):
        infile = open( fileName, 'r' )
        script = "".join( infile.readlines() )
        Template.__init__( self, script )

    def construct( self, parameters ):
        self.output = self.substitute( parameters )

    def save( self, fileName=None ):
        if not fileName:
            h, fileName = tempfile.mkstemp()
            os.close( h )
        outfile = open( fileName, 'w' )
        print >>outfile, self.output
        outfile.close()
        self.fileName = fileName
        return fileName

    def delete( self ):
        if self.fileName:
            os.remove( self.fileName )
