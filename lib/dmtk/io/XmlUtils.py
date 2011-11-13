__doc__="""Simple utility classes for hand-coded XML marshalling and un-marshalling."""
import sys
import os
import cStringIO
#import libxml2
#import libxslt

from xml.sax import *
from xml.sax.handler import ContentHandler
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl

# using Saxon for XSLT2 support
SAXON_PATH = 'saxonb9'

class PrettyXMLGenerator( XMLGenerator ):
    def __init__( self, file ):
        XMLGenerator.__init__( self, file )
        self.indent = 0

    def startElement( self, name, attrs ):
        self.ignorableWhitespace( " " * self.indent )
        XMLGenerator.startElement( self, name, attrs )
        self.ignorableWhitespace( os.linesep )
        self.indent += 2

    def endElement( self, name ):
        self.indent -= 2
        self.ignorableWhitespace( " " * self.indent )
        XMLGenerator.endElement( self, name )
        self.ignorableWhitespace( os.linesep )

class XmlMarshaller:
    EMPTY_ATTRS = AttributesImpl( {} )

    def __init__( self ):
        self.out = cStringIO.StringIO()
        self.generator = PrettyXMLGenerator( self.out )
        self.xmlString = None
        self.nameStack = []

    def __str__( self ):
        return self.xmlString

    def emptyElement( self, name, attrs={} ):
        self.generator.startElement( name, AttributesImpl( attrs ) )
        # save space for empty elements by using />
        self.generator.indent -= 2
        self.out.seek( -2, os.SEEK_CUR )
        self.out.write( '/>' + os.linesep )

    def startElement( self, name, attrs={} ):
        if len(attrs)==0:
            a = XmlMarshaller.EMPTY_ATTRS
        else:
            a = AttributesImpl( attrs )
        self.generator.startElement( name, a )
        self.nameStack.append( name )

    def endElement( self ):
        name = self.nameStack.pop()
        self.generator.endElement( name )

    def characters( self, seq ):
        for ch in seq:
            self.generator.characters( ch )

    def start( self ):
        self.generator.startDocument()

    def end( self ):
        self.xmlString = self.out.getvalue()
        self.out.close()

def xslTransform( xmlFileName, xslFileName, outputFileName, dontRun=False ):
    # using saxon for XSLT 2 support
    cmdLine = '%s -xsl:%s -s:%s -o:%s' % ( SAXON_PATH, \
                                           xslFileName, xmlFileName, outputFileName )
    if dontRun:
        return cmdLine
    else:
        os.system( cmdLine )
    
    # The following code works fine with libxml2 bindings
    # ...but can't handle XSLT 2
    #styleDoc = libxml2.parseFile( xslFileName )
    #style = libxslt.parseStylesheetDoc( styleDoc )
    
    #resultsDoc = libxml2.parseFile( xmlFileName )
    #reportDoc = style.applyStylesheet( resultsDoc, None )
    
    #if isinstance( outputFile, str ):
    #    outputFile = open( outputFile, 'w' )
    #style.saveResultToFd( outputFile.fileno(), reportDoc )

    #reportDoc.freeDoc()
    #resultsDoc.freeDoc()
    #styleDoc.freeDoc()
    #if outputFile.fileno()!=sys.stdout.fileno() and \
    #    outputFile.fileno()!=sys.stderr.fileno():
    #    outputFile.close()

if __name__=='__main__':
    pass
