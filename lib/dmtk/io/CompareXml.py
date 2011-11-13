import sys
import os
import traceback

from xml.sax import *
from xml.sax.handler import ContentHandler
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl

from xml.etree.cElementTree import iterparse

from dmtk.model.AlignmentHit import AlignmentHit

class ICompareHitHandler:
    """ Callback interface for gathering hit elements as they are
    parsed from a compare XML file. """
    def __init__(self):
        pass
    
    def processHit( self, compareHit ):
        pass

class DefaultCompareHitHandler(ICompareHitHandler):
    def __init__(self):
        self.hits = []
        
    def processHit( self, compareHit ):
        self.hits.append( compareHit )
    
    def getHits(self):
        return self.hits

class CompareHit( AlignmentHit ):
    def __init__(self):
        AlignmentHit.__init__(self)
        self.alignedQuery = ''
        self.alignedTarget = ''
        self.zScore=-100.0

class CompareXmlHandler( ContentHandler ):
    def __init__( self, handler=None ):
        self.currentString = []
        self.inHit = False
        self.currentHit = None
        self.handler = handler
      
    def startElement( self, name, attrs ):
        if name=='hit':
            self.currentHit = CompareHit()
            self.currentHit.query_id = attrs.getValue('name')
            self.currentHit.query_start = int( attrs.getValue('start' ) )
            self.currentHit.query_end = int( attrs.getValue('end') )
            self.currentHit.target_start = int( attrs.getValue( 'targetStart' ) )
            self.currentHit.target_end = int( attrs.getValue( 'targetEnd' ) )
            self.currentHit.query_strand = attrs.getValue('strand')
            self.currentHit.target_strand = attrs.getValue('targetStrand')
            self.inHit = True
        elif name=='query':
            self.clearString()
        elif name=='target':
            self.clearString()
        elif name=='zScore':
            self.currentHit.zScore = float( attrs.getValue('value') )
        elif name=='score':
            self.currentHit.score = attrs.getValue('value')
        elif name=='nCorrect':
            self.currentHit.pctCorrect = float(attrs.getValue('percent'))
        elif name=='nInsert':
            self.currentHit.nIns = int(attrs.getValue('value'))
        elif name=='nDelete':
            self.currentHit.nDel = int(attrs.getValue('value'))            
        elif name=='nMismatch':
            self.currentHit.nMismatch = int(attrs.getValue('value'))
        
    def endElement( self, name ):
        if name=='hit':
            if self.handler:
                self.handler.processHit( self.currentHit )
            self.inHit = False
        elif name=='query':
            # convert to 8-bit string for efficiency and compatibility with
            # downstream expectations (no loss here)
            self.currentHit.alignedQuery = self.getCurrentString().encode('utf-8')
        elif name=='target':
            # convert to 8-bit string for efficiency and compatibility with
            # downstream expectations (no loss here)
            self.currentHit.alignedTarget = self.getCurrentString().encode('utf-8')
            
    def getCurrentString( self ):
        return "".join(self.currentString).strip()

    def clearString( self ):
        self.currentString = []

    def characters( self, ch ):
        self.currentString.append(ch)
        
def parseCompareFile( compareFile, hitHandler ):
    try:
        parser = make_parser()
        handler = CompareXmlHandler( hitHandler )
        parser.setContentHandler( handler )
        infile = open( compareFile, 'r' )
        parser.parse( infile )
    except SAXException, se:
        sys.stderr.write( se.getMessage() + os.linesep )
    except Exception, e:
        sys.stderr.write( str(e) + os.linesep )
        traceback.print_exc( file=sys.stderr )
    finally:
        if infile: infile.close()
        
def iterparseCompareFile( compareFile ):
    """Returns an iterator over hit Elements (from the ElementTree API).
    Very memory efficient."""
    infile = open( compareFile, 'r' )
    context = iter(iterparse( infile, ('start','end') ))
    event, root = context.next()
    for event, elem in context:
        if event=='end' and elem.tag=='hit':
            yield elem
            elem.clear()
            root.clear()
    infile.close()
