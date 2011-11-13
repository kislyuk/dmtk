import sys
import os
import traceback

from xml.sax import *
from xml.sax.handler import ContentHandler
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl

from dmtk.model.VariantItem import *

class IVariantItemHandler:
    """ Callback interface for gathering variant elements as they are
    parsed from a variant XML file. """
    def __init__(self):
        pass
    
    def processVariant( self, variant ):
        pass

class DefaultVariantItemHandler(IVariantItemHandler):
    def __init__(self):
        self.variants = []
        
    def processVariant( self, variant ):
        self.variants.append( variant )
    
    def getVariants(self):
        return self.variants

class VariantXmlHandler( ContentHandler ):
    def __init__( self, handler=None ):
        self.currentVariant = None
        self.currentEvidence = None
        self.handler = handler
      
    def startElement( self, name, attrs ):
        if name=='variant':
            self.currentVariant = VariantItem()
            self.currentVariant.position = int(attrs.getValue('pos'))
            self.currentVariant.type = attrs.getValue('type')
            self.currentVariant.iub = attrs.getValue('iub')
            self.currentVariant.length = int(attrs.getValue('length'))
        elif name=='evidence':
            self.currentEvidence = VariantEvidence()
            self.currentVariant.evidence = self.currentEvidence
        elif name=='coverage':
            self.currentEvidence.fwdCoverage = int(attrs.getValue('fwd'))
            self.currentEvidence.revCoverage = int(attrs.getValue('rev'))
        elif name=='callComplexity':
            self.currentEvidence.callComplexity = float(attrs.getValue('value'))

    def endElement( self, name ):
        if name=='variant':
            if self.handler:
                self.handler.processVariant( self.currentVariant )
            
def parseVariantsFile( variantFile, variantHandler ):
    try:
        parser = make_parser()
        handler = VariantXmlHandler( variantHandler )
        parser.setContentHandler( handler )
        infile = open( variantFile, 'r' )
        parser.parse( infile )
    except SAXException, se:
        sys.stderr.write( se.getMessage() + os.linesep )
    except Exception, e:
        sys.stderr.write( str(e) + os.linesep )
        traceback.print_exc( file=sys.stderr )
    finally:
        if infile: infile.close()
