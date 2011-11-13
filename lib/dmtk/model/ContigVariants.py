__doc__="""Encapsulates coverage information for per-contig variant plots"""
import os
import numpy as np

class ContigVariants:
    
    def __init__(self, seqId, name=None):
        """Encapsulates variant info relevant to one chart"""
        self.seqid = seqId
        if name == None:
            name = seqId
        self.name = name
        
        self.variants = []

        self.fileName = "variantsPlot_%s%s"%(self.seqid,".png")

        
    def _toDictionary(self, attributeString):
        """convert an attribute line to a dictionary"""
        return dict(item.split("=") for item in attributeString.split(";"))
    
    def addData(self, gff3Record):
        """Append x,y data from this record to the contig graph"""
        
        dict = self._toDictionary(gff3Record._getAttributeString())
        startPos = int(gff3Record.start )
        inse = int(dict['ins'])
        de1e = int( dict['del'])
        snv = int (dict['snv'])
        
        self.variants.append( (startPos, inse, de1e, snv) )
       
     
  