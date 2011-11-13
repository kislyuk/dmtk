#!/usr/bin/env python

class VcfRecord:
    """Models a record in a VCF3.3 file."""
    def __init__(self):
        self.chrom = ''
        self.pos = 1
        self.id = '.'
        self.ref = ''
        self.alt = ''
        self.qual = -1.00
        self.filter = '0'
        self.info = {}
     
    def put(self, key, value):
        self.info[key] = value
   
    @staticmethod
    def getHeader():
        return 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
             
    def _getInfoString(self):
        return ';'.join(['%s=%s' % (k,v) \
            for k,v in self.info.iteritems()])

    def __str__(self):
        return '%s\t%d\t%s\t%s\t%s\t%.2f\t%s\t%s' % \
            (self.chrom, self.pos, self.id, self.ref, \
              self.alt, self.qual, self.filter, self._getInfoString())
              
class VcfWriter:
    """Outputs VCF (1000 Genomes Variant Call Format) 3.3 files"""
    def __init__(self, outfile):
        self._outfile = outfile
        self._start()

    def close(self):
        self._outfile.close()

    def flush(self):
        self._outfile.flush()

    def _start(self):
        self.writeMetaData('fileformat', 'VCFv3.3')

    def writeHeader(self):
        print >> self._outfile, '#%s' % VcfRecord.getHeader()
        
    def writeMetaData(self, key, value):
        print >> self._outfile, '##%s=%s' % (key, value)

    def writeRecord( self, record ):
        print >> self._outfile, str(record)
