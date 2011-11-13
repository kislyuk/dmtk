#!/usr/bin/env python

class BedRecord:
    """Models a record in a BED file format"""
    def __init__(self, chrom='', chromStart = 0, chromEnd = 0, name = '', score = -1.00, strand = '+'):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.name = name
        self.score = score
        self.strand = strand
           
    def __str__(self):
        return '%s\t%d\t%d\t%s\t%.3f\t%s' % \
            (self.chrom, self.chromStart, self.chromEnd, self.name, \
              self.score, self.strand)
              
class BedWriter:
    """Outputs BED annotation track file"""
    def __init__(self, outfile):
        self._outfile = outfile

    def close(self):
        self._outfile.close()

    def flush(self):
        self._outfile.flush()

    def writeHeader(self, name, description, useScore):
        print >> self._outfile, 'track name=%s description="%s" useScore=%d' \
            % (name, description, useScore)

    def writeRecord(self, record):
        print >> self._outfile, str(record)
