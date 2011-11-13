#!/usr/bin/env python

SEQUENCE_SOURCE = 'SMRT'

import gzip
import os
import sys
import time
from dmtk.io.VcfIO import VcfRecord
from dmtk.io.BedIO import BedRecord
from numpy import median, zeros
from dmtk.io.AmosIO import amosQvCharToInt

class Gff3Record:
    """Models a record in a GFF3 file."""
    def __init__( self, hit=None, seqName='' ):
        if not hit:
            self.seqid = seqName
            self.source = '.'
            self.type = '' 
            self.start = 0
            self.end = 0
            self.score = 0.0
            self.strand = '+'
            self.phase = '.'
            self.attributes = {}
        else:
            self.__initialize( hit, seqName )

    def __getitem__( self, key ):
        return self.attributes[key]

    def __initialize( self, hit, seqName ):
        self.seqid = seqName
        self.source = SEQUENCE_SOURCE
        self.type = 'read'
        # In GFF start and end are 1-based
        # start must be less than or equal to end
        # end is inclusive
        self.start = hit.target_start + 1
        self.end = hit.target_end
        if self.start > self.end:
            t = self.end
            self.end = self.start
            self.start = t
        self.score = hit.zScore
        self.strand = hit.target_strand
        self.phase = '.'
        self.attributes = {}

    def put( self, key, value ):
        self.attributes[ key ] = value

    def _getAttributeString( self ):
        return ';'.join( [ '%s=%s' % (k,v) \
            for k,v in self.attributes.iteritems() ] )

    def __str__( self ):
        return '%s\t%s\t%s\t%d\t%d\t%.2f\t%s\t%s\t%s' % \
            ( self.seqid, self.source, self.type, self.start, \
              self.end, self.score, self.strand, self.phase, \
              self._getAttributeString() )

class CoverageGff3Record(Gff3Record):
    """Subclasses the generic Gff3Record to enable utilities specific to genome coverage"""
    def toBedRecord(self):
        bed = BedRecord()
        bed.chrom = self.seqid
        bed.chromStart = self.start - 1
        bed.chromEnd = self.end 
        bed.name = 'meanCov'
        bed.score = float(self.attributes['cov2'].split(',')[0])
        bed.strand = self.strand                   
        
        return bed
        
class VariantGff3Record(Gff3Record):
    """Subclasses the generic Gff3Record to enable utilities specific to variants"""
    def toVcfRecord(self):
        vcf = VcfRecord()
        vcf.chrom = self.seqid
        vcf.id = '.'
        vcf.ref = self.attributes['reference'].upper()
        vcf.qual = float(self.attributes['confidence'])
        vcf.filter = '0'
        vcf.put('NS', 1)
        vcf.put('DP', self.attributes['coverage'])
    
        feature = self.type
        if feature == 'insertion':
            vcf.alt = 'I%s' % self.attributes['genotype'].upper()
            #GFF associates insertions with the base after the event, while
            #VCF associates insertions with the base before the event
            vcf.pos = self.start - 1
        elif feature == 'deletion':
            vcf.alt = 'D%s' % self.attributes['length']
            vcf.pos = self.start
        elif feature == 'SNV':
            vcf.alt = self.attributes['genotype'].upper()
            vcf.pos = self.start
        elif feature == 'gap':
            #There really is no way to represent discontinuous regions
            #in VCF, so we are going to ignore these features in the GFF
            #for now.
            return None
        else:
            print >> sys.stderr, 'Unsupported feature %s found in GFF3 file.' % feature    
        
        return vcf  
    
    def toBedRecord(self):
        bed = BedRecord()
        bed.chrom = self.seqid
        bed.chromStart = self.start - 1
        bed.score = float(self.attributes['confidence'])
        bed.strand = self.strand
        
        feature = self.type
        #GFF3 coordinates are 1-based and inclusive
        #BED coordinates are 0-based and exclusive 
        if feature == 'insertion':
            bed.chromEnd = bed.chromStart + 1
            bed.name = '%d_%dins%s' % (bed.chromStart + 1, bed.chromEnd + 1, self.attributes['genotype'].upper())
        elif feature == 'deletion':
            featureLen = int(self.attributes['length'])
            bed.chromEnd = bed.chromStart + featureLen
            if featureLen == 1:
                bed.name = "%ddel" % (bed.chromStart + 1)
            else:
                bed.name = '%d_%ddel' % (bed.chromStart + 1, bed.chromEnd)
        elif feature == 'SNV':
            bed.chromEnd = bed.chromStart + 1
            bed.name = '%d%s>%s' % (bed.chromStart + 1, self.attributes['reference'].upper(), self.attributes['genotype'].upper())
        elif feature == 'gap':
            return None
        else:
            print >> sys.stderr, 'Unsupported feature %s found in GFF3 file.' % feature
            
        return bed    

class GffDeNovoRecordFactory:
    RECORD_COUNT = 0

    def __init__(self, contig, hits, multiRef=False, offset=0, regionSize=100, parseTarget=False):
        self._contig = contig
        self._hits = hits
        self._regionSize = regionSize
        self._setupCoverageCount()
        self._offset = offset
        self._multiRef = multiRef
        GffDeNovoRecordFactory.RECORD_COUNT += 1
        self._seqName = "ref%06i" % int(hits[0].target_id.split("_")[0])    \
                        if (multiRef and parseTarget) else                  \
                        "ref%06i" % GffDeNovoRecordFactory.RECORD_COUNT     \
                        if multiRef else                                    \
                        "ref000001" 

    def _setupCoverageCount(self):
        self.coverages = zeros( self._contig.getUngappedLength() )
        for h in self._hits:
            self.coverages[h.target_start:h.target_end] += 1

    def getRecords(self):

        for regionStart in range( 0, self._contig.getUngappedLength(), self._regionSize ):
            regionEnd = min(self._contig.getUngappedLength(), regionStart + self._regionSize) 
            record = GffDeNovoRegionRecord( self._contig, regionStart, regionEnd, \
                       self.coverages, multiRef=self._multiRef, seqName=self._seqName )
            record.start += self._offset
            record.end   += self._offset
            yield record

        yield GffDeNovoContigRecord( self._contig, self.coverages, seqName=self._seqName )

class GffDeNovoRegionRecord(Gff3Record):

    QV_CORRECTION = 10.0
    def __init__(self, contig, start, end, coverages, multiRef=False, seqName="ref000001" ):
        Gff3Record.__init__(self, seqName=seqName)

        self._name = "%i_contig" % contig.getIID()
        coverages = coverages[start:end]
        self.put("cov", "%i,%i,%i" % (min(coverages), median(coverages), max(coverages)))

        qualities = map(amosQvCharToInt, list(contig.getUngappedQualString()) )
        qualities = qualities[start:end]
        div = GffDeNovoRegionRecord.QV_CORRECTION
        self.put("cQv", "%.3g,%.3g,%.3g" % ( min(qualities)/div, median(qualities)/div, max(qualities)/div) )
        self.put("id", self._name)

        regionString = contig.getSeqString()[ contig.ungap2gap(start) + 1 : contig.ungap2gap(end) + 1 ]
        # note this it the number of apparent deletions (or incorrect insertion in reads) in the contig consensus
        self.put( "ins", "%i" % regionString.count("-") )
        # TODO put these somewhere?
        # self.put("del", "%i" % contig.getNumTiles() )
        # self.put("svn", "%i" %  \
        #        sum([ tile.range.end - tile.range.begin for tile in             \
        #              map(contig.getNthReadTile, range(contig.getNumTiles())) ]))
        self.start = start + 1 # 1-based, inclusive
        self.end = end
        self.type = "region"

class GffDeNovoContigRecord(GffDeNovoRegionRecord):
    def __init__(self, contig, coverages, seqName="ref000001" ):
        GffDeNovoRegionRecord.__init__(self, contig, 0, contig.getUngappedLength(), coverages, \
                                       multiRef=True, seqName=seqName)
        self.put("metrics", "%s,%s,%i" % (self.attributes["cov"], self.attributes["cQv"], contig.getNumTiles()) )
        self.type = "contig"

class GapGffDeNovoRecord(Gff3Record):
    """Creates a gap gff record starting at start and ending at end"""

    RECORD_COUNT = 0

    def __init__(self, start, end, multiRef=False):
        GapGffDeNovoRecord.RECORD_COUNT += 1
        seqName = "ref%06i" % GapGffDeNovoRecord.RECORD_COUNT \
                    if multiRef else "ref000001" 
        Gff3Record.__init__(self, seqName=seqName)
        self.start = start
        self.end   = end
        self.put("gaps", "1,%i" % (end - start + 1) )
        self.type = "region"

class GffDeNovoRecord(Gff3Record):
    """Creates a gff record from a list of hits from a single de novo contig."""

    RECORD_COUNT = 0
    # for now we divide quality values by 10.0 (since primaries QVs are off)
    QV_CORRECTION = 10.0

    def __init__(self, contig, hits, multiRef=False ):
        GffDeNovoRecord.RECORD_COUNT += 1
        seqName = "ref%06i" % GffDeNovoRecord.RECORD_COUNT \
                    if multiRef else "ref000001" 
        Gff3Record.__init__(self, seqName=seqName)
        self._hits = hits
        self._contig = contig
        self.start = min(map(lambda h: h.target_start, self._hits)) + 1
        self.end   = max(map(lambda h: h.target_end, self._hits))

        self._name = "%i_contig" % self._contig.getIID()
        self.type = "region"

        self._setAttributes()

    def _setAttributes(self):

        # TODO inefficient, but works for now
        coverages = [0] * (self.end - self.start)
        for h in self._hits:
            for idx in range(h.target_start, h.target_end):
                coverages[idx - self.start] += 1

        self.put("cov", "%i,%i,%i" % (min(coverages), sum(coverages) / 
                                      len(coverages), max(coverages)))

        qualities = map(amosQvCharToInt, list(self._contig.getUngappedQualString()) ) 
        div = GffDeNovoRecord.QV_CORRECTION
        self.put("cQv", "%.3g,%.3g,%.3g" % ( min(qualities)/div, median(qualities)/div, max(qualities)/div) )
        self.put("id", self._name)

    
def parseGffLine( line ):
    values = line.split('\t')
    record = Gff3Record()
    try:
        
        record.seqid = values[0]
        record.source = values[1]
        record.type  = values[2]
        record.start = int(values[3])
        record.end = int(values[4])
        if values[5]=='.':
            # this doesn't really keep with the semantics of GFF
            # (but what kind of format specification allows float and null?)
            record.score = 0.0
        else:
            record.score = float(values[5])
        record.strand = values[6]
        record.phase = values[7]
        record.attributes = {}
        for kvPair in values[8].split(';'):
            vals = kvPair.split('=')
            if len(vals)==2:
                record.attributes[ vals[0] ] = vals[1]
        return record
    except IndexError, e:
        print >> sys.stderr, 'Unexpected GFF line format: %s' % values


class GffReader:
    def __init__( self, fileName ):
        self.fileName = fileName
        if self.fileName.endswith(".gz"):
            self.infile = gzip.open( self.fileName, 'r' )
        else:
            self.infile = open( self.fileName, 'r' )
        self.seqMap = {}
        
    def __iter__( self ):
        for line in self.infile:
            if line[0]=='#':
                splitFields  = line.replace('#', '').split(' ')
                field = splitFields[0]
                value = " ".join( splitFields[1:] )
                if field == 'sequence-header':
                    [internalTag, delim, externalTag] = value.strip().partition(' ')
                    self.seqMap[internalTag] = externalTag
                continue
            record = parseGffLine( line[:-1] )
            if record is None:
                continue
            yield record

    def close( self ):
        self.infile.close()

class GffWriter:
    """Simple class for producing a GFF3 file from python code"""
    def __init__( self, outfile ):
        self._outfile = outfile
        self._start()

    def close( self ):
        self._outfile.close()

    def flush( self ):
        self._outfile.flush()

    def _start( self ):
        self.writeMetaData( 'gff-version', '3' )

    def writeMetaData( self, key, value ):
        # per the GFF3 spec: meta-data should be written as a
        # space-delimited key-value pair, not tab-delimited
        # (addresses bug 11270)
        print >>self._outfile, '##%s %s' % ( key, value )

    def writeRecord( self, record ):
        print >>self._outfile, str(record)

class DeNovoGffWriter(GffWriter):

    def __init__( self, outfile, commandLine="" ):
        self._columnLabels = [ "Min. Cov.", "Med. Cov.", "Max. Cov.", "Min. QV", "Med. QV", "Max. QV", "Num. Reads" ] 
        self._commandLine = commandLine
        GffWriter.__init__( self, outfile )

    def _start( self ):
        GffWriter._start( self )
        self.writeMetaData( 'date', time.ctime().replace(" ", "-") )
        self.writeMetaData( 'source', 'Allora v0.1' )
        self.writeMetaData( 'source-commandline', self._commandLine )
        self.writeMetaData( 'contig-metrics', '"%s"' % '","'.join(self._columnLabels) )
