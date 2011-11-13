__doc__="""Collection of classes and routines for simple FASTA parsing.
    (Consider BioPython for more serious needs)."""
import sys
import os

from itertools import *
from optparse import OptionParser

from dmtk.model.AlignedSequence import AlignedSequence
from dmtk.model.MSA import MSA
from dmtk.model.Range import Range


class FastaStreamReader:
    """Useful for parsing fasta streams as opposed to files"""
    DEFAULT_DELIMITER = '>'
    
    def __init__(self, iterator):
        self.iterator = iterator
        self.delimiter = FastaStreamReader.DEFAULT_DELIMITER

    def setDelimiter(self, delimiter):
        self.delimiter = delimiter

    def __iter__( self ):
        name = ''
        seqBuffer = []
        for line in self.iterator:
            if len(line)<1: continue
            if line[0]=='#' and self.delimiter!='#': continue
            line = line.rstrip()
            if len(line)==0: continue
            if line[0]==self.delimiter:
                if len(seqBuffer)>0:
                    yield FastaEntry( name, "".join(seqBuffer) )
                seqBuffer = []
                if len(line)==1:
                    name = ''
                else:
                    name = line[1:]
            else:
                seqBuffer.append( line )
        if len(seqBuffer)>0:
            yield FastaEntry( name, "".join(seqBuffer) )


class SimpleFastaReader:
    
    def __init__( self, fileName ):
        self.fileName = fileName
        if not os.path.exists( self.fileName ):
            sys.stderr.write( "Can't find file %s\n" % fileName )
            raise IOError, "FastaIO: can't find file %s" % fileName
        self.infile = open( self.fileName, 'r' )
        self.streamReader = FastaStreamReader(self.infile)

    def __iter__(self):
        return self.streamReader.__iter__()

    def setDelimiter(self, delimiter):
        self.streamReader.setDelimiter(delimiter)

    def close(self):
        self.infile.close()

class MsaReader:
    def __init__(self, fileName, ignoreFirst=False, iterate=False, keepTermGaps=False, referenceFirst=False):
        self._fileName = fileName
        self.ignoreFirst = ignoreFirst
        if not iterate:
            self._loadFromFile(self._fileName,ignoreFirst=ignoreFirst, keepTermGaps=keepTermGaps, referenceFirst=referenceFirst)
        
    def _loadFromFile(self, fileName, ignoreFirst=False, keepTermGaps=False, referenceFirst=False):
        reader = SimpleFastaReader(fileName)
        self.msa = MSA()
        first = True
        for entry in reader:
            if ignoreFirst and first:
                first = False
                continue

            aSeq = AlignedSequence( entry.sequence, aligned=True, keepTermGaps=keepTermGaps )

            sStart = entry.getAnnotation('start')
            if sStart: aSeq.alignedRange.addDelta( int(sStart) )

            strand = entry.getAnnotation('strand')
            if not strand: strand = '+'

            isReference = False
            if referenceFirst == True and first:
                first = False
                isReference = True

            self.msa.add( aSeq, name=entry.name, strand=strand, isReference=isReference  )

        reader.close()

    def nameStrandSeqIter( self ):
        reader = SimpleFastaReader( self._fileName )
        first = True
        for entry in reader:
            if self.ignoreFirst and first:
                first = False
                continue
            aSeq = AlignedSequence( entry.sequence, aligned=True )
            sStart = entry.getAnnotation('start')
            if sStart:
                aSeq.alignedRange.addDelta( int(sStart) )
            strand = entry.getAnnotation('strand')
            if not strand:
                strand = '+'
            yield entry.name, strand, aSeq
        reader.close()

    def nameStrandCrSeqIter( self ):
        reader = SimpleFastaReader( self._fileName )
        first = True
        for entry in reader:
            if self.ignoreFirst and first:
                first = False
                continue
            aSeq = AlignedSequence( entry.sequence, aligned=True )
            sStart = entry.getAnnotation('start')
            if sStart:
                aSeq.alignedRange.addDelta( int(sStart) )
            strand = entry.getAnnotation('strand')
            if not strand:
                strand = '+'
            cr = entry.getAnnotation('cr')
            yield entry.name, strand, cr, aSeq
        reader.close()
        
    def getMsa(self):
        return self.msa
    
class FastaEntry:
    """Storage class for modeling a named sequence stored in a FASTA file.
    Supports 'extended-FASTA' notation for key-value annotations."""
    def __init__( self, name, sequence ):
        self.sequence = sequence
        self.raw_name = name
        self.__processAnnotations( name )

    def getAnnotation(self,key):
        if key in self._annotations:
            return self._annotations[key]
        return None

    def getTag(self):
        if len(self._annotations)==0:
            return self.name
        tag = '%s|%s' % ( self.name, \
            '|'.join(['%s=%s'%(k,v) for k,v in self._annotations.iteritems()] ) )
        return tag

    def __processAnnotations(self, tag):
        """Processes extended syntax for entry names of the form
        >readName|key1=value1|key2=value2|..."""
        self._annotations = {}
        if tag.find('|')<0:
            self.name = tag
            return
        pairs = tag.split( '|' )
        self.name = pairs[0]
        for pair in pairs[1:]:
            if '=' not in pair:
                self.name = '%s|%s' % ( self.name, pair )
                continue
            values = pair.split('=')
            if len(values)==2:
                self._annotations[ values[0] ] = values[1]
        # revert to traditional model if this tag doesn't have kv pairs
        if len(self._annotations)==0:
            self.name = tag

    def __str__( self ):
        buffer = []
        buffer.append( ">" + self.getTag() )
        buffer.append( prettyprint(self.sequence) )
        return os.linesep.join(buffer)
    
    def subseq(self, seqRange, name=None):
        if not name:
            name = "%s_%i_%i" % (self.name, seqRange.getStart(), seqRange.getEnd())
        return FastaEntry(name, self.sequence[seqRange.getStart():seqRange.getEnd()] )

class SimpleFastqReader:
    def __init__( self, fileName ):
        self.fileName = fileName
        if not os.path.exists( self.fileName ):
            sys.stderr.write( "Can't find file %s\n" % fileName )
            raise IOError, "Fastq IO: can't find file %s" % fileName
        self.infile = open( self.fileName, 'r' )
        
    def __iter__( self ):
        qualFlag = False
        seqBuffer = []
        qualBuffer = []
        for line in self.infile:
            if len(line)<1: continue
            if line[0]=='#' and self.delimiter!='#': continue
            line = line.rstrip()
            if len(line)==0: continue
            if line[0]=="@" and len(seqBuffer) == len(qualBuffer):
                if len(seqBuffer)>0:
                    yield FastqEntry( name, "".join(seqBuffer), "".join(qualBuffer) )
                qualFlag = False
                seqBuffer = []
                qualBuffer = []
                if len(line)==1:
                    name = ''
                else:
                    name = line[1:]
            elif (line[0] == "+" and len(qualBuffer) == 0):
                qualFlag = True
            else:
                if (qualFlag):
                    qualBuffer.append( line )
                else:
                    seqBuffer.append( line )
        if len(seqBuffer)>0:
            yield FastqEntry( name, "".join(seqBuffer), "".join(qualBuffer) )

    def close( self ):
        self.infile.close()

def qvCharToInt( ch ):
    return ord(ch) - 33

def qvIntToChar( qv ):
    return chr(qv+33)

class FastqEntry:
    def __init__( self, name, sequence, quality ):
        self.name = name
        self.sequence = sequence
        self.quality = quality

    def __str__( self ):
        buffer = []
        buffer.append( "@" + self.name )
        buffer.append( self.sequence )
        buffer.append( "+" + self.name )
        buffer.append( self.quality )
        return os.linesep.join(buffer)
    
    def subseq(self, seqRange, name=None):
        if not name: 
            name = "%s_%i_%i" % (self.name, seqRange.start, seqRange.end)
        return FastqEntry(name, self.sequence[seqRange.start : seqRange.end], self.quality[seqRange.start : seqRange.end])

def prettyprint( sequence, width=70 ):
    return os.linesep.join( \
        [ sequence[i:i+width] for i in xrange(0,len(sequence),width) ] )

def writeFastaEntry( file, entry, width=70 ):
    file.write( '>%s%s%s%s' % ( entry.getTag(), os.linesep, \
                               prettyprint(entry.sequence,width=width), \
                               os.linesep ) )

def writeMsa( file, msa, fullFormat=False ):
    if msa.hasNames():
        names = msa.nameIter()
    else:
        names = imap( lambda x: str(x), count() )
    if fullFormat:
        for s, n in izip( iter(msa), names ):
            writeFastaEntry( file, FastaEntry( n, s.getFullString() ) )
    else:
        for e in msa.entryIter():
            header = '%s|start=%d|strand=%s' % ( e.name, \
                e.sequence.alignedRange.start, e.strand )
            seq = str(e.sequence)
            writeFastaEntry( file, FastaEntry( header, seq ) )
    #    for s, n in izip( iter(msa), names ):
    #        header = '%s|start=%d' % ( n, s.alignedRange.start )
    #        seq = str(s)
    #        writeFastaEntry( file, FastaEntry( header, seq ) )

class FastaSubseq:
    """ Command-line application for generating subsequences
        from reads files.
    """
    def __init__( self, argv ):
        self.__parseOptions( argv )

    def __parseOptions( self, argv ):
        usage = 'Usage: %prog [--help] [options] fastaFile > newFastaFile'
        parser = OptionParser( usage=usage )
        parser.add_option( '--prefix', type="int", \
            help='Number of bases to keep from the prefix of each sequence' )
        parser.add_option( '--suffix', type="int", \
            help='Number of bases to keep from the suffix of each sequence' )
        parser.add_option( '--windowstep', type="int", \
            help='How far to slide a window over for generating a set of subsequences' )
        parser.add_option( '--nwindows', type="int", \
            help='Number of windows to use when generating a set of subsequences (0 for whole sequence)' )
        
        parser.set_defaults( prefix=0, suffix=0, nwindows=1, windowstep=0 )

        self.options, args = parser.parse_args( argv )

        if len(args)!=2:
            parser.error( 'Expected 1 argument' )
        self.fastaName = args[1]

    def run( self ):
        reader = SimpleFastaReader( self.fastaName )
        for entry in reader:
            if self.options.prefix:
                start = 0
                end = self.options.prefix
                if self.options.nwindows==0:
                    nwindows = (len(entry.sequence)-self.options.prefix) \
                        / self.options.windowstep
                else:
                    nwindows = self.options.nwindows
                if nwindows<0:
                    continue
                for i in xrange(nwindows):
                    newSequence = entry.sequence[ start:end ]
                    print '>%d:%s\n%s' % ( start, \
                                           entry.name, prettyprint( newSequence ) )
                    start += self.options.windowstep
                    end += self.options.windowstep
                # print the last subsequence
                newSequence = entry.sequence[ start:len(entry.sequence) ]
                print '>%d:%s\n%s' % ( start, \
                                       entry.name, prettyprint(newSequence) )
                continue
            elif self.options.suffix:
                if len(entry.sequence)>self.options.suffix:
                    newSequence = entry.sequence[ -self.options.suffix: ]
            print '>%s\n%s' % ( entry.name, prettyprint( newSequence ) )
        reader.close()

def splitFasta(fasta, nSplits):
    """Splits fasta into nSplits files and return their absolute paths"""

    reader = SimpleFastaReader(fasta)
    n = int(nSplits)

    basePath = os.path.basename( fasta )
    fileNames = [ "%i.%s" % (i, basePath) for i in range(n) ]
    fhs = map(lambda x: open(x,"w"), fileNames)
    for idx, entry in enumerate(reader):
        fileIdx = idx % n
        print >>fhs[fileIdx], str( entry )
    reader.close()
    for fh in fhs: fh.close()
    splits = map(lambda x: os.path.abspath(x), fileNames)
    return splits

class ReaderFactory:

    def openFile(self, file):
        if file.endswith(".fstq") or file.endswith(".fastq"):
            return SimpleFastqReader(file)
        elif file.endswith(".fsta") or file.endswith(".fasta") or file.endswith(".fa"):
            return SimpleFastaReader(file)
        else:
            raise SystemExit, "Unsupported file type %s" % file

class Converter:

    def convertFastqFasta(self, fastqFile, fastaFile):
        fh = open(fastaFile, "w")
        for e in SimpleFastqReader(fastqFile): fh.write( "%s\n" % str( FastaEntry( e.name, e.sequence ) ) )
        fh.close()

# simple tests
if __name__=='__main__':
    if len(sys.argv)!=2:
        sys.stderr.write( "Usage: %s fastaFileName\n" % sys.argv[0] )
        sys.exit(1)

    # print prettyprint( "".join([ 'A' for i in xrange(100) ]) )

    fastaReader = FastaStreamReader( open(sys.argv[1]) )

    print [ entry.sequence for entry in fastaReader ]
