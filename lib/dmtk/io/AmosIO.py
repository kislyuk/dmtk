"""
Collection of routines for simplfying access to AMOS-formatted messages and data banks. 
Main entry points are the parseBankReport and iParseBankReport methods.
"""
import sys
import os
from tempfile import mkstemp
from dmtk.io.FastaIO import prettyprint
from dmtk.model.AlignedSequence import AlignedSequence
from dmtk.model.SeqUtils import revcomp
from dmtk.model.AlignmentHit import AlignmentHit 

# uses PATH for executable lookup
BANK_REPORT_PATH = 'bank-report'
BANK_TRANSACT_PATH = 'bank-transact'
BANK_UNLOCK_PATH = 'bank-unlock'

class AmosIOException( Exception ):
    """Generic exception thrown by the AmosIO API"""
    pass

class AmosMessageParser:
    """ Base class for a generic event parser for extracting AMOS messages.
        TODO: add more error handling, especially to subprocess calls.
        
        To use:  subclass and override startMessage, keyValuePair, and/or
        endMessage
        
        Parsing is activated by calling the parse() method. 
    """
    def __init__(self):
        self.buffer = []
        self.currentKey = None
        self.inMessage = False
        self.inField = False
        self.messageType = ""
        self.currentMessage = None

    def parseFile(self, fileName):
        infile = open( fileName, 'r' )
        for line in infile:
            self.__parseLine( line )
        infile.close()
        
    def __parseLine(self, line):
        line = line.rstrip('\n')
        if len(line)<1: return
        if line[0]=='{':
            self.startMessage( line[1:] )
        elif line[0]=='}':
            self.endMessage()
        # the 3:4 below allows line to be less than 4 chars long w/o error
        elif line[3:4] == ":"  and not self.inField: 
            if (len(line) == 4): # multiline field have only the key 
                self.inField = True
            else:
                self.inField = False # this field ends and begins here
            values = line.strip().split(':')
            if len(values)==2 and len(values[1])>0:
                self.keyValuePair( values[0].strip(), values[1].strip() )
            elif len(values)==2 or len(values)==1:
                self.currentKey = values[0].strip()
            else:
                return
        elif line[0]=='.':
            self.inField = False
            value = self.getDelimiter().join(self.buffer)
            self.keyValuePair( self.currentKey, value )
            self.buffer = []
            self.currentKey = None
        else:
            self.buffer.append( line )

    def parse(self, lines):
        """
          Very simple (and fragile) parser for event parsing an AMOS-formated
          message file.
          lines: iterable sequence of lines (with or without '\n')
          Returns iterator over messages.
        """
        first = True
        for l in lines:
            self.__parseLine( l )
            if not self.inMessage and not first:
                yield self.currentMessage
            first = False

    def startMessage(self, messageType):
        self.currentMessage = AmosGenericMessage(messageType)
        self.messageType = messageType
        self.inMessage = True
        self.inField = False
    
    def endMessage(self):
        self.inField = False
        self.inMessage = False
    
    def keyValuePair( self, key, value):
        if not self.currentMessage: return
        self.currentMessage[key] = value
    

    # for some multiline records it makes sense to have a line delimiter (i.e. gap locations)
    # for others it doesn't (i.e. sequence data)
    def getDelimiter(self):
        if self.currentKey == "gap":
            return " "
        else:
            return ""
    
    def addListener(self, listener):
        self.listeners.append( listener )

class AmosTestParser( AmosMessageParser ):
    """ Test parser which outputs state to stdout as it parses. """
    def __init__(self):
        AmosMessageParser.__init__( self )
        self.currentMessageType = None
        
    def startMessage(self, messageType):
        print 'Starting %s' % messageType
        self.currentMessageType = messageType
        
    def endMessage(self):
        print 'Ending %s' % self.currentMessageType
        
    def keyValuePair(self, key, value):
        print 'At %s = %s' % ( key, value ) 

class AmosReadListener:
    """ Interface for classes which are interested in reads as they
      become available from the AmosReadParser """
    def processRead(self, read):
        " Subclasses override to receive read events "
        pass

class AmosContigListener:
    """ Interface for classes which are interested in reads as they
      become available from the AmosContigParser """
    def processContig(self, read):
        " Subclasses override to receive contig events "
        pass

class AmosSequenceAccumulator( AmosReadListener ):
    def __init__(self):
        self.id2seq =  {}

    def processRead(self, read):
        " Subclasses override to receive read events "
        self.id2seq[str(read.iid)] = read.seq

    def getSeqsByIID(self):
        return self.id2seq

class AmosReadAccumulator( AmosReadListener ):
    def __init__(self):
        self.reads = []
        
    def getReads(self):
        return self.reads
        
    def processRead(self, read):
        self.reads.append( read )

class AmosContigAccumulator( AmosContigListener ):
    def __init__(self):
        self.contigs = []
        
    def getContigs(self):
        return self.contigs
        
    def processContig(self, contig):
        self.contigs.append( contig )

class AmosLayoutAccumulator( ):
    def __init__(self):
        self.layouts = []
        
    def getLayouts(self):
        return self.layouts
        
    def processLayout(self, layout):
        self.layouts.append( layout )
        
class AmosLayoutParser( AmosMessageParser ):
    def __init__(self):
        AmosMessageParser.__init__(self)
        self.nestLevel = 0
        self.layoutAccumulator = None
        self.listeners = []
    
    def getLayouts(self):
        return self.layoutAccumulator.getLayouts()
    
    def setAccumulateLayouts(self, accumulateLayouts ):
        if accumulateLayouts:
            self.layoutAccumulator = AmosLayoutAccumulator()
            self.listeners.append( self.layoutAccumulator )
    
    def startMessage(self, messageType):
        if messageType=='LAY': 
            self.inMessage = True
            self.currentMessage = AmosLayoutMessage()
        else:
            assert( self.inMessage )
            self.currentMessage.children.append( AmosTileMessage() )
        if self.inMessage:
            self.nestLevel += 1

    def endMessage(self):
        if self.inMessage:
            self.nestLevel -= 1
        if self.nestLevel==0:
            if self.currentMessage:
                for l in self.listeners:
                    l.processLayout( self.currentMessage )
            self.inMessage = False
            
    def keyValuePair(self, key, value):
        
        if not self.currentMessage:
            return
        if self.nestLevel == 1:
            if key=='iid':
                self.currentMessage['iid'] = value
        elif self.nestLevel == 2: 
            currentChild = self.currentMessage.children[-1]
            if key == 'gap':
                currentChild['gap'] = value
            elif key == 'clr':
                currentChild['clr'] = value
            elif key == 'off':
                currentChild['off'] = value
            elif key == 'src':    
                currentChild['src'] = value        
        else:
            raise SystemExit, "Unsupported nesting level!"

class AmosContigParser( AmosMessageParser ):

    def __init__(self):
        AmosMessageParser.__init__(self)
        self.inMessage = False
        self.nestLevel = 0
        self.contigsAccumulator = None
        self.listeners = []
    
    def getContigs(self):
        return self.contigsAccumulator.getContigs()
    
    def setAccumulateContigs(self, accumulateContigs ):
        if accumulateContigs:
            self.contigsAccumulator = AmosContigAccumulator()
            self.listeners.append( self.contigsAccumulator )
    
    def startMessage(self, messageType):
        if messageType=='CTG': 
            self.inMessage = True
            self.currentMessage = AmosContigMessage()
        else:
            assert( self.inMessage )
            self.currentMessage.children.append( AmosTileMessage() )
        if self.inMessage:
            self.nestLevel += 1

    def endMessage(self):
        if self.inMessage:
            self.nestLevel -= 1
        if self.nestLevel==0:
            if self.currentMessage:
                for l in self.listeners:
                    l.processContig( self.currentMessage )
            self.inMessage = False
            
    def keyValuePair(self, key, value):
        
        if not self.currentMessage:
            return
        if self.nestLevel == 1:
            if key=='iid':
                self.currentMessage['iid'] = value
            elif key=='eid':
                self.currentMessage['eid'] = value
            elif key=='seq':
                self.currentMessage['seq'] = value
            elif key=='qlt':
                self.currentMessage['qlt'] = value
        elif self.nestLevel == 2: 
            currentChild = self.currentMessage.children[-1]
            if key == 'gap':
                currentChild['gap'] = value
            elif key == 'clr':
                currentChild['clr'] = value
            elif key == 'off':
                currentChild['off'] = value
            elif key == 'src':    
                currentChild['src'] = value        
        else:
            raise SystemExit, "Unsupported nesting level!"

class AmosContigLightParser( AmosContigParser ):
    """Ignores the most memory intensive data fields (seq, qlt and gap)"""

    def __init__(self):
        AmosContigParser.__init__(self)

    def keyValuePair(self, key, value):
        
        if not self.currentMessage:
            return
        if self.nestLevel == 1:
            if key=='iid':
                self.currentMessage['iid'] = value
            elif key=='eid':
                self.currentMessage['eid'] = value
            elif key=='seq':
                pass
            elif key=='qlt':
                pass
        elif self.nestLevel == 2: 
            currentChild = self.currentMessage.children[-1]
            if key == 'gap':
                pass
            elif key == 'clr':
                currentChild['clr'] = value
            elif key == 'off':
                currentChild['off'] = value
            elif key == 'src':    
                currentChild['src'] = value        
        else:
            raise SystemExit, "Unsupported nesting level!"

class AmosReadParser( AmosMessageParser ):
    def __init__(self):
        AmosMessageParser.__init__(self)
        self.readsAccumulator = None
        self.listeners = []
    
    def getReads(self):
        return self.readsAccumulator.getReads()
    
    def setAccumulateReads(self, accumulateReads ):
        if accumulateReads:
            self.readsAccumulator = AmosReadAccumulator()
            self.listeners.append( self.readsAccumulator )
    
    def startMessage(self, messageType):
        if messageType!='RED':
            return
        AmosMessageParser.startMessage(self, messageType)
        self.currentMessage = AmosReadMessage()
        
    def endMessage(self):
        if self.currentMessage:
            for l in self.listeners:
                l.processRead( self.currentMessage )
        AmosMessageParser.endMessage(self)
        # self.currentMessage = None
        
    def keyValuePair(self, key, value):
        if not self.currentMessage:
            return
        if key=='iid':
            self.currentMessage.iid = int(value)
        elif key=='eid':
            self.currentMessage.eid = value
        elif key=='seq':
            value = value.rstrip("\0").lstrip("\0")
            self.currentMessage.seq = value
        elif key=='qlt':
            self.currentMessage.qlt = value
        elif key=='frg':
            self.currentMessage.frg = int(value)
        elif key=='clr':
            self.currentMessage.clr = value

class AmosGenericMessage:
    """Allows any type of AMOS message to be modeled."""
    def __init__(self, messageType, properties=None):
        self.messageType = messageType
        self.properties = {}
        if properties:
            self.properties.update( properties )
        self.children = []
        self.lineLength = 70
        
    def __setitem__(self, key, value):    
        self.setProperty(key, value)
        
    def __getitem__(self, key):
        return self.getProperty(key)
        
    def setProperty(self, key, value):
        self.properties[ key ] = str(value)
        
    def getProperty(self, key):
        return self.properties[key]
        
    def addChild(self, child):
        self.children.append( child )
        
    def format(self, value):
        if len(str(value)) < self.lineLength: return value
        # if the value contains newlines assume it's already formatted
        if os.linesep in value:
            return value
        buffer = []
        for i in xrange( 0, len(value), self.lineLength ):
            buffer.append( value[ i : i+self.lineLength ] )
        buffer.append('.')
        return os.linesep + os.linesep.join(buffer)
        
    def __str__(self):
        """Note that the string representation does *not* have a newline at the end. To make it
        a valid message, it should have a newline appended, as in default print statements"""
        buffer = []
        buffer.append( '{%s' % self.messageType )
        for k,v in self.properties.iteritems():
                
            if (v != '' and v != None):
                if k == "eid": # eid is special seems to break if it goes more than one line
                    buffer.append( '%s:%s' % (k,v) )
                else:
                    buffer.append( '%s:%s' % (k,self.format(v)) )
        for child in self.children:
            buffer.append( str(child) )
        buffer.append('}')
        return os.linesep.join( buffer )
    
class AmosReadMessage(AmosGenericMessage):
    """ I/O storage class for representing an AMOS RED message """
    def __init__(self):
        AmosGenericMessage.__init__(self, "RED")
        self.iid = 0
        self.eid = ''
        self.seq = ''
        self.qlt = ''
        self.frg = 0
        self.clr = ''
        
    def __str__(self): 
        """For backwards compatibility. TODO should eventually be replaced by 
        AmosGenericMessage.__str__()"""
        buffer = [ '{RED' ]
        buffer.append( 'iid:%d' % self.iid )
        if (self.eid != ''): # tigger doesn't like a null eid field
            buffer.append( 'eid:%s' % self.eid )
        buffer.append( 'seq:\n%s\n.' % prettyprint( self.seq ) ) 
        buffer.append( 'qlt:\n%s\n.' % prettyprint( self.qlt ) )
        buffer.append( 'frg:%d' % self.frg )
        buffer.append( 'clr:%s' % self.clr )
        buffer.append( '}' )
        return '\n'.join(buffer)
    
    def toFasta(self):
        lineLength = 60
        id = self.eid 
        if id == '':
            id = self.iid
        buffer = ">%s\n" % id
        buffer += prettyprint(self.seq, lineLength)
        return buffer

class AmosContigMessage(AmosGenericMessage):
    """ I/O storage class for representing an AMOS CTG message """

    def __init__(self):
        AmosGenericMessage.__init__(self," CTG")
        self.properties['iid'] = 0
        self.properties['eid'] = ''
        self.properties['seq'] = ''
        self.properties['qlt'] = ''
        self.messageType = "CTG"

    def toAmosReadMessage(self):
        read = AmosReadMessage()
        read.iid = self.properties['iid']
        read.eid = "%s-contig" % self.properties['iid'] 
        read.frg = 0

        # amos contig seqs will sometimes contain gap characters
        # but these break when used in reads. unfortunately, the 
        # easiest way to deal with them is to remove them and their
        # associated quality scores, even though this results in 
        # information loss.
        seqlist = list( self.properties['seq'] )
        qltlist = list( self.properties['qlt'] )
        assert(len(seqlist) == len(qltlist))
        newseq = []
        newqlt = []
        for idx in range(len(seqlist)):
            if (seqlist[idx] != '-'):
                newseq.append(seqlist[idx])
                newqlt.append(qltlist[idx])
        read.seq = ''.join(newseq)
        read.qlt = ''.join(newqlt)
        read.clr = "0,%i" % len(read.seq)

        return read

class AmosOverlapMessage(AmosGenericMessage):
    """ I/O storage class for representing an AMOS OVL message """

    def __init__(self):
        AmosGenericMessage.__init__(self, "OVL")
        self.properties = {}
        self.properties['ahg'] = ''
        self.properties['bhg'] = ''
        self.properties['rds'] = ''
        self.properties['adj'] = ''
        self.properties['scr'] = ''
        self.messageType = "OVL"

    def fromHit(self, hit, bank):
        """Needs a bank to get the iid"""
        assert(hit.target_strand == "+")
        try:
            qIID = bank.getReadIIDfromEID(hit.query_id)
        except AmosIOException, e:
            print >>sys.stderr, "Could not find %s" % hit.query_id
            raise e
        try:
            tIID = bank.getReadIIDfromEID(hit.target_id)
        except AmosIOException, e:
            print >>sys.stderr, "Could not find %s" % hit.target_id
            raise e
        # invert the sense of the overlap because target is always positive
        self.properties['rds'] = "%i,%i" % (tIID, qIID)
        if hit.query_strand == "+":
            self.properties['ahg'] = hit.target_start - hit.query_start
            self.properties['bhg'] = (hit.query_length  - hit.query_end) - \
                                     (hit.target_length - hit.target_end)
        else:
            self.properties['ahg'] = hit.target_start -  \
                                     (hit.query_length  - hit.query_end) 
            self.properties['bhg'] = hit.query_start - \
                                     (hit.target_length - hit.target_end)

        self.properties['uid1'] = tIID
        self.properties['uid2'] = qIID
        self.properties['scr'] = int(hit.zScore * 100)
        self.properties['adj'] = "N" if hit.query_strand == "+" else "I"
    
class AmosLayoutMessage(AmosGenericMessage):
    """ I/O storage class for representing an AMOS LAY message """

    def __init__(self):
        AmosGenericMessage.__init__(self, "LAY")
        self.properties = {}
        self.properties['iid'] = ''
        self.messageType = "LAY"

    def tileOffsetZero(self):
        
        offsets = map(lambda x: int(x['off']), self.children)
        offsets.sort()
        adjustment = offsets[0]
        for child in self.children:
            child['off'] = str( int(child['off']) - adjustment )
        
class AmosTileMessage(AmosGenericMessage):
    """ I/O storage class for representing an AMOS TLE message """

    def __init__(self):
        AmosGenericMessage.__init__(self, "TLE")
        self.properties = {}
        self.properties['clr'] = ''
        self.properties['gap'] = ''
        self.properties['off'] = ''
        self.properties['src'] = ''
        self.messageType = "TLE"

    def length(self):
        return max( map(int, self.properties['clr'].split(",") ))

    def gaps(self):
        if self.properties['gap'] == '': return []
        return map(int, self['gap'].split(' '))

    def toAlignmentHit(self, querySequence, contig):
        """Convert this AmosTileMessage to an AlignmentHit object. 
           The message does not contain the read or target sequences, 
           so we need to pass these in. Note this currently does not maintain the target
           sequence."""

        alignedSeq = self.toAlignedSequence(querySequence)

        hit = AlignmentHit()

        hit.query_id = alignedSeq.seqName
        hit.target_id = "%s_contig" % contig['iid']

        clearRange = map(int, self['clr'].split(','))
        hit.query_start = min(clearRange)
        hit.query_end = max(clearRange)

        hit.target_start = alignedSeq.getAlignedStart()
        hit.target_end   = alignedSeq.getAlignedEnd()

        alignedQuery = str(alignedSeq)
        contigAlignedSeq = "A" * (hit.target_end - hit.target_start)

        if clearRange[0] < clearRange[1]:
            hit.target_strand = "+" 
            hit.alignedQuery = alignedQuery
            hit.alignedTarget = contigAlignedSeq
        else:
            hit.target_strand = "-"
            hit.alignedQuery = revcomp(str(alignedQuery))
            hit.alignedTarget = revcomp(contigAlignedSeq)
        hit.query_strand = "+" # for HDF5, query strand is always positive

        return hit

    def toAlignedSequence(self, sequence):
        """Convert this AmosTileMessage to an AlignedSequence object. 
           The message does not have the read sequence, so we need to pass it in."""

        clr = map(int, self['clr'].split(','))
        end = max(clr[0], clr[1])

        if clr[1] < clr[0]: 
            # we need to reverse complement the sequence
            sequence = revcomp( sequence )

        alignedSeq = AlignedSequence( sequence )
        off = int(self['off'])
        alignedSeq.setAlignedStart( off  )

        alignedSeq.setAlignedEnd( off + end ) 
        alignedSeq.seqName = self['src']

        gaps = []
        if (self['gap'] != ""):
            gaps = map(lambda x: int(x) - 1, self['gap'].split(' '))
            alignedSeq.insertGaps(gaps, offset=off )
        alignedSeq.setAlignedEnd( alignedSeq.getAlignedEnd() + len(gaps) )

        return alignedSeq

# Static methods for parsing an existing bank report 
def parseBankReport( bankPath, objectType, parser=AmosTestParser() ):
    """ Applies the given parser to the report generated by bank-report. Returns a list of messages."""

    return [ p for p in iParseBankReport( bankPath, objectType, parser ) ]

def iParseBankReport( bankPath, objectType, parser=AmosTestParser() ):
    """ Applies the given parser to the report generated by bank-report. Returns an iterator over messages."""
    handle, tempOut = mkstemp(".msg")
    os.close(handle)
    cmdLine = '%s -b %s %s > %s' % (BANK_REPORT_PATH, bankPath, objectType, tempOut)
    
    retVal = os.system(cmdLine)
    if retVal: raise AmosIOException, "ERROR: %s failed" % cmdLine
    
    idx = 0
    out = open(tempOut)
    for p in parser.parse( out ): 
        idx += 1
        # if idx % 10000 == 0:
        #    print >>sys.stderr, "Reached record %i" % idx
        yield p
    os.unlink(tempOut)

def amosQvCharToInt( ch ):
    """Static method to convert AMOS quality values to integer representation. 
    Note that it is different than the standard FASTQ which uses 33 instead of 48."""
    return ord(ch) - 48

if __name__=='__main__':
    pass
