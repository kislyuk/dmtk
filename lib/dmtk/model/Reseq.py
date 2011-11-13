__doc__="""Classes for implementing a simulation model of resequencing."""
import os
import sys
import traceback
import numpy
import math
import random

from xml.sax import *
from xml.sax.handler import ContentHandler

def mapPtoX( cumDist, p, epsilon=1.0e-4 ):
    """binary search to find x which gives cumulative probability within
       epsilon. """
    x0 = cumDist.xmin
    x1 = cumDist.xmax
    x = (x0+x1)*0.5
    y = cumDist.value(x)
    while abs(y-p)>epsilon:
        if y>p:
            x1 = x
        else:
            x0 = x
        x = 0.5*(x0+x1)
        y = cumDist.value(x)
    return x

class SamplePrep:
    def __init__(self):
        self.samplePrepBias = None
        
    def setSamplePrepBias(self, samplePrepBias):
        self.samplePrepBias = samplePrepBias
        delta = 1.0e-5
        x = numpy.arange( 0, 1.0+delta, delta )
        y = numpy.repeat( 0.0, len(x) )
        y[0] = self.samplePrepBias.value(0.0) * delta
        for i, xv in enumerate(x[1:]):
            y[i+1] = self.samplePrepBias.value(xv) * delta + y[i]
        self.samplePrepCum = GriddedResponseCurve( 'samplePrepCumulative', 0.0, 1.0, delta )
        #self.samplePrepCum.x = x
        self.samplePrepCum.y = y
        
    def generateStartPoint(self, genomeLength):
        while 1:
            p = random.random()
            x = mapPtoX( self.samplePrepCum, p )
            yield int(round(x*genomeLength))

class Sequencer:
    def __init__(self):
        self.readLengthDistribution = None
        self.qualityValueProfile = None
        self.accuracyBias = None
        self.qualityValueDistribution = None
        
    def setReadLengthDistribution(self, readLengthDistribution):
        self.readLengthDistribution = readLengthDistribution
        self.readLengthCum = self.__calculateCumulative( self.readLengthDistribution )
        
    def setQualityValueProfile(self, qualityValueProfile):
        self.qualityValueProfile = qualityValueProfile
    
    def setQualityValueDistribution(self, qualityValueDistribution):
        self.qualityValueDistribution = qualityValueDistribution
        self.qualityValueCum = self.__calculateCumulative( self.qualityValueDistribution )
        
    def setAccuracyBias(self, accuracyBias):
        self.accuracyBias = accuracyBias
                                                  
    def __calculateCumulative(self, rc):
        n = 1.0e5
        d = (rc.xmax - rc.xmin) / n
        x = numpy.arange( rc.xmin, rc.xmax+d, d )
        y = numpy.repeat( 0.0, len(x) )
        y[0] = rc.value(rc.xmin)
        for i, xv in enumerate(x[1:]):
            y[i+1] = rc.value(xv) + y[i]
        y = y / y[-1]
        cumulative = GriddedResponseCurve( rc.name+"Cumulative", rc.xmin, rc.xmax, d )
        cumulative.y = y
        return cumulative
        
    def generateRead( self, startPoints, genomeLength ):
        """ A read from the sequencer is the following tuple:
        ( startPoint, length, qvs )
        where
        startPoint is the start point on the genome
        length is the read length
        qvs is an array of quality values
        """
        for sp in startPoints:
            p = random.random()
            rl = int(mapPtoX( self.readLengthCum, p ))
            p = random.random()
            qv = mapPtoX( self.qualityValueCum, p )
            bias = self.accuracyBias.value( float(sp) / float(genomeLength) )
            qvs = numpy.repeat( 0, rl )
            for i in xrange(rl):
                x = float(i)/float(rl)
                qvs[i] = int(round(qv * bias * self.qualityValueProfile.value(x)))
            yield sp, rl, qvs
            

class ReseqModelParameters:
    def __init__(self):
        self.genomeLength = 5000
        self.nUsableReads = 5000
        self.curves = {}
        
    def __str__(self):
        " For debugging output "
        buffer = []
        buffer.append( 'RMP {' )
        buffer.append( '  genomeLength = %d; ' % self.genomeLength )
        buffer.append( '  nUsableReads = %d; ' % self.nUsableReads )
        for i, curve in enumerate(self.curves.keys()):
            buffer.append( '  curve %d: %s' % ( i, curve ) )
            buffer.append( str(self.curves[curve]) )
        buffer.append( '}' )
        return os.linesep.join(buffer)
    
class ReseqModelParametersHandler( ContentHandler ):    
    def __init__( self ):
        self.parameters = None
        self.currentString = []
        self.currentCurve = None
      
    def startElement( self, name, attrs ):
        if name=='model':
            self.parameters = ReseqModelParameters()
            return
        if name=='genomeLength':
            self.parameters.genomeLength = int(attrs.getValue('value'))
        elif name=='nUsableReads':
            self.parameters.nUsableReads = int(attrs.getValue('value'))
        elif name=='curve':
            self.currentCurve = ResponseCurve( attrs.getValue('name') )
        elif self.currentCurve:
            if name=='xrange':
                self.currentCurve.setXrange( float(attrs.getValue('min')), float(attrs.getValue('max')) )
            elif name=='file':
                self.currentCurve.fromFile( attrs.getValue('name') )

    def endElement( self, name ):
        if self.currentCurve:
            if name=='expression':
                self.currentCurve.setExpression( self.getCurrentString() )
            elif name=='curve':
                self.parameters.curves[ self.currentCurve.name ] = self.currentCurve
                self.currentCurve = None
            elif name=='class':
                self.currentCurve = eval(self.getCurrentString())
        self.clearString()  

    def getCurrentString( self ):
        return "".join(self.currentString).strip()

    def clearString( self ):
        self.currentString = []

    def characters( self, ch ):
        self.currentString.append(ch)
    
class ResponseCurve:
    """ Models a response curve which is a function of a single variable.
    The curve can be stored as a python expression or modeled as a lookup table with
    linear interpolation.
    In general this class is rather slow.  Use one of the specialized subclasses
    when possible.
    """
    DEFAULT_N_FOR_OUTPUT = 500
    
    def __init__(self, name):
        self.name = name
        self.x = None
        self.xmin = -1.0
        self.xmax = -1.0
        self.y = None
        self.yexpr = None
        
    def setXrange(self, xmin, xmax ):
        self.xmin = xmin
        self.xmax = xmax
    
    def setExpression(self, yexpr):
        """ An expression which evaluates to the response of this function
        for a given value of x.
        
        Should be specified as a formula involving the unknown variable x.
        Functions should be python-callable functions.
        
        Example:
            x * 2 + 2
            1.0 - 0.5 * math.exp( -x / 10.0 )
        """
        # pre-compiling the expression sped up general use of this class by 10x
        self.yexpr = compile( yexpr, '<string>', 'eval' )    
    
    def value(self, x):
        if not self.yexpr:
            return numpy.interp( [x], self.x, self.y )
        return eval( self.yexpr )
    
    def fromFile(self, file):
        infile = open( file, 'r' )
        x, y = [], []
        for line in infile:
            values = line[:-1].split()
            x.append( float(values[0]) )
            y.append( float(values[1]) )
        infile.close()
        self.x = numpy.array( x )
        self.y = numpy.array( y )
        self.xmin = min(self.x)
        self.xmax = max(self.x)
    
    def __str__(self):
        buffer = [ "#  %s" % self.name, "#x\ty" ]
        if self.x is not None:
            for x in self.x:
                buffer.append( '%.2f\t%.3f' % ( x, self.value(x) ) )
        else:
            delta = ( self.xmax - self.xmin ) / ResponseCurve.DEFAULT_N_FOR_OUTPUT
            for x in numpy.arange( self.xmin, self.xmax+delta, delta ):
                buffer.append( '%.2f\t%.3f' % ( x, self.value(x) ) )
        return os.linesep.join(buffer)
    
class GriddedResponseCurve( ResponseCurve ):
    def __init__(self, name, xmin, xmax, xdel ):
        ResponseCurve.__init__(self, name)
        self.xmin = xmin
        self.xmax = xmax
        self.xdel = xdel
        
    def value(self, x):
        """ For gridded data this implementation of linear interpolation
        sped this lookup by 10x. """
        l = (x-self.xmin) / self.xdel
        if l < 0.0:
            return self.y[0]
        i0 = int(math.floor(l))
        if i0 >= len(self.y):
            return self.y[-1]
        d = l - float(i0)
        y = self.y[i0] * (1.0-d) + self.y[i0+1] * d
        return y
        
    def __str__(self):
        buffer = []
        buffer.append( '#  %s' % self.name )
        buffer.append( '#  xmin = %.2f' % self.xmin )
        buffer.append( '#  xmax = %.2f' % self.xmax )
        buffer.append( '#  xdel = %.2f' % self.xdel )
        return os.linesep.join(buffer)
    
class ConstantResponseCurve( ResponseCurve ):
    def __init__(self, name, value ):
        ResponseCurve.__init__(self, name)
        self._value = value
        
    def value(self, x):
        return self._value
    
    def __str__(self):
        return '#  %s\nconstant = %.2f' % ( self.name, self._value )
    
class GaussianResponseCurve( ResponseCurve ):
    def __init__(self, name, mu, sigma):
        ResponseCurve.__init__( self, name )
        self.mu = mu
        self.sigma = sigma
        self.norm = 1.0 / (math.sqrt(math.pi*2.0) * sigma)
        self.xmax = 10.0 * sigma + mu
        self.xmin = -10.0 * sigma + mu
        
    def value(self, x):
        dx = ( x - self.mu ) / self.sigma
        return self.norm * math.exp( - 0.5 * dx * dx )
    
    def __str__(self):
        return '#  %s\nmu = %.2f\nsigma = %.2f\nnorm = %.2f' % \
            ( self.mu, self.sigma, self.norm )