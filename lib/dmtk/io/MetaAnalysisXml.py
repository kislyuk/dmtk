#!/usr/bin/env python
__doc__='''Data model and I/O classes to support the meta-analysis
XML format.'''
import sys
import os
import glob
import copy

from xml.etree.cElementTree import *

from dmtk.io.PacBioIO import *

# Default primary analysis protocol if not specified (should always be
# specified for clarity)
DEFAULT_PROTOCOL = 'DragonRabbit'
#
# Names of tag types that the resolver knows about
#
TAG_MOVIE_TYPE = 'movieType'
TAG_PRIMARY_ANALYSIS_PROTOCOL = 'primaryAnalysisProtocol'
TAG_FILE_TYPE = 'fileType'

class IInputDataUrl:
    """Base interface for class which model data which is present
    at a location (URL)."""
    def toLocal( self, multiplicity=0 ):
        """If multiplicity is set to 1 then returns a single file name
        if the URL resolves to a single file or raises an IOError.

        Otherwise returns a list of files which result from this URL lookup."""
        pass

    def catAll( self, newFileName, data_filter=lambda x:True ):
        """Concatenates the files represented by this input object
        into a single text file.  The optional filter function is
        a function which operates on each data item and returns
        True if the item should be kept in the concatenated output."""
        pass

    def isGlobal( self ):
        """Indicates whether or not this input must be copied to 
        a shared location for processing on compute nodes."""
        # Err on the side of caution
        return False

class MetaAnalysisInputs( IInputDataUrl ):
    """Models the meta-analysis input data model
    See http://usmp-acid/xml/pb_meta_analysis.xsd.
    
    Instantiate from an XML file using the file constructor or load()."""
    def __init__( self, file=None, resolver=None ):
        self.dataReferences = MetaAnalysisReferences( self )
        self.header = None
        self.metadata = None
        if resolver:
            self.resolver = resolver
        else:
            self.resolver = PacBioUrlResolver()
        if file:
            self.load(file)

    def load( self, file ):
        """Instantiates this object from an XML representation.  Supports
        the older smrtpipeInputs format for backwards compatibility."""
        tree = ElementTree( file=file )
        rootNode = tree.getroot()
        if rootNode is None or rootNode.tag!='pacbioAnalysisInputs':
            if rootNode.tag=='smrtpipeInputs':
                self.__loadOldSmrtpipeFormat( rootNode )
                return
            raise IOError( "file %s seems to have an "\
                "incorrect root node; rootNode=%s" % (file,str(rootNode)) )
        headerNode = rootNode.find( 'header' )
        if headerNode is not None:
            self.header = MetaAnalysisHeader( headerNode )
        dataReferencesNode = rootNode.find( 'dataReferences' )
        self.dataReferences.populate( dataReferencesNode )
        metadataNode = rootNode.find( 'metadata' )
        if metadataNode is not None:
            self.metadata = MetaAnalysisMetadata( metadataNode )
            self.__assignAllConditions()

    def __loadOldSmrtpipeFormat( self, rootNode ):
        jobNode = rootNode.find('job')
        if jobNode is not None:
            self.header = MetaAnalysisHeader()
            self.header.id = jobNode.attrib['id']
            dNode = jobNode.find('description')
            if dNode is not None:
                nameNode = dNode.find('name')
                if nameNode is not None:
                    self.header.name = nameNode.text.strip() 
                commentNode = dNode.find('comment')
                if not commentNode is None:
                    if not commentNode.text is None:
                        self.header.description = commentNode.text.strip()
        urlList = rootNode.findall( 'urls/url' )                  
        if urlList:                                         
            for urlNode in urlList:                           
                url = urlNode.attrib['ref']                   
                self.dataReferences.add( url )

    def __assignAllConditions( self ):
        for assignment in self.metadata.conditionMap:
            assignment.apply( self.dataReferences )

    def __str__( self ):
        return self.toXml()

    def __iter__( self ):
        return iter( self.dataReferences )

    def __len__( self ):
        return len( self.dataReferences )

    def getDataReferences( self ):
        return self.dataReferences

    def save( self, fileName ):
        outfile = open( fileName, 'w' )
        print >>outfile, '<?xml version="1.0"?>'
        print >>outfile, self.toXml()
        outfile.close()

    def toXml( self ):
        buffer = []
        buffer.append( '<pacbioAnalysisInputs>' )
        if self.header is not None:
            buffer.append( self.header.toXml() )
        buffer.append( self.dataReferences.toXml() )
        if self.metadata is not None:
            buffer.append( self.metadata.toXml() )
        buffer.append( '</pacbioAnalysisInputs>' )
        return os.linesep.join( buffer )

    def toLocal(self):
        raise SystemExit, "toLocal() method not implemented for %s" % self.__class__

    def catAll( self, newFileName, data_filter=lambda x:True ):
        """Concatenates the files represented by this input object
        into a single text file.  Might want to put logic in 
        later for concatenating XML files (if we have XML inputs)
        or HDF files."""
        outfile = open(newFileName,'w')
        for url in filter(data_filter, self.dataReferences):
            localFiles = url.toLocal()
            for localFile in localFiles:
                infile = open( localFile, 'r' )
                for line in infile:
                    outfile.write( line )
                infile.close()
        outfile.close()

    def getUrlsForCondition( self, condition ):
        pass

    def isGlobal( self ):
        """At the moment, we can assume that the XML specified
        data are global in nature (accessible from any node)"""
        return True

class MetaAnalysisReferences:
    """Models the set of data inputs in this meta-analysis.
    The data inputs are stored as InputDataUrl.
    """
    def __init__( self, parent=None ):
        self.urls = {}
        self.tags = []
        self.parent = parent

    def populate( self, node ):
        tags = node.findall( 'tag' )
        for tagNode in tags:
            type = tagNode.attrib['type']
            value = tagNode.text.strip()
            self.parent.resolver.setTag( type, value )
            self.tags.append( (type,value) )
        urls = node.findall( 'url' )
        if urls is None: urls=[]
        #
        # 12.02.09 - JMS: support newer format
        #
        urls2 = node.findall( 'data' )
        if urls2 is not None:
            urls.extend(urls2)
        for urlNode in urls:
            self._addUrlNode( urlNode )

    def _addUrlNode( self, urlNode ):
        ref = urlNode.attrib['ref']
        if 'rel' in urlNode.attrib:
            rel = urlNode.attrib['rel']
        else:
            rel = None
        resolvedPaths = []
        locations = urlNode.findall('location')
        for location in locations:
            resolvedPaths.append( location.text.strip() )
        self.add( ref, rel, resolvedPaths )

    def add( self, urlRef, urlRel=None, resolvedPaths=[] ):
        self.urls[ urlRef ] = InputDataUrl( urlRef, urlRel, \
            self.parent.resolver, resolvedPaths )

    def __iter__( self ):
        return iter( self.urls.values() )

    def __len__( self ):
        return len( self.urls )

    def toXml( self ):
        buffer = []
        buffer.append( '<dataReferences>' )
        for k,v in self.tags:
            buffer.append( '<tag type="%s">%s</tag>' % ( k, v ) )
        for v in self.urls.values():
            buffer.append( v.toXml() )
        buffer.append( '</dataReferences>' )
        return os.linesep.join(buffer)

class MetaAnalysisHeader:
    """Models the header information regarding the input
    to this meta-analysis."""
    def __init__( self, node=None ):
        self.id = None 
        self.name = ''
        self.description = None
        if node is not None:
            self.__populate( node )

    def __populate( self, node ):
        self.id = node.attrib[ 'id' ]
        self.name = node.attrib[ 'name' ]
        dNode = node.find( 'description' )
        if dNode is not None:
            self.description = dNode.text

    def toXml( self ):
        buffer = []
        if self.id:
            buffer.append( '<header id="%s" name="%s">' % \
                ( self.id, self.name ) )
        else:
            buffer.append( '<header name="%s">' % self.name )
        if self.description:
            buffer.append( '<description>' )
            buffer.append( self.description )
            buffer.append( '</description>' )
        buffer.append( '</header>' )
        return os.linesep.join(buffer)

class MetaAnalysisMetadata:
    """Models the meta-data (annotations) on the input to this
    meta-analysis."""
    def __init__( self, node=None ):
        # associates ConditionGroup name with ConditionGroup
        self.conditionGroups = {}
        # list of ConditionAssignment
        self.conditionMap = []
        if node is not None:
            self.__populate( node )

    def __populate( self, node ):
        groups = node.findall( 'conditionGroup' )
        for groupNode in groups:
            name = groupNode.attrib['name']
            self.conditionGroups[ name ] = ConditionGroup( name )
            conditions = groupNode.findall( 'condition' )
            for conditionNode in conditions:
                id = conditionNode.attrib['id']
                cName = conditionNode.attrib['name']
                group = self.conditionGroups[name]
                group.add( id, Condition(id,cName,group) )
        mapNode = node.find( 'conditionDataMap' )
        if mapNode is not None:
            cdNodes = mapNode.findall( 'conditionData' )
            for cdNode in cdNodes:
                ref = cdNode.attrib[ 'conditionRef' ]
                data = cdNode.attrib[ 'dataRef' ]
                assignment = ConditionAssignment( ref, url=data )
                self.conditionMap.append( assignment )
            crNodes = mapNode.findall( 'conditionRegEx' )
            for crNode in crNodes:
                ref = crNode.attrib[ 'conditionRef' ]
                regex = crNode.text
                self.conditionMap.append( \
                    ConditionAssignment( ref, regex=regex ) )
            # look up the IDs in the condition assignments and resolve
            # these to actual conditions in the data model
            for assignment in self.conditionMap:
                assignment.resolveCondition( self.conditionGroups.values() )

    def toXml( self ):
        buffer = []
        buffer.append( '<metadata>' )
        for group in self.conditionGroups.values():
            buffer.append( group.toXml() )
        buffer.append( '<conditionDataMap>' )
        for cdAssign in self.conditionMap:
            buffer.append( cdAssign.toXml() )
        buffer.append( '</conditionDataMap>' )
        buffer.append( '</metadata>' )
        return os.linesep.join(buffer)

class ConditionGroup:
    """Models a group of related and mutually exclusive conditions.
    
    For example, 'enzyme' might be a condition group---you can have
    one enzyme assigned or another enzyme assigned, but not both.  Or
    if both enyzmes is a valid condition then this is explicitly assigned
    as a separate condition."""
    def __init__( self, name ):
        self.name = name
        self.conditions = {}

    def add( self, key, condition ):
        self.conditions[key] = condition

    def contains( self, condition ):
        return condition.id in self.conditions

    def __iter__( self ):
        return iter( self.conditions.values() )

    def toXml( self ):
        buffer = []
        buffer.append( '<conditionGroup name="%s">' % self.name )
        for v in self.conditions.values():
            buffer.append( v.toXml() )
        buffer.append( '</conditionGroup>' )
        return os.linesep.join( buffer )

class Condition:
    """Models a condition under which a particular set of data
    (referenced by URL) was generated."""
    def __init__( self, id, name, group ):
        self.id = id
        self.name = name
        self.description = ''
        self.group = group

    def __str__( self ):
        "For debugging"
        return 'C { id=%s; name=%s; description=%s; group=%s}' % \
            ( self.id, self.name, self.description, self.group.name )

    def toXml( self ):
        return '<condition id="%s" name="%s">%s</condition>' % \
            ( self.id, self.name, self.description )

class ConditionAssignment:
    """Models the association of a condition with input data."""
    def __init__( self, conditionRef, regex=None, url=None ):
        self.conditionRef = conditionRef
        self.regex = regex
        self.url = url
        self.condition = None

    def apply( self, references ):
        if self.url is not None:
            for url in filter( lambda r: r.ref==self.url, references ):
                url.addCondition( self.condition )
        if self.regex is not None:
            matcher = re.compile( self.regex )
            for url in filter( lambda r: matcher.search(r.ref), references ):
                url.addCondition( self.condition )

    def resolveCondition( self, conditionGroups ):
        for group in conditionGroups:
            for condition in group:
                if condition.id==self.conditionRef:
                    self.condition = condition
                    return

    def toXml( self ):
        if self.regex:
            return '<conditionRegEx conditionRef="%s">%s</conditionRegEx>' %\
                ( self.condition, self.regex )
        if self.url:
            return '<conditionData conditionRef="%s" dataRef="%s">' % \
                ( self.condition, self.url )
        return ''


class InputDataUrl( IInputDataUrl ):
    """Models a URL referring to the location of some data intended for
    input.

    Unlike the conventional notion of a e.g. http URL,
    when resolved to a local file the URL could resolve to multiple files.  
    For example, the URL run:1920005-0004 might resolve to three separate
    FASTA files.

    The heart of the magic is the toLocal() method and the underlying
    IMetaAnalysisUrlResolver.
    """
    def __init__( self, ref, rel, resolver, resolvedPaths=[] ):
        self.ref = ref
        self.rel = rel
        self.resolver = resolver
        self.scheme = None
        self.path = ''
        self.conditions = []
        self.resolvedPaths = resolvedPaths
        self.__parseRef()

    def __parseRef( self ):
        values = self.ref.split(':')
        if len(values)==2:
            self.scheme = values[0]
            self.path = values[1]
        elif len(values)==1:
            self.scheme = 'file'
            self.path = values[0]

    def addCondition( self, condition ):
        # check if a condition from this group has already been
        # added
        for c in self.conditions:
            if c.group.contains( c ):
                raise IOError( \
                    'Trying to assign conflicting conditions to %s' \
                    % str(self) )
        self.conditions.append( condition )

    def hasCondition( self, condition ):
        return condition in self.conditions

    def __str__( self ):
        return '%s/%s' % ( self.scheme, self.path )

    def toLocal( self, multiplicity=0 ):
        """If multiplicity is set to 1 then returns a single file name
        if the URL resolves to a single file or raises an IOError.

        Otherwise returns a list of files which result from this URL lookup."""
        if len(self.resolvedPaths)!=0:
            return self.resolvedPaths
        localPaths = self.resolver.resolve( self )
        if len(localPaths)==0:
            print >>sys.stderr, "Can't resolve URL %s" % str(self)
            print >>sys.stderr, self.resolver.getLastError()
        if multiplicity==1:
            if len(localPaths)!=1:
                raise IOError( "URL %s doesn't resolve to a single file (resolves to %d)" % ( str(self), len(localPaths) ) )
            return localPaths[0]
        return localPaths

    def catAll( self, newFileName, data_filter=lambda x:True ):
        """Concatenates the files represented by this input object
        into a single text file."""

        # if this object does not obey the filter, then do nothing 
        if not data_filter(self): return 

        outfile = open(newFileName,'w')
        for localFile in self.toLocal():
            infile = open( localFile, 'r' )
            for line in infile:
                outfile.write( line )
            infile.close()
        outfile.close()

    def toXml( self ):
        buffer = [ '<data ' ]
        if self.rel:
            buffer.append( 'rel="%s" ' % self.rel )
        buffer.append( 'ref="%s">' % self.ref )
        if len(self.resolvedPaths)!=0:
            buffer.extend( [ '<location>%s</location>' % p for p in self.resolvedPaths ] )
        buffer.append( '</data>' )
        return "".join(buffer)

    def maxChunks( self, context ):
        """Returns the maximum number of chunks we can provide."""
        return 1

    def getChunk( self, chunk ):
        """Returns each file resolvable by toLocal() as a chunk."""
        return self

    def getUrl( self ):
        """Returns a URL constructed from the scheme and path."""
        return self.scheme + ":" + self.path

    def pathReplace( self, find, replace ):
        """Replaces all occurrences of 'find' with 'replace' in the path"""
        self.path = self.path.replace( find, replace )

class FofnUrl( InputDataUrl ):
    """Provides additional chunking power based on the knowledge that
    this is a .fofn file."""
    
    def maxChunks( self, context ):
        return len( [ fn for fn in open( self.toLocal(1), 'r' ) if len(fn.strip()) > 0 ] )

    def getChunk( self, chunk ):
        """Overwrites the existing FOFN with a chunked FOFN, returns self"""
        iChunk, nChunks = map( int, chunk.strip( ).split("/") )
        myFiles = [ fn.strip( ) for i,fn in enumerate( open( self.toLocal(1), 'r' ) ) if i % nChunks == iChunk ] 
        out = open( self.toLocal(1), 'w' )
        out.write( "\n".join( myFiles ) )
        out.close()
        return self
        

class IMetaAnalysisUrlResolver:
    """Interface for classes which resolve input URLs to files."""
    def resolve( self, inputDataUrl ):
        """Returns a sequence of paths on the local file system that
        correspond to the given URL."""
        pass

class PacBioUrlResolver( IMetaAnalysisUrlResolver ):
    """Resolve input URLs to files for the internal PacBio infrastructure."""
    def __init__( self ):
        # implemented as a multi-map: each key has a list of values
        self.tags = {}
        setRootPath( '/mnt' )
        self._error = ''
        self._paths = ['.']

    def setTag( self, key, value ):
        if key not in self.tags:
            self.tags[ key ] = []
        if value not in self.tags[key]:
            self.tags[ key ].append(value)

    def prependPath( self, path ):
        "Prepends a path to the search path list for file resolution."
        self._paths.insert(0,path)

    def appendPath( self, path ):
        "Appends a path to the search path list for file resolution."
        self._paths.append(path)

    def clearTags( self ):
        self.tags = {}

    def getTags( self ):
        return self.tags

    def getTag( self, tag ):
        if tag in self.tags:
            return self.tags[ tag ]
        return None

    def hasTag( self, tag ):
        return tag in self.tags

    def describeTags( self ):
        """Returns a descriptive string describing the tags that are
        set for this resolver."""
        if len(self.tags)==0:
            return '*'
        s = "; ".join([ '%s=%s' % (k,v) for k,v in self.tags.iteritems()])
        return s

    def _setDefaults( self ):
        originalTags = copy.copy( self.tags )
        if TAG_PRIMARY_ANALYSIS_PROTOCOL not in self.tags:
            self.setTag( TAG_PRIMARY_ANALYSIS_PROTOCOL, DEFAULT_PROTOCOL )
        if TAG_FILE_TYPE not in self.tags:
            self.setTag( TAG_FILE_TYPE, 'fasta' )
        return originalTags

    def resolve( self, inputDataUrl ):
        """Resolves either a InputDataUrl or a string which looks 
        like such a URI."""
        originalTags = self._setDefaults()
        result = []
        path = ""
        if isinstance( inputDataUrl, str ):
            # resolve old URI format
            if ':' in inputDataUrl:
                result = self._resolveOldUri( inputDataUrl )
            else:
                path = inputDataUrl
        else:
            path = inputDataUrl.path
        if len(result)==0:
            if path.startswith('Run/'):
                path = path[path.find(':')+1:]
                # convert run code to old format for now
                result = self._resolveRun( path.replace('/','-') )
            else:
                result = self._resolveOldUri( inputDataUrl )
        self.tags = originalTags
        return result

    def _resolveOldUri( self, inputDataUrl ):
        """Resolves either a InputDataUrl or a string which looks 
        like such a URL."""
        if isinstance( inputDataUrl, str ):
            if inputDataUrl.find(':')>=0:
                scheme, path = inputDataUrl.split(':')
            else:
                scheme, path = 'file', inputDataUrl
        else:
            scheme, path = inputDataUrl.scheme, inputDataUrl.path
        result = []
        if scheme=='run':
            result = self._resolveRun( path )
        if scheme=='file':
            result = self._resolveFile( path )
        if scheme=='contigs':
            result = self._resolveFile( path )
        if scheme=='short_read_contigs':
            result = self._resolveFile( path )
        if scheme=='assembled_contigs':
            result = self._resolveFile( path )
        if scheme=='fastq':
            result = self._resolveFile( path )
        if scheme=='fasta':
            result = self._resolveFile( path )
        if scheme=='string':
            result = [ path ]
        return result

    def getLastError( self ):
        return self._error

    def _resolveFile( self, fileName ):
        for searchPath in  self._paths:
            path = os.path.abspath( os.path.join( searchPath, fileName ) )
            if os.path.exists( path ):
                return [ path ]
        self._error = "Can't find data file %s; PATH=%s" % \
            ( fileName, ":".join(self._paths) )
        return []

    def _resolveAnalysisResults( self, runPath ):
        # we have a path to the data so return the requested file(s)
        if TAG_PRIMARY_ANALYSIS_PROTOCOL in self.tags:
            reportsFolders = self.tags[ TAG_PRIMARY_ANALYSIS_PROTOCOL ]
        movieTypes = []
        if TAG_MOVIE_TYPE in self.tags:
            movieTypes = self.tags[ TAG_MOVIE_TYPE ]
        files = set()
        for reportsFolder in reportsFolders:
            if len(movieTypes)>0:
                #
                # special support for p*, wildcard pane selector
                #
                if 'p*' in movieTypes:
                    files |= set( findHdfFiles( \
                        runPath, reportsFolder, movie='_p[0-9]' ) )
                else:
                    for movieType in movieTypes:
                        _files = findHdfFiles( runPath, reportsFolder, movieType )
                        #print >>sys.stderr, 'runPath = %s'%runPath
                        #print >>sys.stderr, 'hdf_files = %s'%str(_files)
                        files |= set( _files )
            else:
                _files = findHdfFiles( runPath, reportsFolder )
                files |= set( _files )
        if len(files)==0:
            self._error = "No data found under reports=%s movieType=%s" % \
                ( reportsFolders, movieTypes )
        else:
            # filter by file type
            allFiles =[]
            if TAG_FILE_TYPE in self.tags:
                for fileType in self.tags[TAG_FILE_TYPE]:
                    allFiles.extend( filter( lambda x: os.path.exists(x), \
                        map( self._findFile(fileType), files ) ) )

            if len(allFiles)==0:
                self._error = "No data of type %s"\
                    " found under reports=%s movieType=%s" % \
                    ( self.tags[TAG_FILE_TYPE], reportsFolders, movieTypes )
            files = allFiles
        return files

    def _findFile(self, fileType):
        def findFile(f):
            path, stem = os.path.split( f )
            if stem.find('.pls.h5') != -1:
                stem = stem[ :stem.find('.pls.h5') ]
            else:
                stem = stem[ :stem.find('.bas.h5') ]
            return os.path.join( path, stem+'.'+fileType )
        return findFile

    def _resolveMovieObject( self, runPath, ext='.avi' ):
        if TAG_MOVIE_TYPE in self.tags:
            files = []
            movieTypes = self.tags[ TAG_MOVIE_TYPE ]
            for movieType in movieTypes:
                files.extend(glob.glob( '%s/*%s*%s' % ( runPath, movieType, ext)))
        else:
            files = glob.glob( '%s/*%s' % (runPath,ext) )
        if len(files)==0:
            self._error = "No %s files found at %s" % ( ext, runPath )
        return files

    def _resolveRun( self, runName ):
        runPath = findPathForChip( runName )
        if runPath is None:
            self._error = "Can't find data path for %s" % runName
            return []
        # we have a path to the data so return the requested file(s)
        files = []
        if TAG_FILE_TYPE in self.tags:
            for fileType in self.tags[TAG_FILE_TYPE]:
                if fileType=='avi' or fileType=='mmd.h5' or fileType=='mcd.h5':
                    files.extend( self._resolveMovieObject( runPath, ext='.'+fileType ) )
                    #print >>sys.stderr, files
        if len(files)==0:
            files.extend( self._resolveAnalysisResults( runPath ) )
        return files
