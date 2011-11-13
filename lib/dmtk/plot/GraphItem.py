__doc__="""Storage classes for simplying XML production of graph
elements (for reporting)."""
import sys
import os
import re
import xml.sax.saxutils as sx
import xml.etree.cElementTree as et

class GraphReportItem:
    """Models an report page containing
    1.  title
    2.  attributes
    3.  graphs
    4.  tables
    """
    def __init__(self,title=''):
        self.title = title
        #Deprecated - should use graphGroups instead
        self.graphs = []
        self.graphGroups = []
        self.tables = []
        self.attributes = {}
        self.__attributesTmp = {}
        self.layout = 'onecolumn'
        self._keyList = []
        self.errors = []
        
        
    #Deprecated - should use graphGroups instead
    def addGraph(self, graphItem):
        self.graphs.append( graphItem )
        
    def addGraphGroup(self, graphGroupItem):
        self.graphGroups.append( graphGroupItem )

    def addErrorMessage( self, errorMsg ):
        self.errors.append( errorMsg )

    def addTable(self, tableItem):
        self.tables.append( tableItem )

    def put(self, key, value, numberValue=None, id=None, hidden=False, hiddenInDetailedReport=False, hiddenInternally=False ):
        """This preserves attribute insertion order."""
        self.__addAttribute( GraphAttribute(key,value,numberValue,id,hidden,hiddenInDetailedReport,hiddenInternally) )

    def __addAttribute(self, graphAttr ):
        self.attributes[ graphAttr._name ] = graphAttr
        self._keyList.append( graphAttr._name )
        
    def putUnordered(self, key, value, numberValue=None, id=None, hidden=False, hiddenInDetailedReport=False, hiddenInternally=False ):
        """Temporarily add attributes to a dictionary. reorder(...) must be called if this method
        is used, or attributes will not make it into the xml. This allows a script to add 
        attributes at any time, whenever the data is available. Order can then be applied just
        before writing out the xml."""
        self.__attributesTmp[ key ] = GraphAttribute(key,value,numberValue,id,hidden,hiddenInDetailedReport,hiddenInternally)
        
        
    def reorder(self, keyOrder ):
        """Define an attribute order """
        #clear list
        self._keyList = []
        for k in keyOrder:
            if k in self.__attributesTmp:
                ga = self.__attributesTmp[k]
                self.__addAttribute( ga )
        
    def toXml(self):
        buffer = []
        buffer.append('<report>')
        buffer.append('<layout %s="true"/>' % self.layout )
        buffer.append('<title>%s</title>' % sx.escape(self.title) )
        buffer.append('<attributes>')
        for k in self._keyList:
            v = self.attributes[k]
            buffer.append( v.toXml() )
#            buffer.append( '<attribute name=%s>%s</attribute>' % \
#                ( sx.quoteattr(k), sx.escape(str(v)) ) )
        buffer.append('</attributes>')
        #TODO - removed references to graph this level. Graphs should only be nested under GroupGroups
        for graph in self.graphs:
            buffer.append( graph.toXml() )
        for graphGroup in self.graphGroups:
            buffer.append( graphGroup.toXml() )
        for table in self.tables:
            buffer.append( table.toXml() )
        if len(self.errors)>0:
            buffer.append( "<errors>%s</errors>" % os.linesep.join(self.errors) )
        buffer.append('</report>')
        return os.linesep.join(buffer)
    
    def saveToFile(self, fileName):
        "Save this report to a file named fileName or a file stream."
        if isinstance(fileName,str):
            outfile = open(fileName,'w')
            print >>outfile, '<?xml version="1.0" encoding="UTF-8"?>'
            print >>outfile, self.toXml()
            outfile.close()
        else:
            print >>fileName, '<?xml version="1.0" encoding="UTF-8"?>'
            print >>fileName, self.toXml()

class GraphAttribute:
    """Models a single attribute of a graph report.
    hidden: refers to whether an attribute is written into metadata.rdf with isForReportsDisplay==true. In this case,
    the overview report shows it.
    hiddenInDetailedReport: refers to whether an attribute is rendered in a detailed report, like variants.xml. 
    
    It's possible to have attributes that you want to show in on, but not the other.
     """
    def __init__( self, name, displayValue, numberValue=None, id=None, hidden=False,hiddenInDetailedReport=False,hiddenInternally=False ):
        if hidden and not hiddenInDetailedReport:
            print >>sys.stderr, "Illegal state: %s cannot be hidden in overview but visible in details!" % name
            #todo - raise exception
        self._name = name
        self._displayValue = str(displayValue)
        self._numberValue = numberValue
        self._id = id
        if self._id is None and self._name:
            self._id = self._defaultIdGenerator()
        if self._numberValue is None and self._displayValue:
            self._numberValue = self._defaultValueGenerator()
        self._hidden = hidden
        self._hiddenInDetailedReport=hiddenInDetailedReport
        self._hiddenInternally=hiddenInternally

    def _defaultIdGenerator( self ):
        "it's bad form to rely on this"
        return self._name\
            .strip()\
            .lower()\
            .replace(' ','_')\
            .replace('%','pct')\
            .replace('#_of_','n_')

    def _defaultValueGenerator( self ):
        "it's bad form to rely on this"
        return re.sub( '[^0-9.]', '', self._displayValue )

    @staticmethod
    def fromElement( element ):
        """Unmarshal an ElementTree.Element into GraphAttribute"""
        name = element.attrib['name']
        displayValue = element.text.strip()
        attr = GraphAttribute( name, displayValue )
        if 'value' in element.attrib:
            attr._numberValue = element.attrib['value']
        if 'id' in element.attrib:
            attr._id = element.attrib['id']
        if 'hidden' in element.attrib:
            attr._hidden = element.attrib['hidden']=='true'
        if 'hiddenInDetailedReport' in element.attrib:
            attr._hiddenInDetailedReport = element.attrib['hiddenInDetailedReport']=='true'
        if 'hiddenInternally' in element.attrib:
            attr._hiddenInternally = element.attrib['hiddenInternally']=='true'
        return attr

    def toXml(self):
        attrs = { 'name': self._name }
        if self._numberValue is not None:
            attrs[ 'value' ] = str(self._numberValue)
        if self._id is not None:
            attrs[ 'id' ] = str(self._id)
        if self._hidden:
            attrs[ 'hidden' ] = 'true'
        if self._hiddenInDetailedReport:
            attrs[ 'hiddenInDetailedReport' ] = 'true'
        if self._hiddenInternally:
            attrs[ 'hiddenInternally' ] = 'true'
        root = et.Element( 'attribute', attrs )
        root.text = self._displayValue
        return et.tostring(root)
            
class GraphGroupItem:
    """Models a logical group of graphs consisting
    1.  title
    2.  legend
    3.  graphs
    """
    def __init__(self,title=None,legend=None):
        self.title = title
        self.legend = legend
        self.graphs = []
        self.thumbnails = []
        
    def addGraph(self,graph):
        self.graphs.append( graph )
        
    def addLegend(self,legend):
        self.legend = legend
        
    def addThumbnailImage(self, path):
        self.thumbnails.append(path)
        
    def toXml(self):
        buffer = []
        buffer.append( '<graphGroup>' )
        if self.title:
            buffer.append( '<title>%s</title>' % sx.escape(self.title) )
        if self.legend:
            buffer.append( '<legend>%s</legend>' % sx.escape(self.legend) )
        for thumb in self.thumbnails:
            buffer.append( '<thumbnail>%s</thumbnail>' % sx.escape(thumb) )
        for graph in self.graphs:
            buffer.append( graph.toXml() )
        buffer.append( '</graphGroup>' )
        return os.linesep.join(buffer)
    
class GraphItem:
    """Models a graph consisting of
    1.  title
    2.  caption
    3.  images
    4.  links to underlying data
    """
    def __init__(self,title='',file=None):
        #deprecated - titles belong to GraphGroupItem
        self.title = title
        self.caption = ''
        self.images = []
        self.links = []
        if file:
            self.__load(file)     
        
    def __load(self, fileName):
        tree = et.ElementTree( file=fileName )
        rootNode = tree
        if rootNode.getroot().tag!='graph':
            print >>sys.stderr, "File %s doesn't look like a graph item" % fileName
            return
        titleNode = rootNode.find('title')
        if titleNode: self.title = titleNode.text
        captionNode = rootNode.find('caption')
        if captionNode: self.caption = captionNode.text
        imageNodes = rootNode.findall('image')
        for imageNode in imageNodes:
            self.images.append( imageNode.text )
        linkNodes = rootNode.findall('link')
        for linkNode in linkNodes:
            format = linkNode.attrib['format']
            path = linkNode.text
            self.addLink( format, path )
        
    def addLink(self,formatType,path):
        self.links.append( ( formatType, path ) )
        
    def addImage(self, path):
        self.images.append(path)
        
    def toXml(self):
        buffer = []
        buffer.append( '<graph>' )
        if self.title:
            buffer.append( '<title>%s</title>' % sx.escape(self.title) )
        if self.caption:
            buffer.append( '<caption>%s</caption>' % sx.escape(self.caption) )
        for image in self.images:
            buffer.append( '<image>%s</image>' % image )
        for formatType, path in self.links:
            buffer.append( '<link format="%s">%s</link>' % ( formatType, path ) )
        buffer.append( '</graph>' )
        return os.linesep.join(buffer)

class TableItem:
    def __init__( self, nrows, ncols ):
        self._nrows = nrows
        self._ncols = ncols
        self.xlabels = [ '' for i in xrange(self._ncols) ]
        self.ylabels = [ '' for i in xrange(self._nrows) ]
        self.values = [ [ '' for i in xrange(self._ncols) ] for j in xrange(self._nrows) ]
        self.xname = ''
        self.yname = ''
        self.title = ''

    def setColumnLabels( self, labels ):
        self.xlabels = labels

    def setColumnName( self, label ):
        self.xname = label

    def setRowLabels( self, labels ):
        self.ylabels = labels

    def setRowName( self, label ):
        self.yname = label

    def setTitle( self, title ):
        self.title = title

    def set( self, row, col, value ):
        self.values[ row ][ col ] = value

    def toXml( self ):
        root = et.Element( 'table' )

        if self.title:
            title = et.SubElement( root, 'title' )
            title.text = self.title

        if len(self.xlabels)>0:
            h1 = et.SubElement( root, 'thead' )
            if self.xname:
                row = et.SubElement( h1, 'tr' )
                if len(self.ylabels)>0:
                    col = et.SubElement( row, 'th' )
                col = et.SubElement( row, 'th', {'colspan':str(self._ncols)} )
                col.text = self.xname
            row = et.SubElement( h1, 'tr' )
            if len(self.ylabels)>0:
                col = et.SubElement( row, 'th' )
                col.text = self.yname
            for col_name in self.xlabels:
                col = et.SubElement( row, 'th' )
                col.text = col_name

        body = et.SubElement( root, 'tbody' )
        for i in xrange( self._nrows ):
            row = et.SubElement( body, 'tr' )
            if len(self.ylabels)>0:
                item = et.SubElement( row, 'td', {'class':'rowname'} )
                item.text = self.ylabels[i]
            for j in xrange( 0, self._ncols ):
                item = et.SubElement( row, 'td' )
                item.text = self.values[ i ][ j ]

        return et.tostring(root)

#
# helper function for legacy report makers
#
def addReportAttribute( report, attrSpec, displayValue, hidden=False, hiddenInDetailedReport=False, useInsertionOrder=True, hiddenInternally=False, convertDisplayToNumber=True ):
    """attrSpec is tuple of form (name,id). useInsertionOrder is the default. If you set it to
    false, you must supply the report with an ordering after all values have been inserted.
    hidden = hidden in overview report (triggers a boolean value in metadata.rdf)
    hiddenInDetailedReport = hidden in an individual report, like variants.xml. In this case, the attibute is written out with attr hidden=true
    """
    # get more sophisicated if this doesn't work
    if convertDisplayToNumber:
        value = re.sub( '[^0-9.]', '', displayValue )
    else:
        value = displayValue
    
    if useInsertionOrder:
        report.put( attrSpec[0], displayValue, numberValue=value,
            id=attrSpec[1], hidden=hidden, hiddenInDetailedReport=hiddenInDetailedReport, hiddenInternally=hiddenInternally )
    else:
        report.putUnordered( attrSpec[0], displayValue, numberValue=value,
            id=attrSpec[1], hidden=hidden, hiddenInDetailedReport=hiddenInDetailedReport, hiddenInternally=hiddenInternally )

def test_table():
    table = TableItem( 3, 4 )
    table.setColumnName("columns")
    table.setRowName("rows")
    table.setColumnLabels( [ 'A', 'B', 'C', 'D' ] )
    table.setRowLabels( [ '0', '1', '2' ] )
    table.set( 1, 1, '3' )
    table.set( 1, 2, 'Hello' )
    table.set( 2, 1, 'test' )
    print table.toXml()

def test_graph_report():
    # feel free to add more test code here for helping verify the XML outputs
    # of this module
    report = GraphReportItem( title='Testing' )
    report.put( 'my&attribute', '0', numberValue=0, id='my.attribute' )
    report.put( 'my&attribute 2', '0', numberValue=0 )
    graphGroup = GraphGroupItem( title='Graph Group 1' )
    report.addGraphGroup( graphGroup )
    print report.toXml()

def test_graph_attr_unmarshal():
    attrElement = et.Element( 'attribute', \
        { 'name':'test_attr', 'id':'test_id', \
          'value':'100', 'hidden':'true' } )
    attrElement.text = '100 bp'
    attr = GraphAttribute.fromElement( attrElement )
    print 'Original: '
    print et.tostring(attrElement)
    print 'Unmarshalled: '
    print attr.toXml()

if __name__=='__main__':
    test_table()
    test_graph_report()
    test_graph_attr_unmarshal()
