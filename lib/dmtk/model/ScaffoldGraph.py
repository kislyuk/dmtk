#!/usr/bin/env python
import sys
import os
import re

from dmtk.model.Range import Range, Ranges

class Node:
    def __init__( self, name, length, orientation, desc ):
        self.name = name
        self.length = int(length)
        self.orientation = orientation
        self.description = desc
        self.color = "lightblue"
        self.shape = "house"
        self.style = "filled"

    def __str__( self ):
        return 'contig_%s [ label="%s (%d)", color="%s", style="%s", shape="%s", orientation=%d ]' % \
                    ( self.name, self.name, self.length, self.color, self.style, self.shape, 90 if self.orientation == "+" else -90 )

class Edge:
    def __init__( self, id1, id2, weight, span, desc ):
        self.id1 = id1
        self.id2 = id2
        self.weight = weight
        self.span = span
        self.description = desc
        self.color = "black"

    def __str__( self ):
        return 'contig_%s -> contig_%s [ label="%d (%d)", color="%s" ]' % ( self.id1, self.id2, self.weight, self.span, self.color )

class SubGraph:
    def __init__( self, id ):
        self.edges = []
        self.nodes = {}
        self.id = id
        self.label = ""

    def setLabel( self, label ):
        self.label = label

    def addNode( self, name, length, orientation, desc ):
        self.nodes[ name ] = Node( name, length, orientation, desc )

    def addEdge( self, n1, n2, weight, span, desc ):
        self.edges.append( Edge(n1,n2,weight,span,desc) )

    def __str__( self ):
        return 'subgraph cluster_scaff_%s {\n' % (self.id) +   \
                        self.label + "\n" + \
                        '\n'.join([str(n) for n in self.nodes.values()]+[str(e) for e in self.edges]) + \
                    '}'

    def startNodes( self ):
        """Returns a set of IDs for 'start nodes', defined as nodes without
        leading input edges."""
        startNodes = set(self.nodes.keys())
        for edge in self.edges:
            startNodes.discard( edge.id2 )
        return startNodes

    def breadthFirst( self, start=None ):
        """Performs a breadth-first search of the graph, starting from
        all startNodes simultaneously. Returns a set of edges at each step."""
        currentNodes = set([start]) if start != None else self.startNodes()
        while len(currentNodes) != 0:
            nextEdges = set()
            nextNodes = set()
            for edge in self.edges:
                if edge.id1 in currentNodes:
                    nextEdges.add( edge )
                    nextNodes.add( edge.id2 )
            yield nextEdges
            currentNodes = nextNodes

    def updateSpans( self ):
        CONTIG_POSNS = re.compile( r'(\d+)\s+(\d+)' )
        SPAN_RE = re.compile( r'\d+\s+\(([-\d]+)\)' )
        for e in self.edges:
            
            # check for our format
            search = SPAN_RE.search( e.description )
            if search:
                e.span = int(search.group(1))
                continue

            # If not, check for Bambus format
            n1 = self.nodes[e.id1]
            n2 = self.nodes[e.id2]
            search1 = CONTIG_POSNS.search( n1.description )
            search2 = CONTIG_POSNS.search( n2.description )
            if not ( search1 and search2 ):
                raise IOError, "Unable to find spans for contigs (%s,%s)" % ( n1.description, n2.description )
            e.span = int(search2.group(1)) - int(search1.group(2))

class Graph:

    def __init__( self ):
        pass

    def load( self, dotFile ):
        infile = open( dotFile, 'r' )
        self.graphs = []
        WEIGHT = re.compile( r'label\s*=\s*"(\d+)' )
        ZERO_WEIGHT = re.compile( r'label\s*=\s*""' )

        NODE_DESC = re.compile( r'label\s*=\s*"([^"]+)"' )
        EDGE_DESC = re.compile( r'label\s*=\s*"([^"]+)"' )
        
        ALT_CONTIG = re.compile( r'([\S]+)' ) #r'contig_([\S]+)' )
        ALT_ORIENT = re.compile( r'orientation\s*=\s*(-?\d+) ' )
        ALT_LENGTH = re.compile( r'label\s*=\s*"[^\(]+\((\d+)\)' )

        SUBGRAPH = re.compile( r'subgraph cluster_scaff_(\d+)' )

        currentGraph = None
        for line in infile:
            
            search = SUBGRAPH.search(line)
            if search:
                currentGraph = SubGraph( search.group(1) )
                self.graphs.append( currentGraph )
                continue    
            
            if line.strip().startswith("label"):
                currentGraph.setLabel( line.strip() )
                continue

            search = ALT_ORIENT.search(line)
            if search:
                o = int( search.group(1) )
                if o<0:
                    orient = '+'
                else:
                    orient = '-'
                search = ALT_CONTIG.search(line)
                name = search.group(1)
                search = ALT_LENGTH.search(line)
                length = search.group(1)
                search = NODE_DESC.search(line)
                desc = search.group(1)

                if currentGraph == None:
                    raise IOError, "Error parsing dot file: No subgraph line found."
                currentGraph.addNode( name, length, orient, desc )    
                continue
            
            if '->' in line:
                values = line.strip().split()
                search = ZERO_WEIGHT.search(line)
                if search:
                    continue
                n1 = values[0]
                n2 = values[2]
                search = ALT_CONTIG.search(n1)
                n1 = search.group(1)
                search = ALT_CONTIG.search(n2)
                n2 = search.group(1)
                
                search = WEIGHT.search( line )
                if search:
                    weight = int(search.group(1))
                else:
                    weight = 0

                search = EDGE_DESC.search(line)
                desc = search.group(1)

                currentGraph.addEdge( n1, n2, weight, 0, desc )

        infile.close()
        for g in self.graphs:
            g.updateSpans()

    def __str__( self ):
        myRep = """digraph ROOT {
                rankdir = LR
                orientation = portrait
                ranksep = 0.3
                nodesep = 0.3
                fontsize = 12
                margin = ".2,.2"
                """
            
        myRep += "\n".join( [ str(g) for g in self.graphs ] )
        myRep += "\n}"
        return myRep

    def writeDot( self, outFN ):
        out = open( outFN, 'w' )
        out.write(str(self))
        out.close()

    def __iter__( self ):
        return iter(self.graphs)

    def __len__( self ):
        return len(self.graphs)


if __name__=="__main__":

    mapping = {}
    for line in open( sys.argv[1], 'r' ):
        line = line.replace("\"","")
        elts = line.strip().split(",")
        contig = elts[0]
        mapping[ contig ] = ( elts[2], elts[7], elts[8] )

    g = Graph( )
    g.load( sys.argv[2] )
    print ",".join([ "c1", "c2", "span", "weight", "minSpan"])
    for sg in g:
        for edge in sg.edges:
            n1 = sg.nodes[ edge.id1 ]
            n2 = sg.nodes[ edge.id2 ]
            
            if n1.orientation == "-" and n2.orientation == "+":
                c1 = edge.id1.split("_")[0]
                c2 = edge.id2.split("_")[0]
                m1 = mapping[c1]
                m2 = mapping[c2]

                score = "NA"
                if m1[0] == m2[0]:
                    score = min(  abs(max(int(m1[1]),int(m1[2])) - min(int(m2[1]),int(m2[2]))),
                             abs(max(int(m2[1]),int(m2[2])) - min(int(m1[1]),int(m1[2]))) )


                print ",".join([ c1, c2, str(edge.span), str(edge.weight), str(score) ])
