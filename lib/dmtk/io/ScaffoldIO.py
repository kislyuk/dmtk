__doc__="""Provides IO operations for Scaffold objects and associated files (Links, ScaffoldGraphs, bambus output)"""
import sys
import os
import time
from dmtk.model.ScaffoldLink import BaseLink
import networkx as nx
from xml.etree.ElementTree import ElementTree
import logging
from dmtk.io.FastaIO import SimpleFastaReader

#--- Bambus uses these constants to indicate edge validity
VALIDITY_VALID          = "VALID"   # a valid edge
VALIDITY_LENGTH         = "LEN"     # invalidated due to length
VALIDITY_UNSEEN         = "UNSEEN"  
VALIDITY_ORIENTATION    = "ORI"     # invalidated due to orientation
VALIDITY_UNKNOWN_REASON = "UNKNOWN"

#--- TODO give the edges a type (i.e. SequenceLink, BambusLink etc.)
class LinkDAO:

    DEFAULT_XCOORD = "-1"

    def __init__(self):
        pass

    def _linkToEdgeTuple(self, link):
        edgeDict = {"seq1" : link.getFirstId(), "seq2" : link.getSecondId(), \
                    "min"  : link.getMinSpan(), "max"  : link.getMaxSpan(),  \
                    "seq1_strand" : link.getFirstStrand(),                   \
                    "seq2_strand" : link.getSecondStrand(),                  \
                    "used" : "True", "type" : link.getTypeString(),          \
                    "validity" : VALIDITY_VALID }

        # TODO do we need these fields to be consistent with below?
        #        edgeDict["seq1_start"]  = seq1loc[2]
        #        edgeDict["seq2_start"]  = seq2loc[2]
        #        edgeDict["seq1_stop"]   = seq1loc[3]
        #        edgeDict["seq2_stop"]   = seq2loc[3]

        return (link.getFirstId(), link.getSecondId(), edgeDict)
        
    def _createGraph(self, links):
        graph = nx.MultiDiGraph()
        for link in links:
            graph.add_node( link.getFirstId(),  length=link.getFirstLength(),   \
                            strand=link.getFirstStrand(), xcoord=LinkDAO.DEFAULT_XCOORD, \
                            repeat="true" if link.firstIsRepeat else "false")
            graph.add_node( link.getSecondId(), length=link.getSecondLength(), \
                            strand=link.getSecondStrand(), xcoord=LinkDAO.DEFAULT_XCOORD, \
                            repeat="true" if link.secondIsRepeat else "false")
            edgeTup = self._linkToEdgeTuple(link) 
            graph.add_edges_from([ edgeTup ])
        return graph

    def outputLinks(self, links, outputFile, outputHeader=True):
        graph = self._createGraph(links)
        nx.write_graphml(graph, outputFile)

    def close(self): 
        pass


# TODO write a read links routine?

class BambusGraphParser:
    """Parses bambus output into a networkx MultiDiGraph. A MultiDiGraph is a directional 
       graph that allows multiple edges between nodes. In this case, contigs are nodes and 
       sequence inserts are edges. It designates edges as to whether or not they are used in 
       the 'SCAFFOLD' region of the out.xml file by the boolean used. Edges (both used and unused) 
       can be invalidated for a particular reason which is noted in the 'validity' attribute."""

    def __init__(self, prefix1, prefix2 = None):
        self.prefix1 = prefix1
        if (prefix2 == None):
            prefix2 = prefix1 + ".out"
        self.prefix2 = prefix2
            
        self.graph = nx.MultiDiGraph()
        self.insert2seqs   = {} # maps insertid -> (seq1, seq2, min, max)
        self.seq2contigloc = {} # maps seqid    -> (contig, strand, start, stop)
        self.contig2length = {} # maps contigid -> length
        self.contig2orientation = {} # maps contigid -> orientation
        self.seqname2seqid = {} # maps seqname  -> seqid

    def _parseDetectiveFile(self, file):

        tree = ElementTree()
        tree.parse(file)
        for library in tree.getiterator(tag="LIBRARY"):
            for insert in library.getiterator(tag="INSERT"):
                seqs = list(insert.getiterator(tag="SEQUENCE"))
                self.insert2seqs[insert.attrib["ID"]] =     \
                                (seqs[0].attrib["NAME"],    \
                                seqs[1].attrib["NAME"],     \
                                library.attrib["MIN"],      \
                                library.attrib["MAX"])

                self.seqname2seqid[ seqs[0].attrib["NAME"] ] = seqs[0].attrib["ID"]
                self.seqname2seqid[ seqs[1].attrib["NAME"] ] = seqs[1].attrib["ID"]

        for contig in tree.getiterator(tag="CONTIG"):
            contigId = self._toContigId(contig.attrib["ID"]) 
            self.contig2length[ contigId ] = contig.attrib["LEN"]
            for seq in contig.getiterator(tag="SEQUENCE"):
                self.seq2contigloc[ seq.attrib["ID"] ] =            \
                       (contigId,                                  \
                        "+" if seq.attrib["ORI"] == "BE" else "-",  \
                        seq.attrib["ASM_LEND"], 
                        seq.attrib["ASM_REND"])

    def _linkToEdgeTuple(self, link, used=True, validity=None):
        tup = self.insert2seqs[ link.attrib["ID"] ]
        edgeDict = {"seq1":tup[0], "seq2":tup[1], "min":tup[2], "max":tup[3]}
        seq1loc = self.seq2contigloc[ self.seqname2seqid[ edgeDict["seq1"] ] ]
        seq2loc = self.seq2contigloc[ self.seqname2seqid[ edgeDict["seq2"] ] ]
        # check for reversed orientation and swap edge direction if necessary
        if (self.contig2orientation[seq1loc[0]] != seq1loc[1]):
            edgeDict = {"seq1":tup[1], "seq2":tup[0], "min":tup[2], "max":tup[3]}
            seq1loc = self.seq2contigloc[ self.seqname2seqid[ edgeDict["seq1"] ] ]
            seq2loc = self.seq2contigloc[ self.seqname2seqid[ edgeDict["seq2"] ] ]

        edgeDict["seq1_strand"] = seq1loc[1]
        edgeDict["seq2_strand"] = seq2loc[1]
        edgeDict["seq1_start"]  = seq1loc[2]
        edgeDict["seq2_start"]  = seq2loc[2]
        edgeDict["seq1_stop"]   = seq1loc[3]
        edgeDict["seq2_stop"]   = seq2loc[3]
        edgeDict["used"]        = "True" if used else "False"
        edgeDict["validity"]    = validity 
        edgeDict["weight"]      = 1
        edgeDict["type"]      = "BambusLink"
        edgeTuple = (seq1loc[0], seq2loc[0], edgeDict)
        return edgeTuple

    def _parseValidity(self, validityString):
            validity = VALIDITY_VALID       if validityString == VALIDITY_VALID       else \
                       VALIDITY_LENGTH      if validityString == VALIDITY_LENGTH      else \
                       VALIDITY_UNSEEN      if validityString == VALIDITY_UNSEEN      else \
                       VALIDITY_ORIENTATION if validityString == VALIDITY_ORIENTATION else \
                       VALIDITY_UNKNOWN_REASON
            return validity

    def _toContigId(self, bambusId):
        # take off the corrupting "contig_" prefix
        if bambusId.startswith("contig_"):
            bambusId = bambusId[7:]
        else:
            logging.warn("Strange bambus id does not begin with contig_:%s" % bambusId)
        return bambusId

    def _parseOutFile(self, file):
        tree = ElementTree()
        tree.parse(file)
        for scaffold in tree.getiterator(tag="SCAFFOLD"):
            for contig in scaffold.getiterator("CONTIG"):
                contigId = self._toContigId(contig.attrib["ID"])
                self.contig2orientation[ contigId ] = "+" if contig.attrib["ORI"] == "BE" else "-"
                self.graph.add_node(contigId,
                                    length=int(self.contig2length[contigId]),
                                    xcoord=contig.attrib["X"],
                                    strand="+" if contig.attrib["ORI"] == "BE" else "-",
                                    repeat="false", 
                                    compositionNodes="%s+" % contigId)
            for link in scaffold.getiterator("LINK"):
                validity = self._parseValidity( link.attrib["VALID"] )
                self.graph.add_edges_from([ self._linkToEdgeTuple(link, used=True, validity=validity) ])

        for link in tree.find("UNUSED").getiterator("LINK"):
            validity = self._parseValidity( link.attrib["VALID"] )
            self.graph.add_edges_from([ self._linkToEdgeTuple(link, used=False, validity=validity) ])


    def parseMultiDiGraph(self):
        """Parses the bambus output files and returns a networkx MultiDiGraph."""
        
        detectiveFile = "%s.detective.xml" % self.prefix1
        outFile = "%s.xml" % self.prefix2
        try:
            self._parseDetectiveFile(detectiveFile)
        except Exception, e:
            raise SystemExit, "Trouble parsing %s: %s" % (detectiveFile, e)
        try:
            self._parseOutFile(outFile)
        except Exception, e:
            raise SystemExit, "Trouble parsing %s: %s" % (outFile, e)
        return self.graph

class AbyssGraphParser:
    """Parses abyss intermediate files into a networkx MultiDiGraph. A MultiDiGraph is a directional 
       graph that allows multiple edges between nodes. In this case, contigs (with directionality!) 
       are nodes and entries in the abyss adj (or dist) file for single (paired end) sequencing
       are edges. Note that the directionality means that each contig is represented as two nodes in the graph."""

    def __init__(self, prefix):
        self.prefix = prefix
        self.singleGraph = nx.DiGraph()
        self.pairedGraph = nx.DiGraph()
        self._contigToComposition = {}

    # TODO move somewhere more sensible?
    def rc(self, id): return id.replace("+", "a").replace("-", "+").replace("a","-")

    def _addSingleEdges(self, id, succ, pred):
        rc = self.rc # for brevity

        posEdges = zip( [id]*len(succ), succ ) + zip(pred, [id]*len(pred) )
        self.singleGraph.add_edges_from(posEdges)

        negEdges = zip( [rc(id)]*len(pred), map(rc, pred) ) + zip( map(rc, succ), [rc(id)]*len(succ) )
        self.singleGraph.add_edges_from(negEdges)

    def _parseSingleContigAdjacencies(self, prefix=None):
        """Assumes the suffix used is the abyss standard (TODO is it?) of prefix-3.adj. 
           Populates self.singleGraph."""
        adjacencyFile = "%s-3.adj" % (prefix if prefix else self.prefix)
        for line in open(adjacencyFile): 
            line = line.rstrip("\n")
            succ, pred = line.split(";")
            pred = pred.rstrip(" ").lstrip(" ").split(" ")
            succ = succ.rstrip(" ").split(" ")
            id, length, succ = "%s+" % succ[0], succ[1], succ[2:]
            succ = [] if len(succ) == 1 and succ[0] == "" else succ
            pred = [] if len(pred) == 1 and pred[0] == "" else pred
            self._addSingleEdges(id, succ, pred)

        singleEndContigsFile = "%s-3.fa" % (prefix if prefix else self.prefix)
        for entry in SimpleFastaReader(singleEndContigsFile):
            prefix = entry.name.split(" ")[0] 
            self.singleGraph.node[ "%s+" % prefix ]["length"] = len(entry.sequence)
            self.singleGraph.node[ "%s-" % prefix ]["length"] = len(entry.sequence)

    def _parsePairedContigAdjacencies(self, prefix=None):
        """Assumes the suffix is the abyss standard (TODO is it?) of prefix-contigs.fa. 
           Populates self.pairedGraph"""
        contigsFile = "%s-contigs.fa" % (prefix if prefix else self.prefix)
        for line in open(contigsFile):
            if not line.startswith(">"): continue
            line = line.rstrip("\n").lstrip(">")
            fields = line.split(" ")
            id = fields[0]
            if len(fields) < 4: 
                self._contigToComposition[id] = [ "%s+" % fields[0] ] # default to single contig comp
                continue 

            self._contigToComposition[id] = fields[-1].split(",")

        # map the single end contig that starts each contig to its respective contig id (including strands for both)
        singleStartToContigs = nx.DiGraph()
        for contig, composition in self._contigToComposition.items():
            start = composition[0]
            singleStartToContigs.add_edge( start, "%s+" % contig )
            end = self.rc( composition[-1] )
            singleStartToContigs.add_edge( end, "%s-" % contig )
        
        # iterate through the contig and determine if their terminal single ends go any
        # other contigs single ends, and note the edge
        def successors(graph, node): return graph.successors(node) if graph.has_node(node) else []

        for contig, composition in self._contigToComposition.items():
            end = composition[-1] # see where we can go from this contig's end
            for singleSuccessor in successors(self.singleGraph, end):
                for contigSuccessor in successors(singleStartToContigs, singleSuccessor):
                    self.pairedGraph.add_edge( "%s+" % contig, contigSuccessor )

            rcStart = self.rc( composition[0] ) # see where we can go from this contigs rc'ed start
            for singleSuccessor in successors(self.singleGraph, rcStart):
                for contigSuccessor in successors(singleStartToContigs, singleSuccessor):
                    self.pairedGraph.add_edge( "%s-" % contig, contigSuccessor )

    def parseMultiDiGraph(self, singleGraph=False):
        """Parses the bambus output files and returns a networkx MultiDiGraph."""
        
        self._parseSingleContigAdjacencies(self.prefix)
        if not singleGraph:
            self._parsePairedContigAdjacencies(self.prefix)

        return self.singleGraph if singleGraph else self.pairedGraph

    # TODO sometimes there are adjacencies internal to other contigs. What do we do with these?

if __name__=="__main__":
    ap = AbyssGraphParser(sys.argv[1])
    graph = ap.parseMultiDiGraph(singleGraph=True)
    gmlFile = "abyss.gml"
    nx.write_graphml(graph, gmlFile)

#    bp = BambusGraphParser(sys.argv[1])
#
#    graph = bp.parseMultiDiGraph()
#    gmlFile = "out.gml"
#    nx.write_graphml(graph, gmlFile)
#
#    newgraph = nx.read_graphml(gmlFile)
#    newGmlFile = "newout.gml"
#    nx.write_graphml(newgraph, newGmlFile)
#
#    newergraph = nx.read_graphml(newGmlFile)
#    newGmlFile = "newerout.gml"
#    nx.write_graphml(newergraph, newGmlFile)
#
