__doc__="""Base class for representing a single pairwise overlap. Inspired by Amos Overlap object. It has slightly more information than an AmosOverlapMessage, but less than an AlignmentHit."""
import os
import sys
import logging
import networkx as nx
from copy import copy


class Overlap:

    def __init__( self, ahg=0, bhg=0, uid1=0, uid2=0, adj="N", eid1="", eid2="", len1=0, len2=0, iid=0 ):
        """Does not inherit from AMOS Overlap_t to avoid direct swig dependency."""
        self.ahg  = ahg
        self.bhg  = bhg
        self.iid  = iid
        self.uid1 = uid1
        self.uid2 = uid2
        self.adj  = adj
        # extra fields not present in Amos Overlap object, but necessary to do anything useful
        self.eid1 = eid1
        self.eid2 = eid2
        self.len1 = len1
        self.len2 = len2

    def fieldsFromAmos(self, o, r1, r2):
        """Populates the object from an Amos Overlap_t and two Read_t Read_t object"""
        assert( o.getFirstRead() == r1.getIID() and o.getSecondRead() == r2.getIID() )
        self.ahg  = o.getAhang()
        self.bhg  = o.getBhang()
        self.uid1 = o.getFirstRead()
        self.uid2 = o.getSecondRead()
        self.adj  = o.getAdjacency()
        self.iid = o.getIID()
        self.eid1 = r1.getEID()
        self.eid2 = r2.getEID()
        self.len1 = r1.getLength()
        self.len2 = r2.getLength()
        self.scr = o.getScore()

    def __str__( self ):
        buffer =  "uid1:%i " % self.uid1
        buffer += "uid2:%i " % self.uid2
        buffer += "iid:%i " % self.iid
        buffer += "ahg:%i " % self.ahg
        buffer += "bhg:%i " % self.bhg
        buffer += "adj:%s " % self.adj
        buffer += "eid1:%s " % self.eid1
        buffer += "eid2:%s " % self.eid2
        buffer += "len1:%s " % self.len1
        buffer += "len2:%s " % self.len2
        return buffer

    def inferOverlap(self, o, minOverlapLength=50):
        """Infers an overlap between self.uid1 and o.uid2 via self.uid2 and o.uid1
        If no actual overlap of at least minOverlapLength returns None. 
        Set minOverlapLength to None to infer a new overlap regardless of whether they actually 
        physically overlap. Useful for coordinate transformations."""
        
        if minOverlapLength == None: minOverlapLength = -sys.maxint 
        assert( self.uid2 == o.uid1 )

        new = Overlap() 

        if self.adj == "N": new.ahg = self.ahg + o.ahg
        else:               new.ahg = self.ahg - o.bhg

        if new.ahg  >= self.len1 - minOverlapLength: return None # off the end of self!
        if -new.ahg >= o.len2    - minOverlapLength: return None # off the end of o!

        new.bhg = o.len2 + new.ahg - self.len1

        new.adj = "N" if self.adj == o.adj else "I"

        new.uid1 = self.uid1
        new.eid1 = self.eid1
        new.len1 = self.len1

        new.eid2 = o.eid2
        new.uid2 = o.uid2
        new.len2 = o.len2

        return new
    
    def reverse(self):
        self.uid1, self.uid2 = self.uid2, self.uid1
        self.len1, self.len2 = self.len2, self.len1
        self.eid1, self.eid2 = self.eid2, self.eid1
        if self.adj == "N": self.ahg, self.bhg = -self.ahg, -self.bhg
        else:               self.ahg, self.bhg =  self.bhg,  self.ahg       


class Overlaps(object):
    """Organizes a list of Overlaps to allow extraction of overlaps by uid"""

    def __init__(self, overlaps):
        """Input is a list of overlaps"""
        self._overlaps = overlaps
        self.weightAdjustment = 0
        self.overlapGraph = self._overlapGraph()
        self.idsToOverlaps = self._calcIdsToOverlaps()

    def __str__(self):
        return "\n".join( map(lambda o: str(o), self._overlaps) )
    
    def _overlapGraph(self):

        g = nx.DiGraph()
        minScore = 0

        # currently using the negative score as the weight. This may not be the right thing to do, 
        # but since we are checking against the spans I think it's okay
        for o in self._overlaps:

            minScore = min(-o.scr, minScore)


            g.add_edge( "%s+"  % o.eid1,                                \
                        "%s%s" % (o.eid2, "+" if o.adj == "N" else "-"),\
                        weight = -o.scr,
                        length =  o.bhg )

            g.add_edge( "%s-" % o.eid1,                                 \
                        "%s%s" % (o.eid2, "-" if o.adj == "N" else "+"),\
                        weight = -o.scr,
                        length = -o.ahg )

        # create the reverse edges 
        allEdges = list( g.edges(data=True) )
        for node1, node2, data in allEdges:
            g.add_edge( node2, node1, attr_dict={ "weight":data["weight"], "length":-data["length"] } )

        # remove edges that don't extend the right-most end
        # allEdges = list( g.edges(data=True) )
        # for node1, node2, data in allEdges:
        #    if g[node1][node2]["length"] <= 0:
        #        g.remove_edge(node1, node2)

        # make sure all of our weights are positive
        minScore -= 1 
        allEdges = list( g.edges(data=True) )
        for node1, node2, data in allEdges:
            g[node1][node2]["weight"] -= minScore

        return g

    def _getLength(self, name):
        for o in self._overlaps:
            if o.eid1 == name: return o.len1
            if o.eid2 == name: return o.len2

    def findConnectingPairsShortestPath(self, contig1=None, contig2=None,       \
            contig1strand="+", contig2strand="+", minSpan=100, maxSpan=200, \
            dontUse=set(), maxDepth=3, **kwargs):

        contig1length = self._getLength(contig1)
        contig2length = self._getLength(contig2)

        previousPath = []

        while True:
            pathLength = -sys.maxint
            path = []
            try:
                path = nx.algorithms.shortest_paths.generic.shortest_path( \
                            self.overlapGraph,                      \
                            "%s%s" % (contig1, contig1strand),      \
                            "%s%s" % (contig2, contig2strand),      \
                            weighted=True)
                if not path: raise nx.exception.NetworkXError
                pathLength = sum([ 
                    self.overlapGraph[path[i]][path[i+1]]["length"] \
                       for i in range(len(path)-1)])
                pathLength -= contig2length
                logging.debug("Path length %i" % pathLength)
            except nx.exception.NetworkXError, e:
                logging.debug("Could not find a path between %s and %s - %s", contig1, contig2, e)
                path = []
                break
            except KeyError, e:
                logging.debug("Could not find a key %s in graph", e)
                path = []
                break


            if pathLength > minSpan:
                if pathLength > maxSpan: 
                    path = []
                break

            logging.debug("path length: %i %i outside: (%i,%i)" \
                    % (len(path), pathLength, minSpan, maxSpan))
            for idx in range(len(path)-1): self.overlapGraph.remove_edge( path[idx], path[idx+1] )

        noStrandPath = map(lambda x: x.rstrip("-+"), path)
        return [tuple(noStrandPath[idx:idx+2]) for idx in range(len(noStrandPath)-1) ]

    def findConnectingPairsFromDict(self, mydict):
        pairs = self.findConnectingPairs(**mydict)
        if pairs: return pairs
        pairs = self.findConnectingPairsShortestPath(**mydict)
        return pairs

    def getOverlapsForEid(self, eid):
        """Returns all overlaps for a given id"""
        return [] if not self.idsToOverlaps.has_key(eid) else self.idsToOverlaps[eid].values()

    def findConnectingPairs(self, contig1=None, contig2=None, contig1strand="+", minSpan=100, maxSpan=200, dontUse=set(), maxDepth=3, **kwargs):
        """Returns a set of reads that connect contig1 to contig2, with a span between minSpan and maxSpan 
           and not using any reads in dontUse."""
        if contig1 == None or contig2 == None: raise SystemExit, "Null contig id"

        observed  = set(dontUse) - set([contig2])
        #--- starting at the first contig, do a DFS from that contig using all overlaps
        connectingOverlaps = []

        def searchOverlaps(currentOverlaps, connectingOverlaps, depth): 

            returnVal = False
            if depth > maxDepth: return returnVal
            for current in currentOverlaps:

                #--- stop the descent when you reach
                #--- a) a read that has already been visited
                if current.eid2 in observed: continue
                if current.eid2 != contig2: observed.add( current.eid2 )

                #--- note that the origin is the end of contig1
                start = current.ahg - current.len1 if contig1strand == "+" else - current.bhg - current.len1

                #--- b) a read with a min coord > contig max span. it's possible that we could 
                # still loop back around and find the second contig, but maybe unlikely?
                if maxSpan < start: continue

                #--- c) a read with a max coord < 0. again, it's possible that we could 
                # still loop back around (this time in the other direction) but ... 
                maxCoord = current.ahg + current.len2 if contig1strand == "+" else - current.ahg 
                if maxCoord < 0: continue

                #--- d) (our goal) you get to the second contig within the spans
                if (current.eid2 == contig2):
                    logging.debug("contig1: %s contig2: %s minSpan:%i maxSpan:%i start:%i" \
                                % (contig1, contig2, minSpan, maxSpan, start))
                    if (minSpan <= start and start <= maxSpan):
                        connectingOverlaps.append( current )
                        returnVal = True 
                else:
                    toInferFrom = self.getOverlapsForEid( current.eid2 ) 
                    newOverlaps = [ current.inferOverlap(o, minOverlapLength=None) for o in toInferFrom ] 
                    # newOverlaps = self._subsetOverlaps( newOverlaps, contig2 )
                    #--- if we found the contig
                    if searchOverlaps( newOverlaps, connectingOverlaps, depth + 1 ):
                        connectingOverlaps.append( current )
                        returnVal = True 

            # connectingOverlaps ends up with overlap in the form (contig1, *) where * is everyone else
            return returnVal
        
        startingOverlaps = self.getOverlapsForEid( contig1 )
        self._reorderOverlaps( startingOverlaps, contig2 )
        # logging.debug("Overlaps for %s are %s" % (contig1, " ".join(map(str, startingOverlaps))) )
        searchOverlaps(startingOverlaps, connectingOverlaps, 0)

        # construct connecting pairs by examining adjacent overlaps, 
        # starting with contig1 and ending with contig1
        connectingPairs = []
        lastRead = None
        connectingOverlaps.reverse()
        for o in connectingOverlaps:
            if not lastRead: lastRead = contig1
            connectingPairs.append( (lastRead, o.eid2) )
            if o.eid2 != contig2:
                lastRead = o.eid2
            else:
                lastRead = None
        return connectingPairs

    def _subsetOverlaps(self, overlaps, contig2):
        """Subsets overlaps so that we only visit the most likely to be useful. If we can get to the contig
        in one go, we do."""
        return filter( lambda x: x.eid2 == contig2, overlaps )


    def _reorderOverlaps(self, overlaps, contig2):
        """Reorders overlaps so that we visit those with distance 1 to contig2 first. 
        At some point we should use Dijkstra's algorithm."""
        # relying on stable sort below

        overlaps.sort( lambda x,y: -1 if self.idsToOverlaps[x.eid2].has_key( contig2 ) else 1 )
        overlaps.sort( lambda x,y: -1 if x.eid2 == contig2 else 1 )

    def _calcIdsToOverlaps( self ):
        """Creates a dict that maps ids to their respective overlaps."""

        pairs = {}
        for overlap in self._overlaps:

            if pairs.has_key(overlap.eid1): pairs[overlap.eid1][overlap.eid2] = overlap 
            else:                           pairs[overlap.eid1] = { overlap.eid2: overlap }

            reversed = copy(overlap)
            reversed.reverse()

            if pairs.has_key(reversed.eid1): pairs[reversed.eid1][reversed.eid2] = reversed
            else:                            pairs[reversed.eid1] = { reversed.eid2: reversed }

        return pairs

 
if __name__ == "__main__":
    o1 = Overlap(ahg=100, bhg=100, uid1=1, uid2=2, adj="N", eid1="1", eid2="2", len1=200, len2=200 )
    o2 = Overlap(ahg=100, bhg=100, uid1=2, uid2=3, adj="N", eid1="2", eid2="3", len1=200, len2=200 )
    o3 = o1.inferOverlap(o2, minOverlapLength=-sys.maxint)
    print str(o3)
    o3 = o1.inferOverlap(o2)
    print str(o3)

    o1 = Overlap(ahg=100, bhg=100, uid1=1, uid2=2, adj="I", eid1="1", eid2="2", len1=200, len2=200 )
    o2 = Overlap(ahg=-100, bhg=-100, uid1=2, uid2=3, adj="I", eid1="2", eid2="3", len1=200, len2=200 )
    o3 = o1.inferOverlap(o2, minOverlapLength=-sys.maxint)
    print str(o3)
    o3 = o1.inferOverlap(o2)
    print str(o3)
