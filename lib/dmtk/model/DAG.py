"""Implements a lightweight DAG model."""
import sys

class Node( object ):
    
    def __init__( self, obj ):
        self.inNodes    = set( )
        self.outNodes   = set( )
        self.obj = obj
        
    @property
    def depth( self ):
        if len( self.inNodes ) == 0:
            return 1
        return 1 + max([ node.depth for node in self.inNodes ])
        
    @property
    def URL( self ):
        """Delegated to the underlying object"""
        return self.obj.URL
            
class DAG( object ):
    
    def __init__( self, URL ):
        """Initializes an empty graph""" 
        self.nodes = set()
        self.edges = set()
        self._url2Node = { }
        self.URL = URL
        
    def newNode( self, obj, type ):
        """Creates a new node, adds it to the graph, and returns it."""
        node = self._url2Node.setdefault( obj.URL, type( obj ))
        self.nodes.add( node )
        return node
    
    def addLink( self, parent, child ):
        """Links this node to the specified child node and vice versa."""
        parent.outNodes.add( child )
        child.inNodes.add( parent )
        self.edges.add(( parent, child ))

    def __getitem__( self, url ):
        """Graph["URL"] ==> Node"""
        return self.url2Node[ url ]
    
    def longestPath( self, edge2distance=lambda x: 1 ):
        """Returns a list of nodes representing the longest path within
        this DAG given a function that computes distance from an edge."""
        longestPaths = [ self._longestPath( n, edge2distance ) for n in self.leafNodes ]
        longestPaths.sort( key=lambda x: x[0], reverse = True )
        return longestPaths[0][1] if len(longestPaths) > 0 else [ ]

    
    def _longestPath( self, node, edge2distance=lambda x: 1 ):
        """Returns the (length, longest path) which ends at the given node."""
        if len( node.inNodes ) == 0:
            return ( 0, [ node ] )
        allPaths = [ self._longestPath( n, edge2distance ) for n in node.inNodes ]
        allPaths.sort( key = lambda x: x[0], reverse = True )
        distance = edge2distance( allPaths[0][1][-1], node )
        return ( allPaths[0][0] + distance, allPaths[0][1] + [ node ] )
        
    @property
    def rootNodes( self ):
        """Returns the set of nodes which have no edges leading into them."""
        return set([ node for node in self.nodes if len(node.inNodes) == 0 ])
    
    @property
    def leafNodes( self ):
        """Returns the set of nodes which have no edges leading out of them."""
        return set([ node for node in self.nodes if len(node.outNodes) == 0 ])
    
    @property
    def orphanedNodes( self ):
        """Searches for and removes nodes which have no other nodes associated with them."""
        return self.rootNodes & self.leafNodes
        