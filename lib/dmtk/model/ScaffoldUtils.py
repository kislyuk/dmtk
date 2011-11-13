import sys
import numpy
import networkx as nx
import logging
from copy import deepcopy
import Range

__doc__="""Class for manipulating scaffold graphs"""

VALIDITY_VALID          = "VALID"   # a valid edge
VALIDITY_LENGTH         = "LEN"     # invalidated due to length
VALIDITY_UNSEEN         = "UNSEEN"  
VALIDITY_ORIENTATION    = "ORI"     # invalidated due to orientation
VALIDITY_UNKNOWN_REASON = "UNKNOWN"



class ScaffoldUtils:
    
    def __init__(self, graph, edgeBoundMethod="Median"):
        self.G = graph
        if (edgeBoundMethod == "strict"):            
            self.edgeBound = StrictEdgeBound(self.G)        
        elif (edgeBoundMethod == "Myer"):
            self.edgeBound = MyerEdgeBound()
        else:
            self.edgeBound = MedianEdgeBound(self.G)
        self.nodeTriangleDict = {}

    def convertPathsToMultiGraph (self, paths):
        """Change a path list to a graph with data attributes"""
        pathG = nx.MultiDiGraph()
        pathG.add_nodes_from(self.G.nodes(data=True))
        for path in paths:
            for i in range (len(path)-1):
                if (path[i+1] in self.G[path[i]]):
                    edgeDict = self.G[path[i]][path[i+1]]
                    for edgeAttr in edgeDict.values():
                        pathG.add_edge(path[i], path[i+1], attr_dict=edgeAttr)
                else:
                    #This is not ideal logic - need to fix! 
                    #Tries to find a shared link along path to infer the edge parameters
                    linked = False
                    backI = i-1
                    forwardI = i+2
                    edge = None
                    while (backI >= 0 and not(linked)):
                        if (path[i+1] in self.G[path[backI]] and 
                            path[i] in self.G[path[backI]]):
                            edge = self._fillEdgeMissingWithPredecessor (path[backI], path[i], path[i+1])
                            linked = True
                        backI += -1
                    while (forwardI < len(path) and not(linked)):
                        #if (path[i+1] in self.G[path[forwardI]] and 
                        #    path[i+1] in self.G[path[forwardI]]):
                        if (path[forwardI] in self.G[path[i+1]] and
                            path[forwardI] in self.G[path[i]]):
                            edge = self._fillEdgeMissingWithSuccessor (path[i], path[i+1], path[forwardI])
                            linked = True
                        forwardI += 1
                    if edge:
                        pathG.add_edge(edge[0], edge[1], attr_dict=edge[2])
                    
        return pathG

    def modifyComplexTriangleNodesInGraph (self):
        """Convert nodes in triangles new Nodes.
        This is proof of concept code.  It only handles 
        the case for which the predecessor to the triangle
        only has two successors and which a single node is 
        jumped over."""
        finished = False
        nodeTriangleDict = self.nodeTriangleDict
        while (not(finished)):
            finished = True
            # Iterate through nodes looking for triangles
            for node in self.G.nodes():
                ss = self.G.successors(node)
                # Iterate through successors check and see if they contain links to the 
                # same set of repeat nodes
                if (len(ss) < 3):
                    continue
                for i in range (len(ss)):
                    sharedNodes = ss[0:i]
                    sharedNodes.extend(ss[i+1:])
                    sNode = ss[i]
                    pred_s = self.G.predecessors(sNode)
                    # Check candidate node
                    if (len(sharedNodes) == len(pred_s)-1):
                        allShared = True
                        for sharedNode in sharedNodes:
                            if (not (sharedNode in pred_s)):
                                allShared = False
                        # Check if the set of shared nodes is the same
                        if (allShared):

                            self.sortNodesByLeftCoord (sharedNodes)
                            # Note, for now assume xcoord position sorting works
                            sharedNodesB = []

                            for sharedNodeIndex in range (len(sharedNodes)):
                                sharedNode = sharedNodes[sharedNodeIndex]
                                count = 0
                                if (sharedNode in nodeTriangleDict):
                                    count = nodeTriangleDict[sharedNode]
                                    nodeTriangleDict[sharedNode] += 1
                                else:
                                    nodeTriangleDict[sharedNode] = 1
                                sharedNodeNew = "dup%i/%s" % (count, sharedNode)
                                sharedNodesB.append(sharedNodeNew)
                                newAttrs = self.G.node[sharedNode]
                                self.G.add_node(sharedNodeNew, attr_dict=newAttrs)                        
                                
                            newNodes = []
                            prevNode = node
                            prevNodeB = node
                            for sharedNodeIndex in range (len(sharedNodes)):
                                sharedNode = sharedNodes[sharedNodeIndex]
                                sharedNodeNew = sharedNodesB[sharedNodeIndex]
                                if (sharedNodeIndex+1 < len(sharedNodes)):
                                    nextNodeB = sharedNodesB[sharedNodeIndex+1]
                                    nextNode = sharedNodes[sharedNodeIndex+1]
                                else:
                                    nextNode = sNode
                                newAttrs = self.G.node[sharedNode]
                                # Add edges to new node copy, delete edges from prev node
                                if (prevNode == node):
                                    edgeDictprevToshared = self.G[node][sharedNode]
                                    for edgeAttr in edgeDictprevToshared.values():
                                        self.G.add_edge(node, sharedNodeNew, attr_dict=edgeAttr)
                                else:
                                    nb, nc, prevEdgeAttr = self._fillEdgeMissingWithPredecessor (node, prevNode, sharedNode)
                                    self.G.add_edge(prevNodeB, sharedNodeNew, attr_dict=prevEdgeAttr)
                                
                                # Add edges to new node copy, delete edges from next node
                                if (nextNode == sNode):
                                    edgeDictsharedTosNode = self.G[sharedNode][sNode]
                                    for edgeAttr in edgeDictsharedTosNode.values():
                                        self.G.add_edge(sharedNodeNew, sNode, attr_dict=edgeAttr)
                                else:
                                    na, nb, nextEdgeAttr = self._fillEdgeMissingWithSuccessor (sharedNode, nextNode, sNode)
                                    self.G.add_edge(sharedNodeNew, nextNodeB, attr_dict=nextEdgeAttr)
                                prevNode = sharedNode
                                prevNodeB = sharedNodeNew
                                newNodes.append(sharedNodeNew)
                            
                            # Remove internal edges with the current sharedNode
                            for tempNode in sharedNodes:
                                for temp2Node in sharedNodes:
                                    while (self.G.has_edge(tempNode, temp2Node)):
                                        self.G.remove_edge (tempNode, temp2Node)
                                    while (self.G.has_edge(temp2Node, tempNode)):
                                        self.G.remove_edge (temp2Node, tempNode)
                                
# keeping edges which give evidence for linkage outside the resolved repeat
#                                while (self.G.has_edge(tempNode, node)):
#                                    self.G.remove_edge (tempNode, node)
#                                while (self.G.has_edge(sNode, tempNode)):
#                                    self.G.remove_edge (sNode, tempNode)
 
                                while (self.G.has_edge(tempNode, sNode)):
                                    self.G.remove_edge (tempNode, sNode)
                                    
                                while (self.G.has_edge(node, tempNode)):
                                    self.G.remove_edge (node, tempNode)
                               
                                if (len(self.G.successors(tempNode)) == 0 and
                                    len(self.G.predecessors(tempNode)) == 0):
                                    self.G.remove_node(tempNode)
                            # remove outter edges
                            while (self.G.has_edge(node, sNode)):
                                self.G.remove_edge (node, sNode)
                            finished = False
                            break
                if (not(finished)):
                    break

    def modifySimpleTriangleNodesInGraph (self):
        """Convert nodes in triangles new Nodes.
        This is proof of concept code.  It only handles 
        the case for which the predecessor to the triangle
        only has two successors and which a single node is 
        jumped over."""
        finished = False
        nodeTriangleDict = self.nodeTriangleDict
        while (not(finished)):
            finished = True
            #Iterate through nodes looking for triangles
            for node in self.G.nodes():
                ss = self.G.successors(node)
                nodeCoord = self._getNodeLeftCoord(node)
                if (len(ss) != 2):
                    continue

                ss0Coord = self._getNodeLeftCoord(ss[0])
                ss1Coord = self._getNodeLeftCoord(ss[1])
                #Determine which node is the beginnig/end of the triangle
                if (ss0Coord < nodeCoord or 
                    ss1Coord < nodeCoord):
                    continue
                if (ss0Coord < ss1Coord):
                    ss0 = ss[0]
                    ss1 = ss[1]
                else:
                    ss0 = ss[1]
                    ss1 = ss[0]

                # Ensure that the "triangle" node fits within the span
                if (not(self.IsValidLengthTriangle(node, ss0, ss1)) or 
                    len(self.G.predecessors(ss1)) != 2):
                    edges_n_ss1 = self.G[node][ss1].values()
                    continue
                    
                s0ss = self.G.successors(ss0)
                s1ss = self.G.successors(ss1)

                if (ss1 in s0ss):
                    newAttrs = self.G.node[ss0]
                    if (len(self.G.predecessors(ss0)) > 1 or len(s0ss) > 1): 
                        count = 0
                        if (ss0 in nodeTriangleDict):
                            nodeTriangleDict[ss0] += 1
                            count = nodeTriangleDict[ss0]
                        else:
                            nodeTriangleDict[ss0] = 0
                        ss0new = "dup%i/%s" % (count, ss0)
                        # Add edges to new node copy, delete edges from original node
                        self.G.add_node(ss0new, attr_dict=newAttrs)
                        edgeDictnTos0 = self.G[node][ss0]
                        for edgeAttr in edgeDictnTos0.values():
                            self.G.add_edge(node, ss0new, attr_dict=edgeAttr)
                        while (self.G.has_edge(node, ss0)):
                            self.G.remove_edge (node, ss0)
                        edgeDicts0Tos1 = self.G[ss0][ss1]
                        for edgeAttr in edgeDicts0Tos1.values():
                            self.G.add_edge(ss0new, ss1, attr_dict=edgeAttr)
                        while (self.G.has_edge(ss0, ss1)):
                            self.G.remove_edge (ss0, ss1)
                    # remove outter edges
                    while (self.G.has_edge(node, ss1)):
                        self.G.remove_edge (node, ss1)
                    finished = False
                    break


    def simplifyPathsAndTriangles (self):
        """Linearize Graph and remove triangles"""
        paths = self.getLinearPaths()
        trianglePaths = self.simplifyTriangles (paths)
        return self.convertPathsToMultiGraph(trianglePaths)
    
    def simplifyTriangles (self, paths):
        """Join paths with simple triangles"""
        triangles = self._getSimpleTriangles ()
        return self._joinPathsWithTriangles (paths, triangles)

    def simplifyOpenTriangles (self, paths):
        """Join paths with open triangles"""
        triangles = self._getOpenTriangles ()
        return self._joinPathsWithOpenTriangles (paths, triangles)
    
    def _getOpenTriangles (self):
        """Get depth 1 open triangles that agree with span/scaffold position"""    
        triangles = []
        visited = []
        for n, nbrs in self.G.adjacency_iter():
            if (n in visited):
                continue   
            triangle = self._getOpenTriangle(n)
            if (len(triangle) > 0):
                triangles.append(triangle)
        simpleTriangles = []
        for i in range (len(triangles)):
            simpleCheck = True
            for j in range (len(triangles)):
                if (i != j):
                    for n in triangles[i]:
                        if (n in triangles[j]):
                            simpleCheck = False
            if (simpleCheck):
                simpleTriangles.append(triangles[i])    
        return simpleTriangles

    def _getOpenTriangle (self, n):
        """Check if node is the start of an open triangle
           If triangle, return ordered triangle (as list)
           Otherwise, return empty list"""

        triangleSet = []
        n_ss = self.G.successors(n)
        for n_s1 in n_ss:
            for n_s2 in n_ss:
                if (len(self.G.successors(n_s1)) == 0):
                    if (self.IsValidLengthTriangle(n,n_s1,n_s2)):
                        triangleSet.append([n, n_s1, n_s2])
        if (len(triangleSet) == 1):
            return triangleSet[0]
        for n_s2 in n_ss:
            n_s2_ps = self.G.predecessors(n_s2)
            for n_p1 in n_s2_ps:
                if (len(self.G.successors(n_s2_ps)) == 0):
                    if (self.IsValidLengthTriangle(n,n_s1,n_s2)):
                        triangleSet.append([n, n_s1, n_s2])
        if (len(triangleSet) == 1):
            return triangleSet[0]
        else:
            return []

    def _getSimpleTriangles (self):
        """Get depth 1 triangles that agree with span/scaffold position"""    
        triangles = []
        visited = []
        for n, nbrs in self.G.adjacency_iter():
            if (n in visited):
                continue   
            triangle = self._getSimpleTriangle(n)
            if (len(triangle) > 0):
                triangles.append(triangle)
        simpleTriangles = []
        for i in range (len(triangles)):
            simpleCheck = True
            for j in range (len(triangles)):
                if (i != j):
                    for n in triangles[i]:
                        if (n in triangles[j]):
                            simpleCheck = False
            if (simpleCheck):
                simpleTriangles.append(triangles[i])    
        return simpleTriangles

    def _getSimpleTriangle (self, n):
        """Check if node is the start of a simple triangle
           If triangle, return ordered triangle (as list)
           Otherwise, return empty list"""

        triangleSet = []
        n_ss = self.G.successors(n)
        for n_s1 in n_ss:
            for n_s2 in n_ss:
                if (n_s2 in self.G.successors(n_s1)):
                    if (self.IsValidLengthTriangle(n,n_s1,n_s2)):
                        triangleSet.append([n, n_s1, n_s2])
        if (len(triangleSet) == 1):
            return triangleSet[0]
        else:
            return []

    def IsValidLengthTriangleReverse(self, n, n_s1, n_s2):
        if ('xcoord' in self.G.node[n]):
            coord_n = self._getNodeLeftCoord(n)
            coord_n_end = coord_n + self._getNodeLength(n)

            edges1 = self.G[n_s2][n_s1].values()
            minSpan1 = self.edgeBound.getMinSpan(edges1, n, n_s1)
            
            edges2 = self.G[n][n_s2].values()
            maxSpan2 = self.edgeBound.getMaxSpan(edges2, n, n_s2)

            
            if ('xcoord' in self.G.node[n_s2]):
                coord_n2 = self._getNodeLeftCoord(n_s2)
            else:
                coord_n2 = coord_n_end + maxSpan2
            
            if ('xcoord' in self.G.node[n_s1]):
                coord_n1 = self._getNodeLeftCoord(n_s1)
            else:
                coord_n1 = coord_n2 -  minSpan1
            coord_n1_length = self._getNodeLength(n_s1)

            if (coord_n < coord_n1 and coord_n < coord_n2 and
                maxSpan2 >= (coord_n1_length+minSpan1)):
                return True
        return False


    def IsValidLengthTriangle(self, n, n_s1, n_s2):
        if ('xcoord' in self.G.node[n]):
            coord_n = self._getNodeLeftCoord(n)
            coord_n_end = coord_n + self._getNodeLength(n)

            edges1 = self.G[n][n_s1].values()
            minSpan1 = int(self.edgeBound.getMinSpan(edges1, n, n_s1))
            
            edges2 = self.G[n][n_s2].values()
            maxSpan2 = int(self.edgeBound.getMaxSpan(edges2, n, n_s2))

            if ('xcoord' in self.G.node[n_s1]):
                coord_n1 = self._getNodeLeftCoord(n_s1)
            else:
                coord_n1 = coord_n_end + minSpan1
            coord_n1_length = self._getNodeLength(n_s1)
            if ('xcoord' in self.G.node[n_s2]):
                coord_n2 = self._getNodeLeftCoord(n_s2)
            else:
                coord_n2 = coord_n_end + maxSpan2

            if (coord_n < coord_n1 and coord_n < coord_n2 and
                maxSpan2 >= (coord_n1_length+minSpan1)):
                return True
        return False

            
    def IsStrictTriangle(self, triangle, path, node):
        if (node == triangle[0] and not(path[-1] == node) or 
            (node == triangle[-1] and (not(path[0] == node))) or
            (node in triangle[1:-1] and len(path) != 1)):
            return  False
        else:
            return True

    def _joinPathsWithTriangles (self, paths, triangles):
        validTriangles = []
        selectedPaths = []
        #Loop has two roles
        # 1.) Get all non-triangle paths 
        # 2.) select the triangles that are consistent with linear paths
        for triangle in triangles:
            triangleCheck = False
            strictCheck = True
            pathsInTriangle = []
            for i in range(len(paths)):
                path = paths[i]
                for n in path:
                    if (n in triangle):
                        triangleCheck = True
                        pathsInTriangle.append(i)
                        if(not(self.IsStrictTriangle(triangle, path, n))):
                            strictCheck = False
            if (triangleCheck and strictCheck):
                validTriangles.append(triangle)
                selectedPaths.extend(pathsInTriangle)
        newPaths = []        
        for i in range (len(paths)):
            if (not(i in selectedPaths)):
                newPaths.append(paths[i])
        #Join all paths/singletons incident on valid triangles
        for triangle in validTriangles:
            newPath = []
            for n in triangle:
                for path in paths:
                    if (n in path):
                        newPath.extend(path)
            newPaths.append(newPath)
        return newPaths
    
    def getLinearPathsWithStrand(self):
        """Like getLinearPaths but includes strand of contig at end. e.g. contigName+"""

        withStrand = []
        for path in self.getLinearPaths():
            withStrand.append( map( lambda c: "%s%s" % (c, self.G.node[ c ]["strand"]), path ) )
        return withStrand

    def getLinearPaths(self):
        """Get all non-branching paths in the graph
           Note:  Paths include singletons"""

        visited = []
        paths = []
        for n, nbrs in self.G.adjacency_iter():
            if (n in visited):
                continue   
            #Necessary to avoid circular paths
            walkBackVisited = []
            path_start_n = self._walkBackLinearPath (n, visited, walkBackVisited)
            path = [path_start_n]
            visited.append(path_start_n)
            ss = self.G.successors(path_start_n)
            if (len(ss) > 0 and len(ss) == ss.count(ss[0]) and 
                not(ss[0] in visited)):                
                self._extendLinearPath(ss[0], path, visited)
            paths.append(path)
        return paths

    def _walkBackLinearPath (self, n, visited, walkBackVisited):
        """Get the start of a linear path
           Last predessor node that maintains linearity"""

        walkBackVisited.append(n)
        ps = self.G.predecessors(n)
        if (len(ps) > 0 and len(ps) == ps.count(ps[0])):
            prev = ps[0]
            prev_ss = self.G.successors(prev) 
            if (not(prev in visited) and not(prev in walkBackVisited) and len(prev_ss) == prev_ss.count(n)):
                return self._walkBackLinearPath (prev, visited, walkBackVisited)
        return n

    def _extendLinearPath (self, n, path, visited):
        """Recursively extend path as long as the graph is linear"""
        ps = self.G.predecessors(n)
        ss = self.G.successors(n)
        if (len(ps) > 0 and len(ps) == ps.count(ps[0])):
            path.append(n)
            visited.append(n)            
            if (len(ss) > 0 and len(ss) == ss.count(ss[0]) and
                not(ss[0] in visited)):
                s = ss[0]
                self._extendLinearPath(s, path, visited)


    def getNextNodePossibilities (self, n, visited):
        """Check to see if the placement of the downstream node is ambiguous"""
        ss = self.G.successors(n)
        currSet = []
        # check if node is terminal
        if (len(ss) == 0):
            pass
        # check if next node is unambiguous
        elif (len(ss) == 1 and len(self.G.predecessors(ss[0])) == 1):
            currSet = ss
        # determine ambiguous nodes
        else: 
            possibleNodes = ss
            for s in ss:
                # add all predecessors nodes of the current successor to possibleNodes
                for s_p in self.G.predecessors(s):
                    if (not(s_p == n or s_p in possibleNodes)):
                        possibleNodes.append(s_p)
            currSet = possibleNodes
        outputSet = []
        for node in currSet:
            if (node not in visited): 
                outputSet.append(node)
            else:
                outputSet.append(None)
        return outputSet
                
    def _getNodeLeftCoord(self, node):
        return int(self.G.node[node]["xcoord"]) - \
                (0 if self._nodeOnPositiveStrand(node) else self._getNodeLength(node))

    def sortNodesByLeftCoord (self, nodes):
        """Return nodes sorted by left coord"""
        nodes.sort(lambda x,y: cmp(self._getNodeLeftCoord(x), self._getNodeLeftCoord(y)) )

    def _getMaxLeftCoord(self, nodes):
        """Return max value on a set of nodes"""
        return numpy.max(map(lambda x: self._getNodeLeftCoord, nodes))
        
    def _sortPathsByLeftCoord(self, paths):
        """Return paths sorted by their max xcoord"""
        paths.sort(lambda x,y: cmp(self._getMaxLeftCoord(x), self._getMaxLeftCoord(y)))

    def getEdgeReads(self, edge):
        edgeDataDict = self.G.get_edge_data(edge[0], edge[1])
        if edgeDataDict[0].has_key("reads"): 
            return ",".join(map(lambda x: str(x["reads"]),  edgeDataDict.values())).split(",")
        else:
            return []

    def sumNodeLengths (self, nodes):
        """Return nodes sorted by xcoord"""
        return numpy.sum(map(self._getNodeLength, nodes))
        
    def checkRepeatEdge (self, n):
        leftRepeat = False
        rightRepeat = False
        nextNode = False
        ss = self.G.successors(n)
        for s in ss:
            if (self.G.node[s]['repeat'] == true):
                leftRepeat = s
            else:
                nextNode = s
                for p in self.G.predecessors(s):
                    if (self.G.node[s]['repeat'] == true):
                        rightRepeat = p
        return leftRepeat, rightRepeat, nextNode


    def _repeatFourCaseFillinCheck (self, currNode, leftRepeat, rightRepeat, nextNode, currPath, paths):
        """Repeat Fillin helper for conservativeXcoordUntangle"""
        pass


    def conservativeXcoordUntangle(self, simpleRepeatFill=False):
        """From a simplified graph return a new graph resolving simple 
        tangles (such as sequenceLinks scaffolded by bambus links)"""
        paths = []
        edges = []
        visited = set([])
        self.modifyComplexTriangleNodesInGraph ()
        self.modifySimpleTriangleNodesInGraph ()
        for subGraph in nx.weakly_connected_component_subgraphs(self.G):
            cc = subGraph.nodes()
            currPath= []
            self.sortNodesByLeftCoord(cc)
            # is the path already linear?
            linearPath = ScaffoldUtils(subGraph).getLinearPaths()[0]
            if ( len(linearPath) == len(cc) ):
                visited |= set( cc )
                paths.append( linearPath )
                continue
            for i in range (len(cc)):
                currNode = cc[i]
                # Walk back to the beginning of the path
                currNode = self._walkBackLinearPath(currNode, visited, [])                
                if (currNode in visited):
                    continue
                # break if there are no unvisited nodes
                currPath = []
                while (currNode != None and len(set(self.G.successors(currNode)) - visited ) > 0):
                    nextNodes = self.getNextNodePossibilities(currNode, visited)
                    if (None in nextNodes):
                        visited.add(currNode)
                        currPath.append(currNode)
                        paths.extend([currPath])
                        currNode = None
                    elif (len(nextNodes) == 1):
                        visited.add(currNode)
                        currPath.append(currNode)
                        currNode = nextNodes[0]                        
                    else:
                        pathTuples = self._linearizeLengthConsistentNodes (currNode, nextNodes, currPath, visited)
                        # make sure the input and output nodes are the same
                        inpSet = set(currPath) 
                        retSet = set([ node for path in pathTuples for node in path ]) 
                        assert inpSet == retSet, "%s != %s" % (inpSet, retSet)
                        if (len(pathTuples) > 1):
                            if (currNode in pathTuples[0]):
                                paths.extend([pathTuples[0]])
                                for node in pathTuples[0]: visited.add(node)
                            else:
                                visited.add(currNode)
                                currPath.append(currNode)
                                paths.extend([currPath])
                            currNode = None
                        else:    
                            # update the current node and path
                            currPath = pathTuples[-1][0:-1]
                            for node in currPath: visited.add(node)
                            currNode = pathTuples[-1][-1]
                if (currNode != None):
                    if (not(currNode in visited)):
                        visited.add(currNode)
                        currPath.append(currNode)
                    paths.append(currPath)
        return self.convertPathsToMultiGraph(paths)
    
    def _getEdgesBetweenNodePair (self, n1, n2):
        """Return all edges between nodes n1 and n2
        if no edges exist return and empty list"""
        edges = []
        if (n2 in self.G[n1]):
            edgeDict = self.G[n1][n2]
            for edgeAttr in edgeDict.values():
                edges.append([n1, n2, edgeAttr])
        return edges

    def _getNodeLength(self, n):
        return int(self.G.node[n]["length"])

    def _nodeOnPositiveStrand(self, n):
        return self.G.node[n]["strand"] == "+"

    def _vanillaSummaryLink(self, seq1, seq2, min, max):
        """Creates a SummaryLink with default values for everything but seq1, seq2, min and max"""
        return [seq1, seq2, {   "seq1" : seq1, "seq2" : seq2,   \
                                "min"  : min, "max"  : max,     \
                                "seq1_strand" : "+",            \
                                "seq2_strand" : "+",            \
                                "used" : "True", "type" : "SummaryLink",    \
                                "weight" : 1,  "validity" : VALIDITY_VALID }\
               ]

    def _fillEdgeMissingWithPredecessor (self, nA, nB, nC):
        """Return a new edge between nB->nC when the following edges exist:
        nA->nB, nA->nC."""
        
        edgesAB = self.G[nA][nB].values()
        edgesAC = self.G[nA][nC].values()
        lengthA = self._getNodeLength(nA)
        lengthB = self._getNodeLength(nB)
        lengthC = self._getNodeLength(nC)
        minAB = self.edgeBound.getMinSpan(edgesAB, nA, nB)
        maxAB = self.edgeBound.getMaxSpan(edgesAB, nA, nB)
        minAC = self.edgeBound.getMinSpan(edgesAC, nA, nC)
        maxAC = self.edgeBound.getMaxSpan(edgesAC, nA, nC)
        if (minAC > lengthB + maxAB):
            minBC = minAC - lengthB - maxAB
        else:
            minBC = 0
        maxBC = maxAC-lengthB-max(minAB,0)
        return self._vanillaSummaryLink(nB, nC, minBC, maxBC)

    def _fillEdgeMissingWithSuccessor (self, nA, nB, nC):
        """Return a new edge between nA->nB when the following edges exist:
        nA->nC, nB->nC."""
        edgesBC = self.G[nB][nC].values()
        edgesAC = self.G[nA][nC].values()
        lengthA = self._getNodeLength(nA)
        lengthB = self._getNodeLength(nB)
        lengthC = self._getNodeLength(nC)
        minBC= self.edgeBound.getMinSpan(edgesBC, nB, nC)
        maxBC= self.edgeBound.getMaxSpan(edgesBC, nB, nC)
        minAC = self.edgeBound.getMinSpan(edgesAC, nA, nC)
        maxAC = self.edgeBound.getMaxSpan(edgesAC, nA, nC)
        if (minAC > lengthB + maxBC):
            minAB = minAC - lengthB - maxBC
        else:
            minAB = 0
        maxAB = maxAC-lengthB-max(minBC,0)
        return self._vanillaSummaryLink(nA, nB, minAB, maxAB)

    def _makeSubGraphFromNodes (self, nodes):
        """Return a new graph with the input nodes.  It 
        maintains all edges within "nodes" from the current
        graph"""
        subGraph = nx.MultiDiGraph()
        for node in nodes:
            newAttrs = self.G.node[node]
            subGraph.add_node( node, attr_dict=newAttrs )
        for edge in self.G.edges(data=True):
            if (edge[0] in nodes and edge[1] in nodes):
                subGraph.add_edge(edge[0], edge[1], attr_dict=edge[2])
        return subGraph
    
    def _linearizeLengthConsistentNodes (self, currNode, nextNodes, currPath, visited):
        """From a local graph region return a [path] resolving simple 
        tangles (such as sequenceLinks scaffolded by bambus links) or
        intertwined,unconnected nodes (i.e edges: (a,b), (a,c)).  
        If a single valid path is not available return multiple paths, 
        where the paths are sorted from lowest scaffold graph position."""
        # Get all nodes sorted by xcoord
        nextNodes.append(currNode)
        self.sortNodesByLeftCoord (nextNodes)
        # Check and see if all the nodes fit between the currNode and the last node
        # extend current path with the xcoord sorted nodes

        if (nextNodes[0] == currNode):
            if (nextNodes[-1] in self.G[currNode]):
                sumLengths = self.sumNodeLengths(nextNodes[1:-1])
                n1 = currNode
                n2 = nextNodes[-1]
                edges = self.G[n1][n2].values()
                if (self.edgeBound.checkLengthValidity(edges, n1, n2, sumLengths-1, sumLengths+1) > 0):
                    currPath.extend(nextNodes)
                    return [currPath]

        # Otherwise Split the graph into linear paths (preserve out current path!)
        outPaths = [currPath]
        subGraphUtils = ScaffoldUtils(self._makeSubGraphFromNodes(nextNodes))
        linearPaths = subGraphUtils.getLinearPaths()

        # if its a simple path at the end, return two paths, currentPath and this new path
        if (len(linearPaths) == 1):
            outPaths.extend( linearPaths ) 
        else:
            # break out each separate path as a list and return the path 
            # terminating with the greatest xcoord 
            # as the last path for continued extension
            frontPath = None
            otherPaths = []
            for path in linearPaths:
                if (currNode in path):
                    frontPath = path
                else:
                    otherPaths.append(path)
            self._sortPathsByLeftCoord(otherPaths)
            allPaths = [frontPath]
            allPaths.extend(otherPaths)
            outPaths.extend(allPaths)
        return outPaths


    def _resolveRepeatTangles (self):
        """From a simplified graph return a new graph resolving simple 
        tangles (such as sequenceLinks scaffolded by bambus links)"""
        paths = self.getLinearPaths()
        for edge in self.G.edges(data = True):
            u,v,d = edge
            d
            pass
    
    
        
        
    def enumerateAllPaths (self, u,v, allowedPaths):
        paths = []
        pass

    def _validSequenceLinkLengthConstraints (self, linksGraph):
        """Returns a new graph with broken scaffold links in cases 
        where a Sequence Link Size is invalidated
        i.e. a graph with edge (a,b) (b,c) and maxSpan of 10 where there is a sequence link
        (a,r).  If r is 20 then return a graph only with edges (b,c) [and a singleton node a]"""
        pass
    

    def _splitRepeatAndOtherNodesFromEdge(self, G, edge):
        n1 = edge[0]
        n2 = edge[1]
        repeatNodes = False
        otherNodes = False
        if (G.node[n1]["repeat"]):
            repeatNode = n1            
        else:
            otherNode = n1
        if (G.node[n2]["repeat"]):
            repeatNode = n2            
        else:
            otherNode = n2
        return repeatNode, otherNode

    def breakRepeatLengthInvalidatedLinks (self, linksGraph):
        for edge in linksGraph:
            repeatNode, otherNode = self._splitRepeatAndOtherNodesFromEdge(linkGraph, edge)
            if (repeatNode and otherNode):
                successors = self.G
            pass

                
    def IsValidLengthFit (self, nA, nB, nBetweenLength):
        """Returns a boolean indicating if a length can fit 
        between two other contigs"""
        nA = False
        nB = False
        if (nB in self.G[nA]):
            n1 = nA
            n2 = nB
        elif (nA in self.G[nB]):
            n1 = nB
            n2 = nA
        if (nA and nB):
            edges = self.G[n1][n2].values()
            maxSpan = self.edgeBound.getMaxSpan(edges, n1, n2)
            if (nBetweenLength <= maxSpan):
                return True
        return False

    def addLinksNaive(self, linksGraph):	
        """Adds links to graph, usually from repeat screening, returning the modified graph.
        Note that this currently attempts to maintain linearity, by deleting redundant links.
        i.e. a graph with edges A->B->C, A->C will have A->C deleted."""

        DEFAULT_XCOORD = "-1" 
        # iterate through the edges in the linksGraph
        newEdges = []
        for linkSource, linkTarget, linkAttrs in linksGraph.edges(data=True):

            # find if either the source or the target of the node is present in the self.G
            # if both nodes are present in the self.G, something is wrong
            # check the node strand and whether it is a source or target in the links graph
            if self.G.node.has_key(linkSource) and self.G.node.has_key(linkTarget):
                logging.warning("Weird that self.G has both source %s and target %s for a link" % (linkSource, linkTarget))
                newAttrs = linksGraph.node[linkPartner]
                newAttrs["strand"] = nodeStrand
                self.G.add_node( linkPartner, attr_dict=newAttrs)
                continue
            if self.G.node.has_key(linkSource):
                isLinkSource = True 
                nodeId = linkSource
                linkPartner = linkTarget
            elif self.G.node.has_key(linkTarget):
                isLinkSource = False
                nodeId = linkTarget
                linkPartner = linkSource
            else:
                logging.error("Weird that self.G has neither source %s nor target %s for a link" % (linkSource, linkTarget))
                continue

            nodeStrand = self.G.node[nodeId]["strand"] 

            # add the linkPartner as node
            newAttrs = linksGraph.node[linkPartner]
            newAttrs["strand"] = nodeStrand
            newAttrs["xcoord"] = DEFAULT_XCOORD # Ali add in a new xcoord here
            self.G.add_node( linkPartner, attr_dict=newAttrs)

            linkPartnerLength = int( linksGraph.node[ linkPartner ]["length"] )

            #    if the node is positive and a link source
            # or if the node is negative and a link target
            if (nodeStrand == "+" and isLinkSource)       or (nodeStrand == "-" and not isLinkSource):

                # add edges to the "right" of the node in the self.G - add an edge from the node to its link partner, 
                newEdges.append( (nodeId, linkPartner, linkAttrs) )
                self.G.node[linkPartner]["xcoord"] = int(self.G.node[nodeId]["xcoord"]) + int(self.G.node[nodeId]["length"])

            #    if the node is positive and a link target
            # or if the node is negative and a link source
            elif (nodeStrand == "+" and not isLinkSource) or (nodeStrand == "-" and isLinkSource):

                # add edges to the "left" of the node in the self.G - add an edge from the link partner to the node 
                newEdges.append( (linkPartner, nodeId, linkAttrs) )
                self.G.node[linkPartner]["xcoord"] = int(self.G.node[nodeId]["xcoord"]) - int(self.G.node[linkPartner]["length"])

            else:
                logging.error("SW error: this should never occur!")

        # actually add in the edges
        for edge in newEdges: 
            self.G.add_edge( edge[0], edge[1], attr_dict = edge[2] )
        return self.G

    def isRepeatNode(self, node):
        return self.G.node.has_key(node) and self.G.node[node]["repeat"] == "true"

    def addRepeatLinks(self, linksGraph, allowDuplicateNodePairs=False):	
        """Adds links to graph, usually from repeat screening, returning the modified graph.
        Note that this currently attempts to maintain linearity, by deleting redundant links.
        i.e. a graph with edges A->B->C, A->C will have A->C deleted."""

        DEFAULT_XCOORD = "-1" 
        MIN_MAX_DELTA = 50
        # iterate through the edges in the linksGraph
        newEdges = []
        edgesToDelete = {}
        for linkSource, linkTarget, linkAttrs in linksGraph.edges(data=True):

            # find if either the source or the target of the link is a non-repeat
            # if both nodes are repeats in the self.G, something is wrong
            # check the non-repeat node strand and whether it is a source or target in the links graph
            if self.G.node.has_key(linkSource) and self.G.node.has_key(linkTarget):
                logging.warning("Weird that self.G has both source %s and target %s for a link" % (linkSource, linkTarget))
                newAttrs = linksGraph.node[linkPartner]
                newAttrs["strand"] = nodeStrand
                self.G.add_node( linkPartner, attr_dict=newAttrs)
                continue
            # note that the node is always on the positive strand. the linkPartner can be negative or positive.
            if self.G.node.has_key(linkSource):
                isLinkSource = True 
                nodeId = linkSource
                linkPartner = linkTarget
            elif self.G.node.has_key(linkTarget):
                isLinkSource = False
                nodeId = linkTarget
                linkPartner = linkSource
            else:
                logging.error("Weird that self.G has neither source %s nor target %s for a link" % (linkSource, linkTarget))
                continue

            nodeStrand = self.G.node[nodeId]["strand"] 
            partnerStrand = "+" if nodeStrand == linksGraph.node[linkPartner]["strand"] else "-"

            # add the linkPartner as node
            newAttrs = linksGraph.node[linkPartner]
            newAttrs["strand"] = partnerStrand
            self.G.add_node( linkPartner, attr_dict=newAttrs)

            linkPartnerLength = int( linksGraph.node[ linkPartner ]["length"] )

            #    if the node is positive and a link source
            # or if the node is negative and a link target
            if (nodeStrand == "+" and isLinkSource)       or (nodeStrand == "-" and not isLinkSource):

                successors = self.G.successors( nodeId )

                # add edges to the "right" of the node in the self.G
                # add edges from the link partner to the node's self.G successor if it exists
                assert len(successors) <= 1, \
                        "More than one successor: %s" % str(successors)

                if successors: 

                    # note the old edge for later deletion
                    edgePair = (nodeId, successors[0])
                    # continue if we had already seen this edgePair and we do not allow duplicateNodePairs
                    if edgesToDelete.has_key( edgePair ) and not allowDuplicateNodePairs: continue
                    edgesToDelete[ edgePair ] = 1

                    edgeData = deepcopy( self.G.get_edge_data( nodeId, successors[0] ) )
                    assert len(edgeData.keys()) == 1, \
                            "Strange edgeData: %s " % str(edgeData.keys())
                    # correct for link length
                    edgeData[0]["max"] = str( int(edgeData[0]["max"]) - linkPartnerLength)
                    if int(edgeData[0]["min"]) > int(edgeData[0]["max"]): # ensure min is still less than max
                        edgeData[0]["min"] = str(int(edgeData[0]["max"]) - MIN_MAX_DELTA) 
                    newEdges.append( (linkPartner, successors[0], edgeData[0]) )

                # add an edge from the node to its link partner, 
                newEdges.append( (nodeId, linkPartner, linkAttrs) )

            #    if the node is positive and a link target
            # or if the node is negative and a link source
            elif (nodeStrand == "+" and not isLinkSource) or (nodeStrand == "-" and isLinkSource):

                # add edges to the "left" of the node in the self.G
                # add an edge from the node's self.G predecessor to its link partner if the predecessor exists
                predecessors = self.G.predecessors( nodeId )
                assert len(predecessors) <= 1, \
                        "More than one predecessor: %s" % str(predecessors)
                if predecessors: 

                    # note the old edge for later deletion
                    edgePair = (predecessors[0], nodeId)
                    # continue if we had already seen this edgePair and we do not allow duplicateNodePairs
                    if edgesToDelete.has_key( edgePair ) and not allowDuplicateNodePairs: continue
                    edgesToDelete[ edgePair ] = 1

                    edgeData = deepcopy( self.G.get_edge_data( predecessors[0], nodeId ) )
                    assert len(edgeData.keys()) == 1,\
                            "Strange edgeData: %s " % str(edgeData.keys())
                    edgeData[0]["max"] = str( int(edgeData[0]["max"]) - linkPartnerLength )
                    if int(edgeData[0]["min"]) > int(edgeData[0]["max"]): # ensure min is still less than max
                        edgeData[0]["min"] = str(int(edgeData[0]["max"]) - MIN_MAX_DELTA) 
                    newEdges.append( (predecessors[0], linkPartner, edgeData[0]) )

                # and from the link partner to the node 
                newEdges.append( (linkPartner, nodeId, linkAttrs) )

            else:
                logging.error("SW error: this should never occur!")

        # actually add in the edges
        # only add 
        for edge in newEdges: 
            self.G.add_edge( edge[0], edge[1], attr_dict = edge[2] )
        for edge in edgesToDelete.keys():
            self.G.remove_edge( edge[0], edge[1] )
        return self.G

    def edgeSeq(self, node1, node2):
        """Assumes only one edge between node1 and node2 and returns the N's that should be between them."""
        edgeDict = self.G.get_edge_data( node1, node2 )[0]
        if edgeDict.has_key("mid"):
            numNs = int(edgeDict["mid"])
        else:
            numNs = (int(edgeDict["max"]) + int(edgeDict["min"]) ) / 2
        if numNs < 1:
            logging.warning("Found a non-positive number of N's = %i for edge %s - %s. Putting in one N." % (numNs, node1, node2))
            numNs = 1
        return "N" * numNs
    
    def _checkEdgesForType (self, n1, n2, edgeType):
        """Checks to see if an edge of type=edgeType exists.  
        Returns the first such edge it detects or False otherwise"""
        edges = self.G[n1][n2].values()
        for edge in edges:
            if ("type" in edge and edge["type"] == edgeType):
                return edge
        return False
 
    def filterEdges(self, redundancy=1, absoluteMinSpan=-500):
        """ Build a new graph in which a single edge summarizes each edge
        in the original multiDigraph.  The type of edge is summaryLink"""
        filteredG = nx.MultiDiGraph()
        filteredG.add_nodes_from(self.G.nodes(data=True))
        #iterate over nodes instead of edges to graph all edges incident on n1,n2
        for n1 in self.G.nodes_iter():
            for n2 in self.G.successors_iter(n1):
                edgeAttrs = self.G[n1][n2].values()                
                if len(edgeAttrs) < redundancy: continue
                # an absoluteMinSpan allows us to include UNSEEN edges
                maxSpan = self.edgeBound.getMaxSpan(edgeAttrs, n1, n2)
                if maxSpan < absoluteMinSpan: continue
                for edgeAttr in edgeAttrs: 
                    filteredG.add_edge(n1, n2, attr_dict=edgeAttr)

        return filteredG


    def buildSummaryGraph (self):
        """ Build a new graph in which a single edge summarizes each edge
        in the original multiDigraph.  The type of edge is summaryLink"""
        summaryG = nx.MultiDiGraph()
        summaryG.add_nodes_from(self.G.nodes(data=True))
        #iterate over nodes instead of edges to graph all edges incident on n1,n2
        for n1 in self.G.nodes_iter():
            for n2 in self.G.successors_iter(n1):
                if (self._checkEdgesForType(n1,n2,"SummaryLink")):
                    summaryG.add_edge(n1,n2,attr_dict=self._checkEdgesForType(n1,n2,"SummaryLink"))
                else:
                    edges = self.G[n1][n2].values()                
                    #summarize Spans on edges between n1 -> n2
                    minSpan = self.edgeBound.getMinSpan(edges, n1, n2)
                    maxSpan = self.edgeBound.getMaxSpan(edges, n1, n2)
                    summaryG.add_edge(n1,n2, \
                                      seq1_strand = edges[0]["seq1_strand"],  \
                                      seq2_strand = edges[0]["seq2_strand"],  \
                                      min=str(minSpan), max=str(maxSpan), weight=str(len(edges)), \
                                      type="SummaryLink", used=str(True), validity=VALIDITY_VALID)
        return summaryG

    def mergeGraph(self, graph):
         """ Merge input graph with self.G"""
         # Need to do some error checking!
         self.G.add_edges_from(graph.edges(data=True))
        

    def dijkstraWithNodeLengthsMeanSpans (self, source):
        """For use in getting shortest distances in repeat graph from one
        node to any other node """
        # Initialize distances as infinity
        dist = {}
        previous = {}
        for v in self.G:
            dist[v] = float("inf")
            prev[v] = None
            
        # Initialize current node as distance 0 (distance to itself)
        dist[source] = 0
        
        # Create the set of remaining nodes
        R = set(self.G.nodes())
        
        # Update path distances
        while (len(R) > 0):
            minV = None
            minDist = float("inf")
            for v in R:
                if (dist[v] < minDist):
                    minV = v
                    minDist = dist[v]
            R.remove(minV)
            for u in self.G.successors(minV):
                altDist = dist[minV] + self._getNodeLength(minV) + self._getMidSpan(minV, u)
                if (altDist < dist[u]):
                    dist[u] = altDist
                    prev[u] = minV
        return dist, previous
            
    def enumeratePathsBetweenNodes (self, source, sink, maxSpan):
        """Return the set of all paths between two nodes in which
        the source and sink nodes are less than maxSpan apart"""
        # Initialize Vertex Stack:
        pathStack = [(source, [], 0)]        
        # Initialize output paths
        validPaths = []

        # Enumerate all paths up to maxSpan
        while (len(pathStack) > 0):
            v, p, d = pathStack.pop()
            for u, edgeDict in self.G[v].items():
                vuSpan = self._getMinSpan(v, u)
                # Check to see if the currentPath satisfies criteria
                if (u == sink):
                    currSinkSpan = d + vuSpan + vuSpan
                    if (currSinkSpan < maxSpan):
                        validPaths.append(p)
                elif (u in p):
                    # Avoids cycles
                    continue
                else:
                    # Add the current path to the stack if it is < maxSpan
                    uDist = self._getNodeLength(u) + d + vuSpan
                    if (uDist < maxSpan):
                        uPath = deepcopy(p)
                        uPath.append(u)
                        pathStack.append((u, uPath, uDist))
        
        return validPaths

    def subsetPathsContainAllNodes (self, paths, nodes):
        validPaths = []
        for path in paths:
            nodeFlag = True
            if (False in map(lambda node: node in paths, nodes)):
                continue
            else:
                validPaths.append(path)
        return validPaths
    
    def convertEdgeInRepeatGraph (self, n1, n2, maxSpan, repeatGraph):
        """Check to see if suggested edge is valid considering the 
        deBruijn input"""
        repeatScaff = ScaffoldUtils(repeatGraph)
        rGraphDists, rGraphPrev = self.dijkstraWithNodeLengthsMeanSpans (n1)
        

    def _getMidSpan (self, n1, n2):
        edges = self.G[n1][n2].values()
        return self.edgeBound.getMidSpan(edges, n1, n2)

    def _getMinSpan (self, n1, n2):
        edges = self.G[n1][n2].values()
        return self.edgeBound.getMinSpan(edges, n1, n2)

    def _getMaxSpan (self, n1, n2):
        edges = self.G[n1][n2].values()
        return self.edgeBound.getMaxSpan(edges, n1, n2)
    
    def numUsedValidEdges(self):
        return sum( map( lambda x: int(x[2]["weight"]), \
                buildMultiDiGraphFromTypedEdges(self.G).edges(data=True)))

class IEdgeBound:
    """Interface for business logic produces span bounds on edges"""
    def getMinSpan(self, edges, n1, n2):
        """ Returns a lower bound on edge length"""
        pass
        
    def getMidSpan(self, edges, n1, n2):
        """ Returns a mid bound on edge length"""
        pass
        
    def getMaxSpan(self, edges, n1, n2):
        """ Returns an upper bound on edge length"""
        pass
   
class StrictEdgeBound (IEdgeBound):
    """Edge bound policy which strictly enforces the min max bounds
    on each individual edge when making a summary edge"""
    def __init__ (self, G):
        self.G = G
                        
    def _summarizeSpanType (self, edges, n1, n2, spanType):
        values = numpy.zeros(len(edges))
        for i, e in enumerate(edges):
            summarySpan = float(e[spanType])
            if ("type" in e and e["type"] == "SummaryLink"):
                values[i] = summarySpan
                continue
            if (self.G.node[n1]["strand"] == "+"):
                n1overlap = float(self.G.node[n1]["length"]) - float(e["seq1_start"])
            else:
                n1overlap =  float(e["seq1_start"])
            if (self.G.node[n2]["strand"] == "-"):
                n2overlap = float(self.G.node[n2]["length"]) - float(e["seq2_start"])
            else:
                n2overlap =  float(e["seq2_start"])
            values[i] = summarySpan - n1overlap - n2overlap
        return values

    def getMidSpan (self, edges, n1, n2):
        minSpans =self._summarizeSpanType (edges, n1, n2, "min")
        maxSpans =self._summarizeSpanType (edges, n1, n2, "max")
        midSpans = numpy.zeros(len(minSpans))
        for i in range (len(minSpans)):
            midSpans[i] = (maxSpans[i]-minSpans[i])/2.0
        return numpy.median(midSpans)
            
    def getMinSpan (self, edges, n1, n2):
        minSpans =self._summarizeSpanType (edges, n1, n2, "min")
        return int(numpy.max(minSpans))
        
    def getMaxSpan (self, edges, n1, n2):
        maxSpans =self._summarizeSpanType (edges, n1, n2, "max")
        return int(numpy.min(maxSpans))

class MedianEdgeBound (IEdgeBound):
    """Edge bound policy which strictly enforces the min max bounds
    on each individual edge when making a summary edge"""
    
    def __init__ (self, G):
        self.G = G
                        
    def _summarizeSpanType (self, edges, n1, n2, spanType):
        values = numpy.zeros(len(edges))
        for i, e in enumerate(edges):
            summarySpan = float(e[spanType])
            if ("type" in e and e["type"] == "SummaryLink"):
                values[i] = summarySpan
                continue            
            if (self.G.node[n1]["strand"] == "+"):
                n1overlap = float(self.G.node[n1]["length"]) - float(e["seq1_start"])
            else:
                n1overlap =  float(e["seq1_stop"])
            if (self.G.node[n2]["strand"] == "+"):
                n2overlap =  float(e["seq2_stop"])
            else:
                n2overlap = float(self.G.node[n2]["length"]) - float(e["seq2_start"])
            values[i] = summarySpan - n1overlap - n2overlap
        return values

    def getMidSpan (self, edges, n1, n2):
        minSpans =self._summarizeSpanType (edges, n1, n2, "min")
        maxSpans =self._summarizeSpanType (edges, n1, n2, "max")
        midSpans = numpy.zeros(len(minSpans))
        for i in range (len(minSpans)):
            midSpans[i] = (maxSpans[i]-minSpans[i])/2.0
        return numpy.median(midSpans)
        
    def getMinSpan (self, edges, n1, n2):
        minSpans =self._summarizeSpanType (edges, n1, n2, "min")
        return int(numpy.median(minSpans))

    def getMaxSpan (self, edges, n1, n2):
        maxSpans =self._summarizeSpanType (edges, n1, n2, "max")
        return int(numpy.median(maxSpans))

    def getMaximumValueInRange (self, rangei, meanValue, minValue, maxValue):
        if (rangei.contains(meanValue)):
            return 1
        else:
            if (rangei.getEnd()-meanValue > 0):
                dist = float(rangei.getStart()-meanValue)
                return dist/float(maxValue-meanValue)
            else:
                dist = meanValue - rangei.getEnd()
                return dist/float(meanValue-minValue)

    def checkLengthValidity (self, edges, n1, n2, betweenLengthMin, betweenLengthMax):
        minSpan = self.getMinSpan(edges, n1, n2)
        maxSpan = self.getMaxSpan(edges, n1, n2)
        meanSpan = self.getMidSpan(edges, n1, n2)
        edgeRange = Range.Range(minSpan, maxSpan)
        lengthRange = Range.Range(betweenLengthMin, betweenLengthMax)
        if (edgeRange.intersects(lengthRange)):
            intersectRange = edgeRange.intersect(lengthRange)
            return self.getMaximumValueInRange(intersectRange, meanSpan, minSpan, maxSpan)
        return 0

class MyerEdgeBound (IEdgeBound):
    """Edge bound policy which uses the Eugene Myers summary edge:
    the mean +/- standard deviation of the mean for each edge"""

    def __init__ (self, stdCut=3.0):
        self.stdCut = stdCut
            
    def _getEdgeDist (self, edges):
        midPointEdges = []
        for edge in edges:
            midPointEdges.append((float(edge["min"])+float(edge["max"]))/2.0)
        return numpy.mean(midPointEdges), numpy.std(midPointEdges)
        
    def getMidSpan (self, edges):
        meanSpan, stdSpan = self._getEdgeDist(edges)
        return int(meanSpan)
        
    def getMinSpan (self, edges):
        meanSpan, stdSpan = self._getEdgeDist(edges)
        return int(max(0, meanSpan-self.stdCut*stdSpan))
        
    def getMaxSpan (self, edges):
        meanSpan, stdSpan = self._getEdgeDist(edges)
        return int(meanSpan+self.stdCut*stdSpan)
    
            


def buildMultiDiGraphFromTypedEdges (G, validityTypes = [VALIDITY_VALID], allowUnused=False):
    subsetG = nx.MultiDiGraph()
    subsetG.add_nodes_from(G.nodes(data=True))
    for edge in G.edges_iter(data=True):
        if ((edge[2]["used"] == "True" or allowUnused) and edge[2]["validity"] in validityTypes):
            subsetG.add_edge (edge[0], edge[1], attr_dict=edge[2])
    return subsetG


            
        
if __name__=='__main__':
    G = nx.read_graphml(sys.argv[1])
    subsetG = buildMultiDiGraphFromTypedEdges(G)
    su = ScaffoldUtils(subsetG)
    summaryG = su.buildSummaryGraph()
    nx.write_graphml(su.G, "origina_tmp.gml")
    nx.write_graphml(summaryG, "summary.gml")
    su.mergeGraph(summaryG)
    nx.write_graphml(su.G, "orig_and_summ.gml")

