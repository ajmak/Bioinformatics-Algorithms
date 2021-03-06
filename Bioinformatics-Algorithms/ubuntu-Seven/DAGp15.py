#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Yianni Anastopulos, Xian Chang, Lori Siao

"""
Finds the longest path in a directed acyclic graph from a given source node to sink node and its associated weight.

Input: The edges of a directed acyclic graph and their weights and the source and sink of the subpath desired
Output: The longest path from source to sink and the cumulative edge weight it took to get there

Parses source and sink nodes from input file when provided at the top of tbe file and separated by newlines
"""

import sys
import math
from collections import defaultdict
from random import choice
import numpy as np
import pdb

class DirectedAcyclicGraph:
    """ """
    def __init__(self):
        self.kMers = None

    def getEdges(self, nodesList):
        """Parses the input list of edges and stores them in a dictionary of tuples with their in/out nodes and associated weights. This dictionary is returned along with a list of nodes with no incoming edges """
        inList = []
        edgesDict = {}
        backDict = {}
        # candidates list is a set of all nodes with no incoming edges
        candidates = []
        
        for each in nodesList:
            nodes = each.replace('->', ' ').replace(':', ' ').split()
            inList.append(nodes[1])
            if nodes[0] in edgesDict.keys():
                edgesDict[nodes[0]].append((nodes[1], nodes[2]))
            else:
                edgesDict[nodes[0]] = [(nodes[1], nodes[2])]
            if nodes[1] in backDict.keys(): 
                backDict[nodes[1]].append((nodes[0], nodes[2]))
            else: 
                backDict[nodes[1]] = [(nodes[0], nodes[2])]
        for outgoing in set(edgesDict.keys()):
            # if has no incoming edges and is not already in candidates
            if outgoing not in inList and outgoing not in candidates: candidates.append(outgoing)
        return edgesDict, backDict, candidates

    def topologicalSort(self, edgeDict, backDict, candidates):
        """ Takes a dictionary of edges and their associated weights, the inverse of that dictionary, and the original list of nodes with no incoming edges and sorts all nodes into topological order """
        finalOrder = []
        while len(candidates) != 0:
            # treat as queue rather than randomizing since want to consider children of each node before considering grandchildren of one node
            curNode = candidates[0]
            candidates.remove(curNode)
            finalOrder.append(curNode)
            if curNode not in edgeDict.keys(): continue
            for outCurNode in edgeDict[curNode]:
                ready = None
                for backCurNode in backDict[outCurNode[0]]:
                    if backCurNode[0] not in finalOrder: ready = False
                # if all incoming nodes have been dealt with
                if ready != False:
                    for backCurNode in backDict[outCurNode[0]]:
                        if outCurNode[0] not in candidates:
                            candidates.append(outCurNode[0])
                else:
                    continue

        return(finalOrder)

    def longestPath(self, topOrder, edgeDict, transitionDict, emissionDict):
        """Finds the longest path in a topologically ordered graph given the source of the subgraph. Returns a dictionary of cumulative weights for a path up to each node and a dictionary of back pointers """
        pointBack = {key: {state: 'parent' for state in edgeDict.keys() if state != 'source' and state != 'sink'}for key in range(len(topOrder)-1)}
        weightDict = {pos: {state: float('-inf') for state in edgeDict.keys() if state != 'source' and state != 'sink'}for pos in range(1, len(topOrder)-1)}
        weightDict[0] = {'source': 1}
        weightDict[len(topOrder)-1] = {'sink': float('-inf')}
        pointBack[len(topOrder)-1] = {'sink': 'parent'}
        for i in range(1, len(topOrder)):
            for child in weightDict[i].keys():
                # finds max weighted path to that child node
                for parent in weightDict[i-1].keys():
                    if child == 'sink': 
                        nodeWeight = 1
                    else:
                        nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[child + topOrder[i]]))
                    if weightDict[i][child] < weightDict[i-1][parent] * nodeWeight:
                        weightDict[i][child] = np.multiply(weightDict[i-1][parent], nodeWeight)
                        bestParent = parent
                pointBack[i][child] = bestParent
        
        finalPath = self.backOut(weightDict, pointBack, 'source', 'sink')
        return finalPath

    def backOut(self, weightDict, pointBack, source, sink):
        """Creates a list of the longest path from the source to the sink given a dictionary of back pointers that outline the route taken to each node and a dictionary of cumulative weights up until each node """
        curNode = sink
        outList = []
        for i in range(len(weightDict.keys())-1, 0, -1):
            if curNode == 'sink' or curNode == 'source':
                curNode = pointBack[i][curNode]
                continue
            else:
                outList.append(curNode)
                curNode = pointBack[i][curNode]
        outList = reversed(outList)
        return ''.join(outList)

        
def main():

    nodes = sys.stdin.read().splitlines()
    source = nodes.pop(0)
    sink = nodes.pop(0)
    graphKmers = DirectedAcyclicGraph(nodes)

    edgeDict, backDict, candidates = graphKmers.getEdges(nodes)

    topologicalOrder = graphKmers.topologicalSort(edgeDict, backDict, candidates)

    weightDict, pointBack = graphKmers.longestPath(topologicalOrder, edgeDict, source)

    finalWeight, finalPath = graphKmers.backOut(weightDict, pointBack, source, sink)
    sys.stdin = open('/dev/tty')
    print(finalWeight)
    finalPath.reverse()
    print('->'.join(finalPath))
    

if __name__== "__main__":
    main()
