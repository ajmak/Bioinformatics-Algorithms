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
import pdb

class DirectedAcyclicGraph:
    """ """
    def __init__(self, kMers):
        self.kMers = kMers


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

    def longestPath(self, topOrder, edgeDict, backDict, source):
        """Finds the longest path in a topologically ordered graph given the source of the subgraph. Returns a dictionary of cumulative weights for a path up to each node and a dictionary of back pointers """
        pointBack = {}
        weightDict = {key: -math.inf for key in topOrder}
        weightDict[source] = 0

        for curNode in topOrder:
            # this ensures that no nodes are considered until the source node and only the children in its path will be considered
            if weightDict[curNode] == -math.inf: 
                continue
            if curNode in edgeDict.keys():
                for child in edgeDict[curNode]:
                    # finds max weighted path to that child node
                    if weightDict[child[0]] < weightDict[curNode] + int(child[1]):
                        weightDict[child[0]] = weightDict[curNode] + int(child[1])
                        pointBack[child[0]] = curNode
        print(weightDict)
        print(pointBack)
        return weightDict, pointBack


    def backOut(self, weightDict, pointBack, source, sink):
        """Creates a list of the longest path from the source to the sink given a dictionary of back pointers that outline the route taken to each node and a dictionary of cumulative weights up until each node """
        curNode = sink
        outList = []
        while curNode != source:
            outList.append(curNode)
            curNode = pointBack[curNode]
        
        outList.append(source)
        return weightDict[sink], outList
        
def main():

    nodes = sys.stdin.read().splitlines()
    source = nodes.pop(0)
    sink = nodes.pop(0)
    graphKmers = DirectedAcyclicGraph(nodes)

    edgeDict, backDict, candidates = graphKmers.getEdges(nodes)

    topologicalOrder = graphKmers.topologicalSort(edgeDict, backDict, candidates)

    weightDict, pointBack = graphKmers.longestPath(topologicalOrder, edgeDict, backDict, source)

    finalWeight, finalPath = graphKmers.backOut(weightDict, pointBack, source, sink)

    print(finalWeight)
    finalPath.reverse()
    print('->'.join(finalPath))
    

if __name__== "__main__":
    main()
