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

    def longestPath(self, start, end, step, topOrder, edgeDict, emissionDict):
        """Finds the longest path in a topologically ordered graph given the source of the subgraph. Returns a dictionary of cumulative weights for a path up to each node and a dictionary of back pointers """

        pointBack = {key: {state: 'parent' for state in edgeDict.keys() if state != 'source' and state != 'sink'}for key in range(start, end, step)}
        weightDict = {key: {state: np.float64(0) for state in edgeDict.keys() if state != 'source' and state != 'sink'}for key in range(start-step, end-step, step)}
        if start == 1:
            weightDict[start-step] = {'source': np.float64(1)}
            weightDict[end-step] = {'sink': np.float64(0)}
            pointBack[end-step] = {'sink': 'parent'}
        else:
            weightDict[end-step] = {state: np.float64(1.0) for state in edgeDict.keys() if state != 'sink'}
            weightDict[start-step] = {'sink': 1}

        for i in range(start, end, step):
            for child in weightDict[i].keys():
                # finds max weighted path to that child node
                for parent in weightDict[i-step].keys():
                    
                    if weightDict[i][child] == 1.0: break
                    if child == 'sink' or parent == 'sink': 
                        nodeWeight = 1
                        weightDict[i][child] += weightDict[i-step][parent] * nodeWeight
                    else:
                        if len(topOrder) - 2 < start:
                            nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[parent + topOrder[i]]))
                            weightDict[i][child] += weightDict[i-step][parent] * nodeWeight
                        else:
                            nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[child + topOrder[i]]))
                            weightDict[i][child] += np.multiply(weightDict[i-step][parent], nodeWeight)
        if end == 0: end = 3
        sumTotal = (sum(weightDict[end-2].values()))

        return weightDict, sumTotal


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
        print(''.join(outList))
#        return weightDict[sink], outList
        
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
