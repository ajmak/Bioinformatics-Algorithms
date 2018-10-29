#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Finds Eulerian path of given directed nodes 
Input: A directed path containing a Eulerian path
Output: The Eulerian path 
Algorithm from: http://www.graph-magics.com/articles/euler.php
Took out argparse because it was unnecessary 
"""

import sys
import math
from collections import defaultdict
from random import choice
import pdb

class DirectedAcyclicGraph:
    """ Constructs the de Bruijn graph of a string for kMer length -k """
    def __init__(self, kMers):
        self.kMers = kMers


    def getEdges(self, nodesList):
        """ Counts incoming and outgoing edges of each node and uses this information do determine if they make an Eulerian path or an Eulerian cycle. If it is an Eulerian path then the starting node is determined """
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
        print(backDict)
        return edgesDict, backDict, candidates

    def topologicalSort(self, edgeDict, backDict, candidates):
        finalOrder = []
        # while dictionary of edges is not empty

        while len(candidates) != 0:
            ready = None
            # treat as queue rather than randomizing since want to consider children of each node before considering grandchildren of one node
            curNode = candidates[0]
            candidates.remove(curNode)
            finalOrder.append(curNode)
            if curNode not in edgeDict.keys(): continue
            for outCurNode in edgeDict[curNode]:
                for backCurNode in backDict[outCurNode[0]]:
                    if backCurNode[0] not in finalOrder: ready = False
                if ready != False:
                    for backCurNode in backDict[outCurNode[0]]:
                        if outCurNode[0] not in candidates:
                            candidates.append(outCurNode[0])
                else:
                    continue

        print(finalOrder)
        return(finalOrder)

    def longestPath(self, topOrder, edgeDict, backDict, candidates, source, sink):
        #dictionary of tuples with node, its predecessor that it came from 
        pointBack = {}
        weightDict = {key: -math.inf for key in edgeDict.keys()}
        weightDict[source] = 0
            #make dictionary with cumulative score from previous node and previous node ; second dictionary for back pointers (points to winning path)
        for curNode in topOrder:
            print(curNode)
            if weightDict[curNode] == -math.inf: 
                continue
            curParents = [each for each in backDict[curNode]]
            print(curParents)
            maxCumulative, bestMom = self.findMax(curParents, weightDict) 
            weightDict[curNode] = maxCumulative
            pointBack[curNode] = bestMom
        print(weightDict)
        return weightDict, pointBack

    def findMax(self, curParents, weightDict):
        maxi = -math.inf
        maxMom = None
        for mom in curParents:
            if weightDict[mom[0]] + int(mom[1]) > maxi:
                maxi = weightDict[mom[0]] + int(mom[1])
                maxMom = mom[0]
        return maxi, maxMom

    def backOut(self, weightDict, pointBack):
        pass
def main():

    nodes = sys.stdin.read().splitlines()
    source = nodes.pop(0)
    sink = nodes.pop(0)
    sys.stdin = open('/dev/tty')
    graphKmers = DirectedAcyclicGraph(nodes)
    edgeDict, backDict, candidates = graphKmers.getEdges(nodes)
    topologicalOrder = graphKmers.topologicalSort(edgeDict, backDict, candidates)
    finalPath = graphKmers.longestPath(topologicalOrder, edgeDict, backDict, candidates, source, sink)
#    print(graphList)

    

if __name__== "__main__":
    main()
