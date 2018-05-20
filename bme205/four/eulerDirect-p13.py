#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Finds Eulerian path of given directed nodes 
Input: A directed path containing a Eulerian path
Output: The Eulerian path 
Algorithm from: http://www.graph-magics.com/articles/euler.php
"""


import sys
from collections import defaultdict


class EulerianPath:
    """ Constructs the de Bruijn graph of a string for kMer length -k """
    def __init__(self, dirGraph):
        self.dirGraph = dirGraph
        self.stack = []

    def getEdges(self):
        """ Counts incoming and outgoing edges of each node and uses this information do determine if they make an Eulerian path or an Eulerian cycle. If it is an Eulerian path then the starting node is determined """
        inCount = {}
        outCount = {}
        edgesDict = {}
        maxOut = totalIn = totalOut = 0

        for each in self.dirGraph:

            outNode = each.split('-')[0].rstrip()
            inNode = each.split('>')[1].strip()
            # stores counts of incoming and outgoing edges for each node          
            if outNode not in outCount.keys(): outCount[outNode] = 0
            if outNode not in inCount.keys(): inCount[outNode] = 0
            # if the node has multiple outgoing edges
            if ',' in inNode:
                nodesIn = inNode.split(',')
                edgesDict[outNode] = nodesIn
                for nodes in nodesIn:
                    if nodes not in inCount.keys(): inCount[nodes] = 1
                    else: inCount[nodes] += 1
                    if nodes not in outCount.keys(): outCount[nodes] = 0
                    outCount[outNode] += 1
            # if there is only one outgoing edge for that node
            else:
                outCount[outNode] += 1
                edgesDict[outNode] = [inNode]
                if inNode not in inCount.keys(): inCount[inNode] = 1
                else: inCount[inNode] += 1
                if inNode not in outCount.keys(): outCount[inNode] = 0
        # finds starting node if Eulerian path
        for key in outCount.keys():
            if outCount[key] > inCount[key]:
                totalOut += 1
                curIn = key
                print(curIn)
            if inCount[key] > outCount[key]:
                totalIn += 1
                print(key)
        print(totalIn)
        print(totalOut)
        # chooses a random node in the dict if it is an Eulerian cycle
        if totalOut == 0 and totalIn == 0:
            curIn = outCount[0]
        # if is not a valid Eulerian cycle or path
        elif totalOut != totalIn:
            print('No Eulerian path or cycle is possible for the given kMers')
            sys.exit(0)
        return edgesDict, curIn

    def traversePath(self, edgeDict, firstIn):
        """ Takes the dictionary of available edges and the defined starting node and uses a stack to store traversed nodes while removing edges from the dictionary of available edges. A node is appended to the final path when it has no outgoing edges """
        curNode = firstIn
        finalPath = []
        # while stack is not empty or dictionary of edges is not empty
        while len(self.stack) != 0 or len(edgeDict[curNode]) != 0:
            if curNode not in edgeDict.keys() or not edgeDict[curNode]:
                finalPath.append(curNode)
                curNode = self.stack.pop()
                continue
            self.stack.append(curNode)
            prevNode = curNode
            curNode = edgeDict[curNode][0]
            del edgeDict[prevNode][0]
        finalPath.append(curNode)
        return finalPath

def main():

    dirGraph = sys.stdin.readlines()
    graphKmers = EulerianPath(dirGraph)
    edgeDict, firstIn = graphKmers.getEdges()
    graphList = graphKmers.traversePath(edgeDict, firstIn)

    graphList.reverse()
    print('->'.join(graphList))
    

if __name__== "__main__":
    main()
