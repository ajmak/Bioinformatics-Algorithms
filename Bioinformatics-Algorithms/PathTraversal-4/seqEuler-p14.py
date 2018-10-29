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
from collections import defaultdict

class EulerianPath:
    """ Constructs the de Bruijn graph of a string for kMer length -k """
    def __init__(self, kMers):
        self.kMers = kMers
        self.stack = []

    def buildaBruijn(self):
        """ Finds and organizes all kMers in a list into an adjacency list and returns as a sorted list """
        nodesDict = {}
        for x in self.kMers:
            pre = x[:-1].strip()
            end = x[1:].strip()
            if pre in nodesDict.keys():
                nodesDict[pre].append(end)
            else:
                nodesDict[pre] = [end]
        # sorts lists within dictionary for nodes with multiple edges 
        sortnodesDict = {x:sorted(nodesDict[x]) for x in nodesDict.keys()}
        # sorts output list based on prefix node                   
        orderList = sorted(sortnodesDict.items())
        return orderList, sortnodesDict


    def getEdges(self, deBruijn):
        """ Counts incoming and outgoing edges of each node and uses this information do determine if they make an Eulerian path or an Eulerian cycle. If it is an Eulerian path then the starting node is determined """
        inCount = {}
        outCount = {}
        edgesDict = {}
        maxOut = totalIn = totalOut = 0

        for each in deBruijn.keys():
            outNode = each
            inNode = deBruijn[each]
            # stores incoming and outgoing edges for each node
            if outNode not in outCount.keys(): outCount[outNode] = 0
            if outNode not in inCount.keys(): inCount[outNode] = 0
            # if the node has multiple outgoing edges
            if len(inNode) > 1:
                edgesDict[outNode] = inNode
                for nodes in inNode:
                    if nodes not in inCount.keys(): inCount[nodes] = 1
                    else: inCount[nodes] += 1
                    if nodes not in outCount.keys(): outCount[nodes] = 0
                    outCount[outNode] += 1
            # if there is only one outgoing edge for that node
            else:
                inNode = inNode[0]
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
            if inCount[key] > outCount[key]:
                totalIn += 1
        
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
            # if doesn't have any neighbors left
            if curNode not in edgeDict.keys() or not edgeDict[curNode]:
                finalPath.append(curNode)
                # current node becomes the node on the top of the stack
                curNode = self.stack.pop()
                continue
            self.stack.append(curNode)
            prevNode = curNode
            curNode = edgeDict[curNode][0]
            # node is deleted from the dictionary of edges
            del edgeDict[prevNode][0]
        finalPath.append(curNode)
        return self.joinPath(finalPath)

    def joinPath(self, pathList):
        """ Given a list of kMers, appends ordered genome graph into a full sequence  """
        sequence = ''
        # appends first letter of each kMer until last kMer and appends entire last kMer
#        print(pathList)

        for node in reversed(pathList):
            if node == pathList[0]:
                sequence += node
            else:
                chars = list(node)
                sequence += chars[0]
        return sequence

def main():

    kMers = sys.stdin.read().splitlines()
    graphKmers = EulerianPath(kMers)
    orderList, deBruijn = graphKmers.buildaBruijn()
    edgeDict, firstIn = graphKmers.getEdges(deBruijn)
    graphList = graphKmers.traversePath(edgeDict, firstIn)
    print(graphList)

    

if __name__== "__main__":
    main()
