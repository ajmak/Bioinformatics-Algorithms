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
    def __init__(self, lenStates, sequence):
        self.maxStateIndex = len(sequence)
        #int(lenStates/3 -1)
        self.lenSeq = len(sequence)

    def longestPath(self, sequence, stateList, edgeDict, emissionDict):
        """Finds the longest path in a topologically ordered graph given the source of the subgraph. Returns a dictionary of cumulative weights for a path up to each node and a dictionary of back pointers """
        pointBack = {key: {state: 'parent' for state in edgeDict.keys() if state != 'S' and state != 'E'}for key in range(len(stateList)-1)}
        weightDict = {stateIdx: {} for stateIdx in range(len(sequence)+1)}
        weightDict[0] = {'S': float(1), 'I0': float(1)}
#        weightDict[0] = {'S': 0}
#        weightDict[len(stateList)-1] = {'sink': float('-inf')}
#        pointBack[len(stateList)-1] = {'sink': 'parent'}

        curState = 'S'
        seqIndex = 0

        while curState != 'E':
            for child in edgeDict[curState].keys():
                # finds max weighted path to that child node
#                pdb.set_trace()
                childIndex = seqIndex
                if child[:-1] != 'D': childIndex = seqIndex + 1
                print(str(child) + ' child')
                print(self.getPossibleParents(child, childIndex))
                for parent in self.getPossibleParents(child, childIndex):
                    print(parent + ' parent')
                    print('******************')
                    print(str(childIndex) + ' childIndex')
                    print(sequence)
                    print(len(sequence))
                    print(sequence[childIndex-1])
                    print('******************')
                    nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[child + sequence[childIndex-1]]))
                    print(nodeWeight)
                    if child not in weightDict[childIndex].keys():
                        weightDict[childIndex][child] = float('-inf')
#                    if parent not in weightDict[seqIndex].keys():
#                        weightDict[seqIndex][parent] = float(1)
                        
                    if weightDict[childIndex][child] < np.multiply(weightDict[seqIndex][parent], nodeWeight):
                        weightDict[childIndex][child] = np.multiply(weightDict[seqIndex][parent], nodeWeight)
                        bestParent = parent
#                    else:
#                        weightDict[seqIndex][parent] = 
                pointBack[childIndex][child] = bestParent
        
            curState, seqIndex = self.getNextState(curState, seqIndex)
            print('next State = ' + curState + ' ' + str(seqIndex))
    
            print(weightDict)
#        print(pointBack)
        
        finalPath = self.backOut(weightDict, pointBack, 'source', 'sink')
#        print(finalPath)
        return weightDict, pointBack


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

    def getNextState(self, currentState, currentStateIndex):
        if currentStateIndex == 0:
            if currentState == 'S':
                return 'D1', currentStateIndex
            elif currentState.startswith('D'):
                if int(currentState[1:]) >= self.maxStateIndex:
                    return 'I0', currentStateIndex
                return ('D' + str(int(currentState[1:])+1)), currentStateIndex
            else:
                currentStateIndex+=1
                return 'M1', currentStateIndex
            
        idx = int(currentState[1:])
        if currentState.startswith('M'):
            return ('D' + str(idx)), currentStateIndex

        if currentState.startswith('D'):
            return ('I' + str(idx)), currentStateIndex
        
        if currentState.startswith('I'):
            currentStateIndex+=1
            idx+=1
            if(idx >= self.maxStateIndex):
                return 'E', currentStateIndex
            return ('M' + str(idx)), currentStateIndex

    def getPossibleParents(self, currentState, currentStateIndex):
        if currentStateIndex == 0:
            if currentState == 'D1':
                return ['S']
            else:

                return [('D' + str( int(currentState[1:])-1 ) )]
        if currentState == 'I0' or currentState == 'M1' or currentState == 'D1':
            return ['S', 'I0']
        if currentState == 'E':
            return ['M'+str(self.maxStateIndex), 'D'+str(self.maxStateIndex), 'I'+str(self.maxStateIndex)]
        
        idx = int(currentState[1:])
        if currentState.startswith('I'):
            return['M'+str(idx), 'D'+str(idx), 'I'+str(idx)]
        else:
            return['M'+str(idx-1), 'D'+str(idx-1), 'I'+str(idx-1)]
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
