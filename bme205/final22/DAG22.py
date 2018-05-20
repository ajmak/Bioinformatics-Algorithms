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
        self.maxStateIndex = int((lenStates - 3) / 3)
        self.lenSeq = len(sequence)

    def longestPath(self, sequence, stateList, edgeDict, emissionDict):
        """Finds the longest path in a topologically ordered graph given the source of the subgraph. Returns a dictionary of cumulative weights for a path up to each node and a dictionary of back pointers """
        pointBack = {key: {state: 'parent' for state in edgeDict.keys()}for key in range(self.lenSeq+2)}
        weightDict = {stateIdx: {} for stateIdx in range(len(sequence)+2)}
        weightDict[0] = {'S': float(1), 'I0': float(1)}
        curState = 'S'
        seqIndex = 0
 
        while curState != 'STOP':
            for child in edgeDict[curState].keys():
                # finds max weighted path to that child node
                childIndex = seqIndex
                if not child.startswith('D'): childIndex += 1

                #can only end in 'E'
                if childIndex == self.lenSeq + 1:
                    if child != 'E':
                        continue
                
                #can only take an edge to E at index == lenSeq
                if child == 'E':
                    if seqIndex != self.lenSeq:
                        continue
                
                for parent in self.getPossibleParents(child, childIndex):
                    #'M0' and 'I0' do not exist. can't be parents
                    if seqIndex == 0:
                        if parent.startswith('M') or parent.startswith('I'):
                            continue
                    #'S' does not exist at columns greater than 0
                    elif parent == 'S':
                        continue

                    if child == 'E':
                        nodeWeight = edgeDict[parent][child]
                    else:
                        nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[child + sequence[childIndex-1]]))

                    if child not in weightDict[childIndex].keys():
                        weightDict[childIndex][child] = float('-inf')
                        
                    if weightDict[childIndex][child] <= np.multiply(weightDict[seqIndex][parent], nodeWeight):
                        weightDict[childIndex][child] = np.multiply(weightDict[seqIndex][parent], nodeWeight)
                        bestParent = parent

                pointBack[childIndex][child] = bestParent
                
            curState, seqIndex = self.getNextState(curState, seqIndex)
            
        finalPath = self.backOut(weightDict, pointBack, 'S', 'E')

        return weightDict, pointBack


    def backOut(self, weightDict, pointBack, source, sink):
        """Creates a list of the longest path from the source to the sink given a dictionary of back pointers that outline the route taken to each node and a dictionary of cumulative weights up until each node """
        
        outList = []       
        i = len(weightDict.keys())-1
        curNode = pointBack[i][sink]
        
        while curNode != 'S':
            outList.append(curNode + ' ')
            curNode = pointBack[i-1][curNode]
            i -= 1
            if len(outList) > 2:
                if outList[-2].startswith('D'):
                    i += 1
        outList = reversed(outList)
        
        print('')
        print(''.join(outList))


    def getNextState(self, currentState, currentStateIndex):
        if currentStateIndex == 0:       
            if currentState == 'S':
                return 'D1', currentStateIndex
            elif currentState.startswith('D'):
                if int(currentState[1:]) >= self.maxStateIndex:
                    currentStateIndex+=1
                    return 'I0', currentStateIndex
                return ('D' + str(int(currentState[1:])+1)), currentStateIndex
            else:
                return 'M1', currentStateIndex
                
        if currentState.startswith('E'):
            return 'STOP', -1
            
        idx = int(currentState[1:])
        
        if currentState.startswith('M'):
            return ('D' + str(idx)), currentStateIndex
            
        if currentState.startswith('D'):
            return ('I' + str(idx)), currentStateIndex
            
        if currentState.startswith('I'):
            idx+=1
            if(idx > self.maxStateIndex):
                currentStateIndex+=1
                if currentStateIndex > self.lenSeq:
                    return 'E', currentStateIndex
                return 'I0', currentStateIndex
            return ('M' + str(idx)), currentStateIndex


    def getPossibleParents(self, currentState, currentStateIndex):
        if currentStateIndex == 0:
            if int(currentState[1:]) == 1:
                return ['S']
            else:
                return [('D' + str(int(currentState[1:])-1))]
                
        if currentState == 'I0':
            if currentStateIndex > 1:
                return ['I0']
            return ['S', 'I0']
            
        if currentState == 'M1' or currentState == 'D1':
            return ['S', 'I0']
            
        if currentState == 'E':
            return ['M'+str(self.maxStateIndex), 'D'+str(self.maxStateIndex), 'I'+str(self.maxStateIndex)]
            
        idx = int(currentState[1:])
        if currentState.startswith('I'):
            return ['M'+str(idx), 'D'+str(idx), 'I'+str(idx)]
        else:
            return ['M'+str(idx-1), 'D'+str(idx-1), 'I'+str(idx-1)]
            
def main():

    nodes = sys.stdin.read().splitlines()
    source = nodes.pop(0)
    sink = nodes.pop(0)
    graphKmers = DirectedAcyclicGraph(nodes)

    edgeDict, backDict, candidates = graphKmers.getEdges(nodes)

    topologicalOrder = graphKmers.topologicalSort(edgeDict, backDict, candidates)

    weightDict, pointBack = graphKmers.longestPath(topologicalOrder, edgeDict, source)

    finalWeight, finalPath = graphKmers.backOut(weightDict, pointBack, source, sink)

    finalPath.reverse()
    

if __name__== "__main__":
    main()