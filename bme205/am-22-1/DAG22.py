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
        pointBack = {key: {state: 'parent' for state in edgeDict.keys()}for key in range(len(stateList)-1)}
        weightDict = {stateIdx: {} for stateIdx in range(len(sequence)+2)}
        weightDict[0] = {'S': float(1), 'I0': float(1)}
        curState = 'S'
        seqIndex = 0
        
        print("!!!!EDGE_DICT!!!!!")
        print(edgeDict)
        while curState != 'STOP':
            for child in edgeDict[curState].keys():
                # finds max weighted path to that child node
                childIndex = seqIndex
                if not child.startswith('D'): childIndex += 1
                #if child[:-1] != 'D': childIndex += 1

                #can only end in 'E'
                if childIndex == self.maxStateIndex + 1:
                    if child != 'E':
                        continue;
                        
                if child == 'E':
                    if seqIndex != self.maxStateIndex:
                        continue
                
                for parent in self.getPossibleParents(child, childIndex):
                    #'M0' and 'I0' do not exist. can't be parents
                    if seqIndex == 0:
                        if parent.startswith('M') or parent.startswith('I'):
                            continue
                    #'S' does not exist at columns greater than 0
                    elif parent == 'S':
                        continue
                    #if curState.startswith('D'): parentIndex += 1
                    print('******************')
                    print(child + ' child' + ' index: ' + str(childIndex))
                    print(parent + ' parent' + ' index: ' + str(seqIndex))
                    print('possible parents:')
                    print(self.getPossibleParents(child, childIndex))
                    print(sequence)
                    print(len(sequence))
                    #print(sequence[childIndex-1])
                    print(str(seqIndex) + ' seqIndex')
                    print('weightDict:')
                    #print(weightDict)
                    print('emission')
                    #print(emissionDict)
                    if child == 'E':
                        nodeWeight = edgeDict[parent][child]
                    else:
                        nodeWeight = np.multiply(np.float64(edgeDict[parent][child]), np.float64(emissionDict[child + sequence[childIndex-1]]))
                    print(nodeWeight)
                    print('******************')
                    if child not in weightDict[childIndex].keys():
                        weightDict[childIndex][child] = float('-inf')
                    #if parent not in weightDict[seqIndex].keys():
                        #weightDict[seqIndex][parent] = float(1)
                        
                    if weightDict[childIndex][child] < np.multiply(weightDict[seqIndex][parent], nodeWeight):
                        weightDict[childIndex][child] = np.multiply(weightDict[seqIndex][parent], nodeWeight)
                        bestParent = parent
#                    else:
#                        weightDict[seqIndex][parent] = 
                pointBack[childIndex][child] = bestParent
        
            print('getNextState(' + curState + ' ' + str(seqIndex) + ')')
            curState, seqIndex = self.getNextState(curState, seqIndex)
            print('curState: ' + curState)
            print('seqIndex: ' + str(seqIndex))
            print('maxStateIndex: ' + str(self.maxStateIndex))
            #print('weightDict:')
            #print(weightDict)
#        print(pointBack)
        
        finalPath = self.backOut(weightDict, pointBack, 'S', 'E')
#        print(finalPath)
        return weightDict, pointBack


    def backOut(self, weightDict, pointBack, source, sink):
        """Creates a list of the longest path from the source to the sink given a dictionary of back pointers that outline the route taken to each node and a dictionary of cumulative weights up until each node """
        curNode = sink
        outList = []
        print(weightDict[0])
        print(pointBack[0])
        
        i = len(weightDict.keys())-1
        while curNode != 'S' and i >= 0:
            if curNode == 'E':
                print('EEEEEEEEEE')
                curNode = pointBack[i][curNode]
                outList.append(curNode + ' ')
            else:
#                outList.append(curNode + ' ')
                curNode = pointBack[i][curNode]
                print('append ' + str(i) + str(curNode))
                outList.append(curNode + ' ')
            print(str(curNode) + ' curNode')
            print(i)
            prevNode = curNode

            i -= 1
            if len(outList) > 2:
                if outList[-2].startswith('D'):
                    i += 1
                
        outList = reversed(outList)
        print('')
        print('outList: ' + ''.join(outList))
        #print(pointBack[1]['D2'])

    def getNextState(self, currentState, currentStateIndex):
        if currentStateIndex == 0:       
            if currentState == 'S':
                return 'D1', currentStateIndex
            elif currentState.startswith('D'):
                if int(currentState[1:]) >= self.maxStateIndex+1:
                    currentStateIndex+=1
                    return 'I0', currentStateIndex
                return ('D' + str(int(currentState[1:])+1)), currentStateIndex
            else:
                return 'M1', currentStateIndex
                
        if currentState.startswith('E'):
            return 'STOP', 69
            
        idx = int(currentState[1:])
        
        if currentState.startswith('M'):
            return ('D' + str(idx)), currentStateIndex
            
        if currentState.startswith('D'):
            return ('I' + str(idx)), currentStateIndex
            
        if currentState.startswith('I'):
            idx+=1
            if(idx > self.maxStateIndex+1):
                currentStateIndex+=1
                if currentStateIndex > self.maxStateIndex:
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
            return ['M'+str(currentStateIndex), 'D'+str(currentStateIndex), 'I'+str(currentStateIndex)]
            
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
    sys.stdin = open('/dev/tty')
    print(finalWeight)
    finalPath.reverse()
    print('->'.join(finalPath))
    

if __name__== "__main__":
    main()
