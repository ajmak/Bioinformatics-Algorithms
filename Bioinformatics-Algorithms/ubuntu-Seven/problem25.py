#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: 

"""
Uses the forward/backward algorithm of HMMs to calculate the probability that each node is in the most probable hidden path given the transition probabilities between hidden states and the emission probabilities of each sequence given a hidden state                                                  

Input: An observed sequence, alphabet, the allowed states, a matrix of transition probabilities between states, and a matrix of emission probabilitiesof an observed sequence given each state.                            
Output: The probability at each state that that state belongs in the most probabl hidden path given the emission based on the HMM profile  

Dependencies:
DAGp19.py
prPip16.py
prOutcomep17.py

"""
import pdb
import sys
import math
import numpy as np
from DAGp19 import DirectedAcyclicGraph
from HMMclasses import HMM
HMMstructs = HMM()
class ForwardBackward():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""
    def parseInput(self, toParse):
        """Takes the input file and parses it the sequence and states into lists and the transition and emission matrices into dictionaries"""
        splitified = toParse.split('--------')
        sequence = list(splitified[0].rstrip().strip())
        sequence.insert(0, 'source')
        sequence.append('sink')
        availableStates = splitified[2].rstrip().strip().split()
        transMatrix = splitified[3].rstrip().strip().splitlines()
        emissionMatrix = splitified[4].rstrip().strip().splitlines()
        
        transitionDict = HMMstructs.getTransitionDict(transMatrix)
        emissionDict = HMMstructs.getEmissionDict(emissionMatrix)

        return(list(sequence), transitionDict, emissionDict,availableStates)

    def makeGraph(self, sequence, availableStates, transitionDict):
        """Constructs a graph as a dictionary of transition edges with their edge weights"""
        dagEdges = {'source':{},'sink':{}}
        reversedEdges = {'source':{}, 'sink':{}}
        for state in availableStates:
            for inState in availableStates:
                if state in dagEdges.keys() and inState not in dagEdges[state]:
                    dagEdges[state][inState] = transitionDict[state + inState]
                else:
                    dagEdges[state] = {inState: np.float64(transitionDict[state + inState])}
            dagEdges['source'][state] = np.float64(1-1/len(availableStates))
            dagEdges['sink'][state] = 1
        return dagEdges

    def reverseGraph(self, availableStates, forwardGraph):
        """Builds a reversed graph of transition edges so that each node can be read backwards without altering its transitional probability """
        reversedGraph = {}
        for k1, subdict in forwardGraph.items():
            for k2, value in subdict.items():
                reversedGraph.setdefault(k2, {})[k1] = value
        reversedGraph['sink'] = {}
        for state in availableStates:
            reversedGraph['sink'][state] = 1
            
        return reversedGraph

    def calcEachNode(self, forwardSum, forwardWeights, backwardWeights, sequence, availableStates):
        """Takes a dictionary of indexed forward probabilities and backward probabilities for each node and uses them to calculate the probability that each node is in the most probable hidden path """
        backwardWeights = {key-1: v for key, v in backwardWeights.items()}
        print(forwardWeights)
        print()
        print(backwardWeights)
        print(len(sequence))
        accurateList = [[] for i in range(len(sequence)-1)]
        printList = [[] for i in range(len(sequence)-1)]
        for i in range(1, len(sequence)-1):
            for state in availableStates:
                preDiv = np.multiply((forwardWeights[i][state]), (backwardWeights[i][state]))
                prob = np.divide(preDiv, forwardSum)
                accurateList[i].append(prob)
                printList[i].append(('%.4f' % prob).rstrip('0').rstrip('.'))
        accurateList[0] = [state for state in availableStates]
        printList[0] = [state for state in availableStates]
        return accurateList, printList

def main():

    toParse = sys.stdin.read()

    HMMstructures = ForwardBackward()
    sequence, transitionDict, emissionDict, availableStates = HMMstructures.parseInput(toParse)

    forwardGraph = HMMstructures.makeGraph(sequence, availableStates, transitionDict)
    reversedGraph = HMMstructures.reverseGraph(availableStates, forwardGraph)

    dagProb = DirectedAcyclicGraph()
    forwardWeights, forwardSum = dagProb.longestPath(1, len(sequence), 1, sequence, forwardGraph,  emissionDict)

    backWeights, backSum = dagProb.longestPath(len(sequence)-1, 0, -1, sequence, reversedGraph, emissionDict)
    accurateProbs, printProbs = HMMstructures.calcEachNode(forwardSum, forwardWeights, backWeights, sequence, availableStates)

    for prob in printProbs:
        print(*prob, sep="\t")

if __name__== "__main__":
    main()

