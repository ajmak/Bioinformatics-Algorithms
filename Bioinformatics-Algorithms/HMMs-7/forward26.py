#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian Chang, Lori Siao

"""
Uses the forward algorithm to calculate the probability that the observed sequence is emitted by the HMM given the transition probabilities between hidden states and the emission probabilities of each sequence given a hidden state                                                  

Input: An observed sequence, the allowed states, a matrix of transition probabilities between states, and a matrix of emission probabilitiesof an observed sequence given each state.                            
Output: The probability that the observed sequence is emitted based on the HMM profile  

Dependencies:
DAGp19.py
prPip16.py
prOutcomep17.py

"""
import pdb
import sys
import math
from DAGp19 import DirectedAcyclicGraph
from prPip16 import ProbabilityHiddenPath
from prOutcomep17 import ProbabilitySeqGivenHiddenPath

class BaulmWelch():
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
        
        transitionDict = ProbabilityHiddenPath.getTransitionDict(self, transMatrix)
        emissionDict = ProbabilitySeqGivenHiddenPath.getEmissionDict(self, emissionMatrix)

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
                    dagEdges[state] = {inState: transitionDict[state + inState]}
            dagEdges['source'][state] = (1-1/len(availableStates))
            dagEdges['sink'][state] = 1
        return dagEdges

    def reverseGraph(self, availableStates, forwardGraph):
        reversedGraph = {}
        for k1, subdict in forwardGraph.items():
            for k2, value in subdict.items():
                reversedGraph.setdefault(k2, {})[k1] = value
        reversedGraph['sink'] = {}
        for state in availableStates:
            reversedGraph['sink'][state] = 1
            
        return reversedGraph

    def calcEachNode(self, forwardSum, forwardWeights, backwardWeights, sequence, availableStates):
        printList = [[] for i in range(len(sequence)-1)]
        for i in range(1, len(sequence)-1):
            for state in availableStates:
                prob = ((forwardWeights[i][state]*backwardWeights[i+1][state])/forwardSum)
                printList[i].append(('%.4f' % prob).rstrip('0').rstrip('.'))
        printList[0] = '\t'.join(availableStates)        
        return printList

def main():

    toParse = sys.stdin.read()

    HMMstructures = BaulmWelch()
    sequence, transitionDict, emissionDict, availableStates = HMMstructures.parseInput(toParse)

    forwardGraph = HMMstructures.makeGraph(sequence, availableStates, transitionDict)
    reversedGraph = HMMstructures.reverseGraph(availableStates, forwardGraph)

    dagProb = DirectedAcyclicGraph()
    forwardWeights, forwardSum = dagProb.longestPath(1, len(sequence), 1, sequence, forwardGraph,  emissionDict)

    backWeights, backSum = dagProb.longestPath(len(sequence)-1, 0, -1, sequence, reversedGraph, emissionDict)
    stateProbabilities = HMMstructures.calcEachNode(forwardSum, forwardWeights, backWeights, sequence, availableStates)

    for prob in stateProbabilities:
        print(*prob, sep="\t")

if __name__== "__main__":
    main()

