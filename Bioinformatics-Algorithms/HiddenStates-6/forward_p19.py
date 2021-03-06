#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Yianni Anastopulos

"""
Uses the forward algorithm to calculate the probability that the observed sequence is emitted by the HMM given the transition probabilities between hidden states and the emission probabilities of each sequence given a hidden state
Input: An observed sequence, the allowed states, a matrix of transition probabilities between states, and a matrix of emission probabilities of an observed sequence given each state.
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

class HMM():
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
        probSeq = ProbabilitySeqGivenHiddenPath

        emissionDict = probSeq.getEmissionDict(self, emissionMatrix)

        return(list(sequence), transitionDict, emissionDict, availableStates)

    def makeGraph(self, sequence, availableStates, transitionDict):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        dagEdges = {'source':{},'sink':{}}
        for state in availableStates:
            for inState in availableStates:
                if state in dagEdges.keys() and inState not in dagEdges[state]:
                    dagEdges[state][inState] = transitionDict[state + inState]

                else:
                    dagEdges[state] = {inState: transitionDict[state + inState]}
            dagEdges['source'][state] = (1-1/len(availableStates))
            dagEdges['sink'][state] = 1
        return dagEdges


class ProbabilitySeqGivenHiddenPath:
    def __init__(self):
        self.totalProbability = 1
         
    def calcProbStatePath(self, sequence, statesPath, transitionDict, availableStates):
        seqList = list(sequence)
        statesList = list(statesPath)
        stateCount = 0
        for i in range(len(seqList)):
            seqState = str(statesList[i]) + str(seqList[i])
            self.totalProbability *= float(transitionDict[seqState])
            prevState = statesList[i]
        
        return self.totalProbability

def main():
    toParse = sys.stdin.read()
    HMMstructures = HMM()

    sequence, transitionDict, emissionDict, availableStates = HMMstructures.parseInput(toParse)

    dagGraph = HMMstructures.makeGraph(sequence, availableStates, transitionDict)
    dagProb = DirectedAcyclicGraph()
    
    getViterbi = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
    

if __name__== "__main__":
    main()

