#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian Chang, Lori Siao

"""
The viterbi algorithm finds the most probable path of hidden states that corresponds to the given observed sequence. It uses the longest path algorithm for a directed acyclic graph to find the path of the hidden states.

Input: An observed sequence, the allowed states, a matrix of transition probabilities between states, and a matrix of emission probabilities of an observed sequence given each state.
Output: The longest (most probable) path of hidden states for the provided sequence

Dependencies:
DAGp15.py
prPip16.py
prOutcomep17.py

"""
import pdb
import sys
import math
from DAGp15 import DirectedAcyclicGraph
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
        emissionDict = ProbabilitySeqGivenHiddenPath.getEmissionDict(self, emissionMatrix)
        transitionDict.update((x, math.log(float(y))) if float(y) != 0.0 else float('-inf') for x, y in transitionDict.items())
        emissionDict.update((x, math.log(float(y))) if float(y) != 0.0 else float('-inf') for x, y in emissionDict.items())

        return(list(sequence), transitionDict, emissionDict,availableStates)

    def makeGraph(self, sequence, availableStates, transitionDict):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        dagEdges = {'source':{},'sink':{}}

        for state in availableStates:
            for inState in availableStates:
                if state in dagEdges.keys() and inState not in dagEdges[state]:
                    dagEdges[state][inState] = transitionDict[state + inState]

                else:
                    dagEdges[state] = {inState: transitionDict[state + inState]}
            dagEdges['source'][state] = math.log(1-1/len(availableStates))
            dagEdges['sink'][state] = math.log(1)
        print(dagEdges)
        return dagEdges



def main():

    toParse = sys.stdin.read()

    HMMstructures = HMM()
    sequence, transitionDict, emissionDict, availableStates = HMMstructures.parseInput(toParse)

    dagGraph = HMMstructures.makeGraph(sequence, availableStates, transitionDict)
    dagProb = DirectedAcyclicGraph()
    
    getViterbi = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
    

if __name__== "__main__":
    main()

