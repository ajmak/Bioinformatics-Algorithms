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
from HMMclasses import HMM

class Viterbi():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""
    
    def makeGraph(self, sequence, availableStates, transitionArray):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        dagEdges = {'source':{}}
        #beginning block of permissible states
        

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

def parseInput(toParse):
    """Parses the input file into the threshold, pseudocount(sigma), alphabet, and a matrix of the provided alignment"""
    splitified = toParse.split('--------')
    sequence = splitified[0].rstrip().strip()
    threshold, sigma = splitified[1].rstrip().strip().split()
    alphabet = splitified[2].rstrip().strip().split()
    alignment = splitified[3].rstrip().strip().splitlines()
    alignList = [list(i) for i in alignment]

    alignArray = np.array(alignList)

    return sequence, threshold, sigma, alphabet, alignArray

def main():

    toParse = sys.stdin.read()

    viterbiStructures = Viterbi()
    HMMs = HMM()

    sequence, threshold, sigma, alphabet, alignArray = parseInput(toParse)

    lenStates, stateList, transitionArray, emissionArray = HMMs.buildTransitionEmissionwithPseudo(alphabet, threshold, sigma, alignArray)
    print(transitionArray)

    edgeDict = viterbiStructures.makeGraph(transitionMatrix)
    
    dagProb = DirectedAcyclicGraph()
    
    getViterbi = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
    

if __name__== "__main__":
    main()

