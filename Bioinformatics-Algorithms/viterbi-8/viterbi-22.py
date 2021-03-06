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
import numpy as np
from DAG22 import DirectedAcyclicGraph
from HMMclasses import HMM

class Viterbi():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""
    
    def makeGraph(self, lenStates, stateList, transitionArray):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        dagEdges = {state: {} for state in stateList}

        
        #inserts permissible states from start into graph
        coli = 1
        coli1 = 4
        rowi = 0
        rowi1 = 2
        for i in range(coli, coli1):
            for j in range(rowi, rowi1):
                dagEdges[stateList[j]][stateList[i]] = transitionArray[j,i]
        #inserts mid-range of permissible states into graph
        col0 = coli1
        col1 = col0 + 3
        row0 = rowi1
        row1 = row0 + 3
        realRow = rowi1
        realCol = 0
        while col1 < lenStates:
            for i in range(col0, col1):
                for j in range(row0, row1):
                    dagEdges[stateList[j]][stateList[i]] = transitionArray[j, i]
            col0 += 3
            col1 += 3
            row0 += 3
            row1 += 3
            realRow +=1
        #inserts end permissible states into graph
        col0 = lenStates - 2
        col1 = lenStates
        row0 = lenStates - 4
        row1 = lenStates - 1
        for i in range(col0, col1):
            for j in range(row0, row1):
                dagEdges[stateList[j]][stateList[i]] = transitionArray[j, i]

        return dagEdges

    def getNextNode(currentState, currentStateIndex, maxStateIndex):
        if currentStateIndex == 0:
            if currentState == 'S':
                return 'D1'
            elif currentState.startswith('D'):
                if int(currentState[1:]) >= maxStateIndex:
                    return 'I0'
                return 'D' + str(int(currentState[1:])+1)
            else:
                return 'M1'
            
        idx = int(currentState[1:])
        if currentState.startswith('M'):
            return 'D' + str(idx)
        
        if currentState.startswith('D'):
            return 'I' + str(idx)
        
        if currentState.startswith('I'):
            idx+=1
            if(idx > maxStateIndex):
                return 'E'
            return 'M' + str(idx)
                        
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

    dagGraph = viterbiStructures.makeGraph(lenStates, stateList, transitionArray)

    emissionDict = HMMs.buildEmissionDict((emissionArray).tolist(), alphabet, stateList)

    dagProb = DirectedAcyclicGraph(lenStates, sequence)
    
    getViterbi = dagProb.longestPath(sequence, stateList, dagGraph, emissionDict)
    

if __name__== "__main__":
    main()

