#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian Chang, Lori Siao

"""
Viterbi learning uses initial random parameters and uses them and the emitted string to find the true parameters of the HMM. It gives the initial parameters to the viterbi algorithm, finds the most probable path for those parameters, then uses parameter estimation to create a new set of parameters. It iterates over these steps a set number of times.

Input: Number of iterations, an observed sequence, the allowed states, a matrix initial transition probabilities between states, and a matrix of initial emission probabilities of a symbol given each state.
Output: Optimized parameters of the HMM 

Dependencies:
DAGp15.py
prPip16.py
prOutcomep17.py
p23.py

"""
import pdb
import sys
import math
import numpy as np
from DAGp15 import DirectedAcyclicGraph
from prPip16 import ProbabilityHiddenPath
from prOutcomep17 import ProbabilitySeqGivenHiddenPath
from problem23 import BuildParameters

probHiddenPath = ProbabilityHiddenPath()
probSeq = ProbabilitySeqGivenHiddenPath()
def parseInput(toParse):
    """Takes the input file and parses out the number of interations,  the sequence and states into lists and the transition and emission matrices into dictionaries"""
    splitified = toParse.split('--------')
    iterations = splitified[0].rstrip().strip()
    sequence = list(splitified[1].rstrip().strip())
    alphabet = splitified[2].rstrip().strip().split()
    availableStates = splitified[3].rstrip().strip().split()
    transMatrix = splitified[4].rstrip().strip().splitlines()
    emissionMatrix = splitified[5].rstrip().strip().splitlines()

    transitionDict = probHiddenPath.getTransitionDict(transMatrix)
    emissionDict = probSeq.getEmissionDict(emissionMatrix)

    return(iterations, list(sequence), alphabet, transitionDict, emissionDict,availableStates)

def reParseMatrix(colHeader, rowHeader, matrix, isTrans):
    """Reformats the matrix into a graph that can be processed by the viterbi algorithm """
    unParsedList = [' '.join(colHeader)]

    for letter, row in zip(rowHeader, matrix):
        line = str(row).strip('[').rstrip(']')
        unParsedList.append(str(letter) + line)
    if isTrans == True:
        parsedDict = probHiddenPath.getTransitionDict(unParsedList)

    else:
        parsedDict = probSeq.getEmissionDict(unParsedList)

    return parsedDict

class ViterbiLearning():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""

    def makeGraph(self, sequence, availableStates, transitionDict):
        """Constructs a graph as a dictionary of transition edges with their edge weights"""
        dagEdges = {'source':{},'sink':{}}
        if sequence[0] != 'source':
            sequence.insert(0, 'source')
            sequence.append('sink')
            
        for state in availableStates:
            for inState in availableStates:
                if state in dagEdges.keys() and inState not in dagEdges[state]:
                    dagEdges[state][inState] = transitionDict[state + inState]
                else:
                    dagEdges[state] = {inState: transitionDict[state + inState]}
            dagEdges['source'][state] = 1-1/len(availableStates)
            dagEdges['sink'][state] = 1

        return dagEdges
    def runViterbiLearning(self, iterations, sequence, alphabet, availableStates, transitionDict, emissionDict):
        """Runs Viterbi Learning for a specified number of iterations using the viterbi algorithm and parameter estimation. Returns transition and emission parameters """
        dagProb = DirectedAcyclicGraph()
        buildParam = BuildParameters()
        permSequence = sequence.copy()
        for i in range(int(iterations)):
            #build viterbi graph from the parameters
            dagGraph = self.makeGraph(sequence, availableStates, transitionDict)
            #run the viterbi to get the most probable hidden path
            viterbiPath = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
            #get new transition parameters from parameter estimation
            transitionMatrix = buildParam.makeTransMatrix(availableStates, viterbiPath)
            #turn paramters into graph for the viterbi
            transitionDict = reParseMatrix(availableStates, availableStates, transitionMatrix, isTrans = True)
            #get new emission parameters from parameter estimation
            emissionMatrix = buildParam.makeEmissionMatrix(alphabet, availableStates, permSequence, viterbiPath)
            #turn parameters into graph for the viterbi
            emissionDict = reParseMatrix(alphabet, availableStates, emissionMatrix, isTrans = False)

        return transitionMatrix, emissionMatrix
        
def main():

    toParse = sys.stdin.read()

    viterbiStructures = ViterbiLearning()
    dagProb = DirectedAcyclicGraph()
    buildParam = BuildParameters()
    iterations, sequence, alphabet, transitionDict, emissionDict, availableStates = parseInput(toParse)
    transitionMatrix, emissionMatrix = viterbiStructures.runViterbiLearning(iterations, sequence, alphabet, availableStates, transitionDict, emissionDict)
    
    buildParam.printMatrix(availableStates, availableStates, transitionMatrix)
    print('--------')
    buildParam.printMatrix(availableStates, alphabet, emissionMatrix)

if __name__== "__main__":
    main()

