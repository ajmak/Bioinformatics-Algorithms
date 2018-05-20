#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: 

"""
Uses the forward/backward algorithm of HMMs to calculate the probability that each node is in the most probable hidden path given the transition probabilities between hidden states and the emission probabilities of each sequence given a hidden state                                                  

Input: An observed sequence, alphabet, the allowed states, a matrix of transition probabilities between states, and a matrix of emission probabilitiesof an observed sequence given each state.                            
Output: The probability at each state that that state belongs in the most probabl hidden path given the emission based on the HMM profile  

Dependencies:
DAGp19.py
HMMclasses

"""
import pdb
import sys
import math
import numpy as np
from HMMclasses import HMM
from DAGp19 import DirectedAcyclicGraph
from problem25 import ForwardBackward


def parseInput(toParse, HMMstructs):
    """Takes the input file and parses it the sequence and states into lists and the transition and emission matrices into dictionaries"""
    splitified = toParse.split('--------')
    iterations = splitified[0].rstrip().strip()
    sequence = list(splitified[1].rstrip().strip())
    alphabet = splitified[2].rstrip().strip().split()
    availableStates = splitified[3].rstrip().strip().split()
    transMatrix = splitified[4].rstrip().strip().splitlines()
    emissionMatrix = splitified[5].rstrip().strip().splitlines()
    
    transitionDict = HMMstructs.getTransitionDict(transMatrix)
    emissionDict = HMMstructs.getEmissionDict(emissionMatrix)

    return(iterations, sequence, alphabet, transitionDict, emissionDict,availableStates)

class BaulmWelch():

    def getPiStar(self, emissionDict, sequence, alphabet, availableStates):
        emisMatrix = np.zeros(shape = (len(availableStates), len(sequence)))
        stateCols = {state: {ltr: float(0) for ltr in alphabet} for state in availableStates}

        for row in range(len(availableStates)):
            for col in range(len(sequence)):
                probEmis = np.float64(emissionDict[str(availableStates[row]) + str(sequence[col])])

                emisMatrix[row, col] = probEmis

                print(str(availableStates[row]) + str(sequence[col])) 
                stateCols[availableStates[row]][sequence[col]] += probEmis
        print(emissionDict)
        print(sequence)
        print(stateCols)
        print(emisMatrix)
        return emisMatrix, stateCols

    def getPiStarStar(self, sequence, availableStates, forward, backward, forwardSum, transitionDict, emissionDict):
        """Builds Pi** matrix by multiplying:
(forward probability of a node, l) x (the backwards probability of the next node, k) x (transition(l,k) x emission(k)"""
        piStSt = np.zeros(shape = (len(transitionDict.keys()), len(sequence)-1))
        rowState = {state:[] for state in availableStates}
        transitions = list(transitionDict.keys())
        for row in range(len(transitions)):
            l = transitions[row][:-1]
            k = transitions[row][1:]
            rowState[l].append(row)
            for col in range(len(sequence)-1):
                #fw/bkwd dicts are 1 based, need col+1
                i = col + 1
                FWxBW = np.multiply(np.float64(forward[i][l]), np.float64(backward[i + 1][k]))
                translkxEmitk = np.multiply(np.float64(transitionDict[transitions[row]]), np.float64(emissionDict[str(k) + sequence[col+1]]))
                prodAll = np.multiply(FWxBW, translkxEmitk)
                normalAll = np.divide(prodAll, forwardSum)
                piStSt[row, col] = normalAll
        return piStSt, rowState, transitions

    def emisMatrixtoDict(self, availableStates, alphabet, piSt, emisCounts):
        """Makes a normalized dictionary of emission probabilities for each state given an emission matrix """
        emissionDict = {}
        rowSums = piSt.sum(axis=0)
        for i in range(len(availableStates)):
            for j in range(len(alphabet)):
                emisProb = np.divide(emisCounts[availableStates[i]][alphabet[j]], rowSums[i])
                emissionDict[str(availableStates[i]) + str(alphabet[j])] = emisProb
        print(emissionDict)
        print('emissiondict')
        return emissionDict
    def transMatrixtoDict(self, availableStates, piStSt, rowState, transitions):
        """Makes a normalized dictionary of transition probabilities based on a given transition matrix """
        transitionDict = {}
        normalSums = self.normalize(piStSt, rowState)
        for i in range(len(transitions)):
            transitionDict[transitions[i]] = normalSums[i]
#        print(transitionDict)
        return transitionDict
        
    def normalize(self, matrix, indexDict):
        """Normalizes a transition by sum of transitions from a given state given the matrix and a dictionary of indices outlining which rows begin at the same transition."""
        rowSums = matrix.sum(axis=1)

        for state in indexDict.keys():
            total = 0
            for row in indexDict[state]:
                total += rowSums[row]
            for row in indexDict[state]:
                rowSums[row] = np.divide(rowSums[row], total)
#        print(rowSums)
        return rowSums
    
def runBaulmWelch(iterations, sequence, alphabet, transitionDict, emissionDict, availableStates):
    """Runs Baulm Welch parameter learning algorithm for given iterations and returns the final transition/emission parameters for the emitted sequence"""

    dagProb = DirectedAcyclicGraph()
    baulmWelch = BaulmWelch()
    forwardBackward = ForwardBackward()

    for i in range(int(iterations)):
        #run forward/backward algorithm on parameters
        accurateProbs, forward, forwardSum, backward = forwardBackward.getForwardBackward(sequence, transitionDict, emissionDict, availableStates, dagProb)

        #build pi* and pi**
        piStar, emisSums = baulmWelch.getPiStar(emissionDict, sequence, alphabet, availableStates)
        piStarStar, rowState, transitions = baulmWelch.getPiStarStar(sequence, availableStates, forward, backward, forwardSum, transitionDict, emissionDict)
        #rebuild normalized emission and transition dictionaries
        emissionDict = baulmWelch.emisMatrixtoDict(availableStates, alphabet, piStar, emisSums)
        transitionDict = baulmWelch.transMatrixtoDict(availableStates, piStarStar, rowState, transitions)

    return emissionDict, transitionDict

def printDictionary(rowHeaders, colHeaders, dictionary):
    """Takes either a transition or emission dictionary and prints as matrix based on row and col headers """
    printMatrix = np.zeros(shape = (len(rowHeaders), len(colHeaders))) 
    for row in range(len(rowHeaders)):
        for col in range(len(colHeaders)):
            printMatrix[row, col] = dictionary[str(rowHeaders[row]) + str(colHeaders[col])]
    print('\t' + '\t'.join(colHeaders))
    for header, row in zip(rowHeaders, printMatrix):
        print('%s \t %s' % (header, '\t'.join(('%.3f' % i).rstrip('0').rstrip('.') if i!=1.0 else ('%.1f' % i) for i in row)))
    
def main():

    toParse = sys.stdin.read()

    HMMstructs = HMM()

    iterations, sequence, alphabet, transitionDict, emissionDict, availableStates = parseInput(toParse, HMMstructs)

    finalEmission, finalTransition = runBaulmWelch(iterations, sequence, alphabet, transitionDict, emissionDict, availableStates)

    printDictionary(availableStates, availableStates, finalTransition)
    print('--------')
    printDictionary(availableStates, alphabet, finalEmission)

if __name__== "__main__":
    main()

