#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Calculates the probability of a hidden path of states 
Input: A path of states and a matrix of transition probabilities between all states
Output: The probability of the given path based on the provided transition probabilities 
"""

import sys

import math
import numpy as np
from problem20 import HMMProfile

class HMM:
    """Takes the hidden state path and the dictionary of transition probabilities between states and calculates the probability of the hidden path"""
    def __init__(self):
        self.totalProbability = 1

    def getTransitionDict(self, probMatrix):
        """Builds a dictionary of transition probabilities based on the transition matrix provided"""
        line = 0
        transitionDict = {}
        for prob in probMatrix:
            prob = prob.split()
            if line == 0:
                orderedStates = prob
            else:
                state1 = prob.pop(0)
                for i in range(len(orderedStates)):
                    transitionDict[str(state1) + str(orderedStates[i])] = prob[i]
            line += 1
        return(transitionDict)
         
    def calcProbTransition(self, statesPath, transitionDict, availableStates):
        """Uses the knonw hidden path of states and the dictionary of transition probabilities to calculate the probability of the path"""

        stateCount = 0
        for state in statesPath:
            if stateCount == 0:
                prevState = state
                stateCount += 1
                continue

            transition = str(prevState) + str(state)
            self.totalProbability *= float(transitionDict[transition])
            prevState = state
        beginningProb = 1/availableStates
        self.totalProbability *= beginningProb
        return self.totalProbability

    def getEmissionDict(self, probMatrix):
        """Builds a dictionary of probabilities that a given sequence occurs based on each available state"""
        line = 0
        emissionDict = {}
        for prob in probMatrix:
            prob = prob.split()
            if line == 0:
                orderedStates = prob

            else:
                state1 = prob.pop(0)
                for i in range(0,len(orderedStates)):
                    emissionDict[str(state1) + str(orderedStates[i])] = prob[i]
            line += 1
        return(emissionDict)

    def calcProbEmission(self, sequence, statesPath, emissionDict):
        """Calculates the probability of an observed sequence given its hidden states"""
        seqList = list(sequence)
        statesList = list(statesPath)

        for i in range(len(seqList)):
            seqState = str(statesList[i]) + str(seqList[i])
            self.totalProbability *= float(emissionDict[seqState])
            prevState = statesList[i]
        
        return self.totalProbability


    def addPseudoTrans(self, sigma, lenStates, transitionArray):
        """Incorporates pseudocounts for allowed states in transition matrix and normalizes """
        sigma = float(sigma)
        sigArray = np.zeros(shape = (lenStates, lenStates), dtype = float)
        
        transitionArray[:, 0:lenStates - 2] = transitionArray[:, 0:lenStates - 2]*(1-3*sigma)

        transitionArray[:,lenStates - 2:] = transitionArray[:,lenStates - 2:]*(1-2*sigma)
        #beginning block of permissible states
        coli = 1
        coli1 = 4
        rowi = 0
        rowi1 = 2
        for i in range(coli, coli1):
            for j in range(rowi, rowi1):
                sigArray[j,i] += sigma
        #Mid permissible states
        col0 = coli1
        row0 = rowi1
        col1 = col0 + 3
        row1 = row0 + 3
        while col1 < lenStates:
            for i in range(col0, col1):
                for j in range(row0, row1):
                    sigArray[j,i] += sigma
            col0 += 3
            row0 += 3
            col1 += 3
            row1 += 3
        #End block of permissible states
        col0 = lenStates - 2
        col1 = lenStates 
        row0 = lenStates - 4
        row1 = lenStates - 1
        for i in range(col0, col1):
            for j in range(row0, row1):
                sigArray[j,i] += sigma
        addMatrix = np.add(transitionArray, sigArray)
        pseudoTrans = self.normalize(sigma, lenStates, addMatrix)
        return pseudoTrans

    def addPseudoEmis(self, sigma, lenStates, alphabet, emissionArray):
        sigma = float(sigma)
        sigArray = np.zeros(shape = (lenStates, len(alphabet)), dtype = float)
        emissionArray = emissionArray*(1-len(alphabet)*sigma)
        sigArray = sigArray + sigma
        sigArray[::3] = 0
        sigArray[lenStates-1] = 0
        addMatrix = np.add(emissionArray, sigArray)
        pseudoTrans = self.normalize(sigma, lenStates, addMatrix)

        return pseudoTrans
    
    def normalize(self, sigma, lenStates, array):
        """Normalizes items in array based on the sum of  """
        rowSums = array.sum(axis=1)
        rowSums = np.vectorize(lambda x: 1.0 if x == 0 else x)(rowSums)
        normalArray = np.divide(array, rowSums.reshape(len(rowSums), 1))
        return normalArray
    
    def printMatrix(self, rowHeaders, colHeaders, matrix):
        """Takes row and column header lists and the matrix to be prin
ted and prints it to 3 decimal places """
        print('\t' + '\t'.join(colHeaders))
        for header, row in zip(rowHeaders, matrix):
            print('%s \t %s' % (header, '\t'.join(('%.3f' % i).rstrip('0').rstrip('.') for i in row)))

def main():

    toParse = sys.stdin.read()

    PseudoProfile = BuildHMMProfileWithPseudo()
    BuildProfile = HMMProfile()
    
    threshold, sigma, alphabet, alignArray = parseInput(toParse)

    isInsertDict, stateList, statesDict, emisDict, numInserts = BuildProfile.insertOrMatch(threshold, alignArray)
    numCols, statesDict, emisDict = BuildProfile.makeGraph(alignArray, isInsertDict, statesDict, emisDict, numInserts)
    lenStates = (numCols - numInserts)*3+3
    transitionArray = BuildProfile.buildTransition(lenStates, stateList, statesDict)
    emissionArray = BuildProfile.buildEmission(alphabet, lenStates, stateList, emisDict)

    normalizedTransition = PseudoProfile.addPseudoTrans(sigma, lenStates, transitionArray)
    normalizedEmission = PseudoProfile.addPseudoEmis(sigma, lenStates, alphabet, emissionArray)

    PseudoProfile.printMatrix(stateList, stateList, normalizedTransition)
    print('--------')
    PseudoProfile.printMatrix(stateList, alphabet, normalizedEmission)
    
if __name__== "__main__":
    main()

