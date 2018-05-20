#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian 

"""
Incorporates pseudocounts into the probabilities of all biologically permissible transitions between states and emissions of symbols at a given state and normalizes
Input: Threshold, pseudocount(sigma), alphabet, multiple alignment
 Output: Normalized transition matrix of probabilities with pseudocounts, Normalized emission matrix of probabilities with pseudocounts

Dependencies:
problem20.py
numpy
"""

import sys
import math
import numpy as np
from problem20 import HMMProfile

def parseInput(toParse):
    """Parses the input file into the threshold, pseudocount(sigma),  alphabet, and a matrix of the provided alignment
Input format: 
0.358   0.01
--------
A   B   C   D   E
--------
ADA
ADA
AAA
ADC
-DA
D-A

"""
    splitified = toParse.split('--------')
    threshold, sigma = splitified[0].rstrip().strip().split()
    alphabet = splitified[1].rstrip().strip().split()
    alignment = splitified[2].rstrip().strip().splitlines()
    alignList = [list(i) for i in alignment]
    
    alignArray = np.array(alignList)
    
    return threshold, sigma, alphabet, alignArray
                                                                

class BuildHMMProfileWithPseudo():
    """Adds pseudocounts to the probabilities of all biologically possible transitions between states and the probabilities of all states that can emit a symbol in the alignment"""

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

