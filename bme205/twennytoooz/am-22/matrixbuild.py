#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian Chang, Lori Siao

"""
Constructs a profile of probabilities of transitions between states and of symbols given states from a provided multiple alignment. Insertion states are determined by the provided threshold.
Input: Threshold, alphabet, multiple alignment
Output: Profile HMM outlined by a transition matrix of probabilities between states and an emission matrix of probabilities of emitted symbols at each state

"""
import pdb
import sys
import math
import numpy as np


def parseInput(toParse):
    """Parses the input file into the threshold, alphabet, and a matrix of the provided alignment"""
    splitified = toParse.split('--------')
    threshold = splitified[0].rstrip().strip()
    alphabet = splitified[1].rstrip().strip().split()
    alignment = splitified[2].rstrip().strip().splitlines()
    alignList = [list(i) for i in alignment]
    alignArray = np.array(alignList)
    
    return threshold, alphabet, alignArray
def printMatrix(rowHeaders, colHeaders, matrix):
    """Takes row and column header lists and the matrix to be printed and prints it to 3 decimal places """
    print('\t' + '\t'.join(colHeaders))
    for header, row in zip(rowHeaders, matrix):
        print('%s \t %s' % (header, '\t'.join(('%.3f' % i).rstrip('0').rstrip('.') if i!=1.0 else ('%.1f' % i) for i in row)))
        

class HMMProfile():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""

    def insertOrMatch(self, threshold, alignArray):
        """"Determines if each column cnsists of insert states or match/delete states and stores this as a dictionary by index. Also createsdictionary of all states to store transitions and emissions (so don't have to reiterate through columns an extra time)"""
        isInsertDict = {}
        stateList = []
        emisList = []
        i = 0
        st = 1
        numInserts = 0
        for column in alignArray.T:
            deletes = list(column).count('-')
            ratioInserts = deletes/len(column)
            if (ratioInserts) < float(threshold):
                isInsertDict[i] = 'isMatch'
                s=str(st)
                stateList.extend(('M'+s,'D'+s,'I'+s))
                emisList.extend(('M'+s,'I'+s))
                st += 1
            else:
                isInsertDict[i] = 'isInsert'
                numInserts += 1
            i += 1
        #add start, end, and initial insert states
        stateList.insert(0,'I0')
        stateList.insert(0,'S')
        stateList.append('E')
        emisList.insert(0,'I0')
        statesDict = {states: {} for states in stateList}
        emisDict = {states: {} for states in emisList}

        return isInsertDict, stateList, statesDict, emisDict, numInserts

    def makeGraph(self, alignArray, isInsertDict, statesDict, emisDict, numInserts):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        i = 0
        numCols = alignArray.shape[1]
        numRows = alignArray.shape[0]
        stateNum = 0

        for index, x in np.ndenumerate(alignArray):
            if index[1] == 0: #add start state
                stateNum = 0
                prevState = 'S'

            if isInsertDict[index[1]] == 'isMatch': #in match/delete column
                stateNum += 1 
                if x != '-': #if match
                    if 'M'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['M'+str(stateNum)] = 0
                    if x not in emisDict['M'+str(stateNum)].keys(): emisDict['M'+str(stateNum)][x] = 0
                    statesDict[prevState]['M'+str(stateNum)] += 1
                    emisDict['M'+str(stateNum)][x] += 1
                    prevState = 'M' + str(stateNum)
                if x == '-': #if delete
                    if 'D'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['D'+str(stateNum)] = 0
                    statesDict[prevState]['D'+str(stateNum)] += 1
                    prevState = 'D' + str(stateNum)
            else: #in insert column
                if x != '-': #if is insert
                    if 'I'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['I'+str(stateNum)] = 0
                    if x not in emisDict['I'+str(stateNum)].keys(): emisDict['I'+str(stateNum)][x] = 0
                    statesDict[prevState]['I'+str(stateNum)] += 1
                    emisDict['I'+str(stateNum)][x] += 1
                    prevState = 'I' +str(stateNum)
            if index[1] == numCols-1: #add end state
                if 'E' not in statesDict[prevState].keys():
                    statesDict[prevState]['E'] = 0
                statesDict[prevState]['E'] += 1
                                                                    
            i+=1
        return numCols, statesDict, emisDict
    
    def buildTransition(self, lenStates, stateList, statesDict):
        """Takes count matrix of transitions and returns a matrix of transition probabilities between states """
        
        transArray = np.zeros(shape = (lenStates, lenStates))
        indexKeys = [i for i in range(lenStates)]        
        indexDict = {state:i for state, i in zip(stateList, indexKeys)}

        for row in statesDict.keys():
            rowNum = indexDict[row]
            sumRow = sum(statesDict[row].values())
            for col in statesDict[row].keys():
                colNum = indexDict[col]
                transArray[rowNum, colNum] = statesDict[row][col]/sumRow
        return transArray
    def buildEmission(self, alphabet, lenStates, stateList, emisDict):
        """Takes count matrix of emissions and returns a matrix of emission probabilities at each state/symbol """
        
        lenEmis = len(alphabet)
        emissionArray = np.zeros(shape = (lenStates, lenEmis))
        alphKeys = [i for i in range(lenEmis)]
        indexKeys = [i for i in range(lenStates)]
        indexDict = {state:i for state, i in zip(stateList, indexKeys)}
        alphIndex = {letter:i for letter, i in zip(alphabet, alphKeys)}
        
        for row in emisDict.keys():
            rowNum = indexDict[row]
            sumRow = sum(emisDict[row].values())
            for col in emisDict[row].keys():
                colNum = alphIndex[col]
                emissionArray[rowNum, colNum] = emisDict[row][col]/sumRow
        return(emissionArray)
    
def main():

    toParse = sys.stdin.read()

    BuildProfileHMM = HMMProfile()
    threshold, alphabet, alignArray = parseInput(toParse)
    isInsertDict, stateList, statesDict, emisDict, numInserts = BuildProfileHMM.insertOrMatch(threshold, alignArray)
    numCols, statesDict, emisDict = BuildProfileHMM.makeGraph(alignArray, isInsertDict, statesDict, emisDict, numInserts)
    lenStates = (numCols - numInserts)*3+3
    transitionArray = BuildProfileHMM.buildTransition(lenStates, stateList, statesDict)
    emissionArray = BuildProfileHMM.buildEmission(alphabet, lenStates, stateList, emisDict)

    printMatrix(stateList, stateList, transitionArray)
    print('--------')
    printMatrix(stateList, alphabet, emissionArray)

if __name__== "__main__":
    main()

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

