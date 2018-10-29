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

    def makeTransitionDictforMatrix(self, alignArray, isInsertDict, statesDict, emisDict, numInserts):
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
    
    def buildTransitionMatrix(self, lenStates, stateList, statesDict):
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
    
    def buildEmissionMatrix(self, alphabet, lenStates, stateList, emisDict):
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

    def buildTransitionEmissionwithPseudo(self, alphabet, threshold, sigma, alignArray):
        isInsertDict, stateList, statesDict, emisDict, numInserts = self.insertOrMatch(threshold, alignArray)
        numCols, statesDict, emisDict = self.makeTransitionDictforMatrix(alignArray, isInsertDict, statesDict, emisDict, numInserts)
        lenStates = (numCols - numInserts)*3+3
        transitionArray = self.buildTransitionMatrix(lenStates, stateList, statesDict)
        emissionArray = self.buildEmissionMatrix(alphabet, lenStates, stateList, emisDict)

        #add pseudo counts
        normalizedTransition = self.addPseudoTrans(sigma, lenStates, transitionArray)
        normalizedEmission = self.addPseudoEmis(sigma, lenStates, alphabet, emissionArray)

        return lenStates, stateList, normalizedTransition, normalizedEmission




