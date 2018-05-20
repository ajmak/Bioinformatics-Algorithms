#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Xian Chang, Lori Siao

"""
"""
import pdb
import sys
import math
import numpy as np


class HMM():
    """Parses the input file and builds/returns a dictionary of transition probabilities between states and a dictionary of emission probabilties (probability of a sequence given its state)"""
    def parseInput(self, toParse):
        """Parses the input file into the threshold, alphabet, and a matrix of the provided alignment"""
        splitified = toParse.split('--------')
        threshold = splitified[0].rstrip().strip()
        alphabet = splitified[1].rstrip().strip().split()
        alignment = splitified[2].rstrip().strip().splitlines()
        print(alignment)
        alignList = [list(i) for i in alignment]

        alignArray = np.array(alignList)

        print(alignArray)
        return threshold, sigma, alphabet, alignArray

    def insertOrMatch(self, threshold, alignArray):
        isInsertDict = {}
        stateList = []
        emisList = []
        i = 0
        st = 1
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
            i += 1

        stateList.insert(0,'I0')
        stateList.insert(0,'S')
        stateList.append('E')
        emisList.insert(0,'I0')
        statesDict = {states: {} for states in stateList}
        emisDict = {states: {} for states in emisList}
        return isInsertDict, stateList, statesDict, emisDict

    
                
    def makeGraph(self, alignArray, isInsertDict, statesDict, emisDict):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        i = 0
        numCols = alignArray.shape[1]
        numRows = alignArray.shape[0]
        stateNum = 0
        for index, x in np.ndenumerate(alignArray): 
            if index[1] == 0:
                stateNum = 0
                prevState = 'S'

            if isInsertDict[index[1]] == 'isMatch':
                if stateNum == 0:
                    stateNum += 1 
                if x != '-':
                    if 'M'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['M'+str(stateNum)] = 0
                    if x not in emisDict['M'+str(stateNum)].keys(): emisDict['M'+str(stateNum)][x] = 0
                    statesDict[prevState]['M'+str(stateNum)] += 1
                    emisDict['M'+str(stateNum)][x] += 1
                    prevState = 'M' + str(stateNum)
                if x == '-':
                    if 'D'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['D'+str(stateNum)] = 0
                    statesDict[prevState]['D'+str(stateNum)] += 1
                    prevState = 'D' + str(stateNum)
                stateNum += 1
            else:
                if x != '-':
                    if 'I'+str(stateNum) not in statesDict[prevState].keys(): statesDict[prevState]['I'+str(stateNum)] = 0
                    if x not in emisDict['I'+str(stateNum)].keys(): emisDict['I'+str(stateNum)][x] = 0
                    statesDict[prevState]['I'+str(stateNum)] += 1
                    emisDict['I'+str(stateNum)][x] += 1
                    prevState = 'I' +str(stateNum)
            if index[1] == numCols-1:
                if 'E' not in statesDict[prevState].keys():
                    statesDict[prevState]['E'] = 0
                statesDict[prevState]['E'] += 1
                                                                    
            i+=1
            
#            if stateNum == numCols + 1: stateNum = 0

            #print(transDict)
        print(statesDict)
        print()
        print(emisDict)
        return numCols, statesDict, emisDict
    
    def buildTransition(self, numCols, stateList, statesDict):
        lenStates = numCols*3+3
        transArray = np.zeros(shape = (lenStates, lenStates))
        indexKeys = [i for i in range(lenStates)]        
        indexDict = {state:i for state, i in zip(stateList, indexKeys)}
        print(indexDict)

        for row in statesDict.keys():
            rowNum = indexDict[row]
            sumRow = sum(statesDict[row].values())
            for col in statesDict[row].keys():
                colNum = indexDict[col]
                
                transArray[rowNum, colNum] = statesDict[row][col]/sumRow
        print(transArray)
        return transArray
    def buildEmission(self, alphabet, numCols, stateList, emisDict):
        lenEmis = len(alphabet)
        lenStates = numCols*3+3
        emissionArray = np.zeros(shape = (lenStates, lenEmis))
        alphKeys = [i for i in range(lenEmis)]
        indexKeys = [i for i in range(lenStates)]
        print(alphKeys)
        indexDict = {state:i for state, i in zip(stateList, indexKeys)}
        alphIndex = {letter:i for letter, i in zip(alphabet, alphKeys)}
        print(alphIndex)
        
        for row in emisDict.keys():
            rowNum = indexDict[row]
            sumRow = sum(emisDict[row].values())
            for col in emisDict[row].keys():
                colNum = alphIndex[col]
                emissionArray[rowNum, colNum] = emisDict[row][col]/sumRow
        print(emissionArray)
def main():

    toParse = sys.stdin.read()

    HMMstructures = HMM()
    threshold, alphabet, alignArray = HMMstructures.parseInput(toParse)
    isInsertDict, stateList, statesDict, emisDict = HMMstructures.insertOrMatch(threshold, alignArray)
    numCols, statesDict, emisDict = HMMstructures.makeGraph(alignArray, isInsertDict, statesDict, emisDict)
    transitionArray = HMMstructures.buildTransition(numCols, stateList, statesDict)
    emissionArray = HMMstructures.buildEmission(alphabet, numCols, stateList, emisDict)
#    dagProb = DirectedAcyclicGraph()
    
#    getViterbi = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
    

if __name__== "__main__":
    main()

