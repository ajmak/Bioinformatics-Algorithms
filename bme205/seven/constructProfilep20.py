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
        alignList = [list(i+'S') for i in alignment]

        alignArray = np.array(alignList)
#        alignArray = np.append(alignArray
        print(alignArray)
        return threshold, alphabet, alignArray
    
    def makeGraph(self, threshold, alphabet, alignArray):
        """Constructs a graph as a dictionary of transition edges with the log of their edge weights"""
        
        lenEmis = len(alphabet)
        lenStates = (len(alignArray[0])-1) * 3 + 3
        transList = [[0]*lenStates for t in range(lenStates)]
        transArray = np.zeros(shape = (lenStates, lenStates))
        emissionArray = np.zeros(shape = (lenStates, lenEmis))
        print(transArray)
        emissionList = []
        #i keeps track of columns
        i  = rowT = 0
        colT = 1
        for column in alignArray.T:
            #j keeps track of rows
            j = 0
            deletes = list(column).count('-')
            emisTotal = len(column) - deletes
            tempDict = dict(prevDict)
            ratioInserts = deletes/len(column)
            if i == 0 :
                prevDict = {key: 'source' for key in range(len(column))}
                i+=1
                colT = rowT = i
                continue

            for j in range(len(column)):
                #only consider columns below insertion threshold
                prevState = prevDict[j]
                curState = column[j]
                #appends insert, match, and delete for each row
                if (ratioInserts) < float(threshold):
                    transList = self.addMatchDeletes(rowT, colT, prevState, curState, transList)
                    prevDict[j] = curState
                    prevIsInsert = False
                    
                else:
                    if prevState != '-' and curState == '-': im +=1
                    if prevState == '-' and curState != '-': mi +=1
                    if prevState != '-' and curState != '-': ii +=1
                    prevIsInsert = True
                if prevIsInsert and prevDict[j] == '-': continue
                prevDict[j] = curState
            delCount = list(tempDict.values()).count('-')
            matchCount = len(tempDict.values()) - delCount
            transList[rowT+1] = [round(m/matchCount,1) for m in transList[rowT+1]]
            transList[rowT+2] = [round(d/delCount,1) if delCount !=0 else 0 for d in transList[rowT+2]]
                

            i += 1
            if i == 1: rowT +=1; colT += 1
            if prevIsInsert != True and i != 1:
                rowT += 3
                colT += 3
                    
        print(np.array(transList))
    def addMatchDeletes(self, rowT, colT, prevState, curState, transList):
        if prevIsInsert == True:
            #im
            if prevState != '-' and curState != '-': transList[rowT][colT+1] += 1
            #id
            if prevState != '-' and curState == '-': transList[rowT][colT+2] += 1
        else:
            #mm
            if prevState != '-' and curState != '-': transList[rowT+1][colT+1] += 1
            #md
            if prevState != '-' and curState == '-': transList[rowT+1][colT+2] += 1
            #dm
            if prevState == '-' and curState != '-': transList[rowT+2][colT+1] += 1
            #dd
            if prevState == '-' and curState == '-': transList[rowT+2][colT+2] += 1
        return transList

    def addInserts(self, rowT, colT, prevState, curState, transList):
        if prevIsInsert == True:
            #ii
            if prevState != '-' and curState != '-': transList[rowT][colT] += 1
        else:
            #mi
            if prevState != '-' and curState != '-': transList[rowT+1][colT+3] += 1
            #di
            if prevState == '-' and curState != '-': transList[rowT+2][colT+3] += 1
def main():

    toParse = sys.stdin.read()

    HMMstructures = HMM()
    threshold, alphabet, alignArray = HMMstructures.parseInput(toParse)

    dagGraph = HMMstructures.makeGraph(threshold, alphabet, alignArray)
#    dagProb = DirectedAcyclicGraph()
    
#    getViterbi = dagProb.longestPath(sequence, dagGraph, transitionDict, emissionDict)
    

if __name__== "__main__":
    main()

