#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Calculates the probability of the observed sequence given a hidden path of states
Input: An observed sequence, the alphabet of the sequence, the available states, the hidden path of states, and a matrix of probabilities of the observed sequence given a hidden state
Output: The probability of the outcome of an observed sequence given a hidden path of states 
"""

import sys

def parseInput(toParse):
    """Parses the input file and returns the path of the observed sequence, the alphabet of the sequence, the path of hidden states, the available states, and the transition matrix of probabilities"""

    splitified = toParse.split('--------')
    sequence = splitified[0].rstrip().strip()
    statesPath = splitified[2].rstrip().strip()
    probMatrix = splitified[4].rstrip().strip().splitlines()
    return(sequence, statesPath, probMatrix)

class ProbabilitySeqGivenHiddenPath():
    """Takes the path of the observed sequence, hidden paths and the emission matrix to calculate the probability of the observed sequence"""
    def __init__(self):
        self.totalProbability = 1

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

    def calcProbStatePath(self, sequence, statesPath, emissionDict):
        """Calculates the probability of an observed sequence given its hidden states"""
        seqList = list(sequence)
        statesList = list(statesPath)

        for i in range(len(seqList)):
            seqState = str(statesList[i]) + str(seqList[i])
            self.totalProbability *= float(emissionDict[seqState])
            prevState = statesList[i]
        
        return self.totalProbability

def main():

    toParse = sys.stdin.read()
    sequence, statesPath, probMatrix = parseInput(toParse)
    probPath = ProbabilitySeqGivenHiddenPath()
    emissionDict = probPath.getEmissionDict(probMatrix)
    finalProb = probPath.calcProbStatePath(sequence, statesPath, emissionDict)
    print(finalProb)

if __name__== "__main__":
    main()

