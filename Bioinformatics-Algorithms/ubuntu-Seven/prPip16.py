#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Calculates the probability of a hidden path of states 
Input: A path of states and a matrix of transition probabilities between all states
Output: The probability of the given path based on the provided transition probabilities 
"""

import sys

def parseInput(toParse):
    """Parses the input file and returns the path of hidden states, the available states, and the transition matrix of probabilities"""
    splitified = toParse.split('--------')
    statesPath = splitified[0].rstrip().strip()
    availableStates = len(splitified[1].rstrip().strip().split())
    probMatrix = splitified[2].rstrip().strip().splitlines()

    return(probMatrix, statesPath, availableStates)

class ProbabilityHiddenPath:
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
         
    def calcProbStatePath(self, statesPath, transitionDict, availableStates):
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

def main():

    toParse = sys.stdin.read()
    probMatrix, statesPath, availableStates = parseInput(toParse)
    probPath = ProbabilityHiddenPath()
    transitionDict = probPath.getTransitionDict(probMatrix)
    finalProb = probPath.calcProbStatePath(statesPath, transitionDict, availableStates)
    print(finalProb)

if __name__== "__main__":
    main()

