#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: None

"""
Determines the parameters of an HMM given the emitted sequence and hidden path and prints them as matrices

Input: An observed sequence, the alphabet of possible observed sequences, the hidden path, and the allowed states, 
Output: The transition and emission matrices outlining the parameters of the HMM
"""
import sys
import numpy as np

def parseInput(toParse):
    """Takes the input file and parses the sequence, corresponding\
 alphabet, states, and available states into lists and the """
    splitified = toParse.split('--------')
    sequence = list(splitified[0].rstrip().strip())
    alphabet = splitified[1].rstrip().strip().split()
    hiddenPath = splitified[2].rstrip().strip()
    availableStates = splitified[3].rstrip().strip().split()
    
    return(list(sequence), alphabet, hiddenPath, availableStates)
                                                        

class BuildParameters:
    """Builds/prints matrices outlining the parameters of an HMM. Consists a matrix of transition probabilities between states and a matrix of emission probabilties (probability of a sequence given its state)"""

    def makeTransMatrix(self, availableStates, hiddenPath):
        """Counts all transitions in the hiddenPath, calculates the occurence of each based on the total transitions from each state and stores them in a transition matrix """
        transCounts = {state: {inState:0 for inState in availableStates} for state in availableStates}
        transMatrix = np.zeros(shape = (len(availableStates),len(availableStates)))
        indexKeys = [i for i in range(len(availableStates))]
        indexDict = {state:i for state, i in zip(availableStates, indexKeys)}
        #populate count dictionary for all transitions
        for i in range(len(hiddenPath)-1):
            transCounts[hiddenPath[i]][hiddenPath[i+1]] += 1

            #populate transition matrix with probabilities based on count dict
        for row in transCounts.keys():
            rowNum = indexDict[row]
            sumRow = sum(transCounts[row].values())
            for col in transCounts[row].keys():
                colNum = indexDict[col]
                if sumRow == 0:
                    transMatrix[rowNum, colNum] = 1/len(availableStates)
                else:
                    transMatrix[rowNum, colNum] = transCounts[row][col]/sumRow
        #print formatted transition matrix
        return transMatrix

    def makeEmissionMatrix(self, alphabet, availableStates, sequence, hiddenPath):
        """Counts all emitted symbols in sequence for all given states, calculates the occurrence of each symbol in each state and stores them in an emission matrix """
        emisCounts = {state: {letter:0 for letter in alphabet} for state in availableStates}
        emisMatrix = np.zeros(shape = (len(availableStates),len(alphabet)))
        alphKeys = [i for i in range(len(alphabet))]
        indexKeys = [i for i in range(len(availableStates))]
        alphDict = {letter:i for letter, i in zip(alphabet, alphKeys)}
        indexDict = {state:i for state, i in zip(availableStates, indexKeys)}
        #populate emission count dictionary for each state/emitted symbol
        for i in range(len(sequence)):
            emisCounts[hiddenPath[i]][sequence[i]] += 1
        #Print probability based on count dict into matrix
        for row in emisCounts.keys():
            rowNum = indexDict[row]
            sumRow = sum(emisCounts[row].values())
            for col in emisCounts[row].keys():
                colNum = alphDict[col]
                if sumRow == 0:
                    emisMatrix[rowNum, colNum] = 1/len(alphabet)
                else:
                    emisMatrix[rowNum, colNum] = emisCounts[row][col]/sumRow
        return(emisMatrix)
        
    def printMatrix(self, rowHeaders, colHeaders, matrix):
        """Takes row and column header lists and the matrix to be printed and prints it to 3 decimal places """
        print('\t' + '\t'.join(colHeaders))
        for header, row in zip(rowHeaders, matrix):
              print('%s \t %s' % (header, '\t'.join(str(round(i,3)) for i in row)))
def main():

    toParse = sys.stdin.read()

    HMMstructures = BuildParameters()
    sequence, alphabet, hiddenPath, availableStates = parseInput(toParse)

    transMatrix = HMMstructures.makeTransMatrix(availableStates, hiddenPath)
    emissionMatrix = HMMstructures.makeEmissionMatrix(alphabet, availableStates, sequence, hiddenPath)
    HMMstructures.printMatrix(availableStates, availableStates, transMatrix)
    print('--------')
    HMMstructures.printMatrix(availableStates, alphabet, emissionMatrix)
    
if __name__== "__main__":
    main()

