#File: TopologicalOrderGenerator.py
#Authors: CammyG and AMakiZZle
#Description: Generates a topological order for the viterbi algorithm

import sys

def GetNextTopologicalState(currentState, currentStateIndex, maxStateIndex):
    if currentStateIndex == 0:       
        if currentState == 'S':
            return 'D1', currentStateIndex
        elif currentState.startswith('D'):
            if int(currentState[1:]) >= maxStateIndex:
                return 'I0', currentStateIndex
            return ('D' + str(int(currentState[1:])+1)), currentStateIndex
        else:
            currentStateIndex+=1
            return 'M1', currentStateIndex
    
    idx = int(currentState[1:])
    if currentState.startswith('M'):
        return ('D' + str(idx)), currentStateIndex
        
    if currentState.startswith('D'):
        return ('I' + str(idx)), currentStateIndex
        
    if currentState.startswith('I'):
        currentStateIndex+=1
        idx+=1
        if(idx > maxStateIndex):
            return 'E', currentStateIndex
        return ('M' + str(idx)), currentStateIndex

def main():

    currentState = 'S'
    currentStateIndex = 0
    maxStateIndex = 8
    
    # currentState index is 0 for S through I0
    i=0
    while i <= maxStateIndex+1:
        print(str(currentStateIndex) + ' ' + currentState)
        currentState,currentStateIndex = GetNextTopologicalState(currentState, currentStateIndex, maxStateIndex)
        i+=1
    
    #loop until we hit 'E'
    while (currentState != 'E'):
        currentState,currentStateIndex = GetNextTopologicalState(currentState, currentStateIndex, maxStateIndex)
        print(str(currentStateIndex) + ' ' + currentState)

if __name__== "__main__":
    main()

