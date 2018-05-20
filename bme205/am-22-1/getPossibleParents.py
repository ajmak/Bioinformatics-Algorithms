#File: TopologicalOrderGenerator.py
#Authors: CammyG and AMakiZZle
#Description: Generates a topological order for the viterbi algorithm

import sys
        
def GetPossibleParents(currentState, currentStateIndex):
    if currentStateIndex == 0:
        if int(currentState[1:]) == 1:
            return ['S']
        else:
            return [('D' + str(int(currentState[1:])-1))]
            
    if currentState == 'I0':
        if currentStateIndex > 1:
            return ['I0']
        return ['S', 'I0']
        
    if currentState == 'M1' or currentState == 'D1':
        return ['S', 'I0']
        
    if currentState == 'E':
        return ['M'+str(currentStateIndex), 'D'+str(currentStateIndex), 'I'+str(currentStateIndex)]
        
    idx = int(currentState[1:])
    if currentState.startswith('I'):
        return ['M'+str(idx), 'D'+str(idx), 'I'+str(idx)]
    else:
        return ['M'+str(idx-1), 'D'+str(idx-1), 'I'+str(idx-1)]
        
def main():
    print('E(9)')
    print(GetPossibleParents('E',8))
    print('M2(5)')
    print(GetPossibleParents('M2',5))
    print('D2(4)')
    print(GetPossibleParents('D2',4))
    print('I2(3)')
    print(GetPossibleParents('I2',3))
    print('I0(5)')
    print(GetPossibleParents('I0',5))
    print('I0(1)')
    print(GetPossibleParents('I0',1))
    print('I8(2)')
    print(GetPossibleParents('I8',2))
    

if __name__== "__main__":
    main()

