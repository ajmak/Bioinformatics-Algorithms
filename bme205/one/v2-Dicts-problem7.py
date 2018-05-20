
#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Given a string this program will return a count of each word in that string
Input: a string separated by spaces
Output: each substring with the count of the substring in the given string
"""
import sys
import os.path

class DictionaryCount():
    def __init__(self, reads):
        CountDict = {}
        maintainOrder=[]
        self.maintainOrder = maintainOrder
        self.CountDict = CountDict
        self.reads = reads

    #removed argparse and unused arguments    
    def CountWords(self):
        """Counts number of times a word was used in a sentence and returns a dictionary of those counts and the order in which they were first seen """
        
        for word in self.reads.split():
            if word not in self.CountDict.keys():
                self.CountDict[word] = 1
                ##################################################
                # to address comment "what if order mattered?"   #
                # while still using Dictionary as assignment     #
                # specified. Chose to not use list comprehension #
                # to fill order so that words and their counts   #
                # would not be printed twice                     # 
                self.maintainOrder.append(word)
            else:
                self.CountDict[word] += 1
        return self.CountDict, self.maintainOrder


def main():
    """ Reads file from stdin and prints ordered words in sentence with counts """
    #changed to reading file from stdin since comment indicated accepting a string from commandline was not preferred
    reads = sys.stdin.read().rstrip()
    Counter = DictionaryCount(reads)
    countDict, order = Counter.CountWords()
    #moved output to main
    for item in order:
            print('{0}\t{1}'.format(item, countDict[item]))

if __name__== "__main__":
    main()
