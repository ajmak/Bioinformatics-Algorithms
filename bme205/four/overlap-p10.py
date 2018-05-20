#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Builds an overlap graph in the form of an adjacency list based on component kMers of a sequence
Input: List of kMers from a sequence
Output: A lexicographically ordered adjacency list of the component kMers based on their overlap
"""

import sys

class OverlapGraph:
    """ Reconstruct sequence based on overlap of given list of kMers """
    def __init__(self, kMerList):
        self.kMerList = kMerList
        self.edgesDict = {x:None for x in self.kMerList}

    def setAdjacency(self):
        """ Compares the suffix of every provided kMer to the prefix of every kMer and returns a dictionary of adjacency based on sequence overlap and an ordered list of kMers for printing """
        for key in self.edgesDict.keys():
            keyEnd = key[1:]
            for x in self.kMerList:
                kPre = x[:-1]
                if kPre == keyEnd:
                    self.edgesDict[key] = x
        sortedkList = sorted(self.edgesDict.keys())
        return sortedkList, self.edgesDict


def main():
    """ Takes input file of list of ordered kMers from stdin and prints the overlap graph in the form of an adjacency list in lexicographic order"""
    kMerList = sys.stdin.read().splitlines()
    graphKmers = OverlapGraph(kMerList)
    orderedList, adjacencyDict = graphKmers.setAdjacency()
    for kMer in orderedList:
        if adjacencyDict[kMer] == None: continue
        else:
            print('{0} -> {1}'.format(kMer, adjacencyDict[kMer]))

if __name__== "__main__":
    main()
