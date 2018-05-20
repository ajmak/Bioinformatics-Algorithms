#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Creates an adjacency list in the form of a De Bruijn graph from a list of kMers. The De Bruijn graph illustrates the overlap of the given kMers
Input: A list of all possible kMers from a sequence from stdin
Output: A De Bruijn graph outlining the overlap of the kMers to stdout
"""

import sys

class DeBruijn:
    """ Constructs the de Bruijn graph of a list of kMers """
    def __init__(self, kMers):
        self.kMers = kMers
 
    def buildaBruijn(self):
        """ Finds and organizes all kMers in a list into an adjacency list and returns as a sorted list """
        nodesDict = {}
        for x in self.kMers:
            pre = x[:-1].strip()
            end = x[1:].strip()
            if pre in nodesDict.keys():
                nodesDict[pre].append(end)
            else:
                nodesDict[pre] = [end]
        #sorts lists within dictionary for nodes with multiple edges
        sortnodesDict = {x:sorted(nodesDict[x]) for x in nodesDict.keys()}
        #sorts output list based on prefix node
        orderList = sorted(sortnodesDict.items())
        return orderList, sortnodesDict


def main():
    """ Given a list of kMers comprising a sequence, prints the prefix/suffix nodes showing the relationship between the nodes in order """
    kList = sys.stdin.read().splitlines()
    graphKmers = DeBruijn(kList)
    deBruijnList, bruijnDict = graphKmers.buildaBruijn()
    for kMer in deBruijnList:

        print('{0} -> {1}'.format(kMer[0], ','.join(kMer[1])))

if __name__== "__main__":
    main()

