#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Builds a De Bruijn graph in the form of an adjacency list given a sequence and a kMer length -k
Input: -k kMer length and a sequence from stdin
Output: Adjacency list of component kMers of sequence based on overlap
"""

import argparse
import sys

class CommandLine():
    """ Parses command line argument for kmer length"""
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-k', nargs='?', default=3, type=int, help='Set kmer length, default=3')
        self.args = self.parser.parse_args()
        
class DeBruijn:
    """ Constructs the de Bruijn graph of a string for kMer length -k """
    def __init__(self, seqStr, args):
        self.seqStr = seqStr
        self.args = args
        self.k = args.k
        
    def getkMers(self):
        """ Finds and organizes all kMers in a sequence into an adjacency list and returns as a sorted list """

        adjacencyDict = self.slideSeq()
        
        #https://stackoverflow.com/questions/32081729/how-to-sort-dictionary-on-first-element-of-the-key-tuple
#        sortedkList = sorted(adjacencyDict.items())
        sortedkDict = {x:sorted(adjacencyDict[x]) for x in adjacencyDict.keys()}
        dsortedkList = sorted(sortedkDict.items())
        return dsortedkList

    def slideSeq(self):
        """ Slides across sequence stores beginning and end of each kMer as an adjacency list """
        nodesDict = {}
        for x in range(len(self.seqStr)-self.k+1):
            pre = self.seqStr[x:self.k+x-1]
            end = self.seqStr[x+1:self.k+x]
            if pre in nodesDict.keys():
                nodesDict[pre].append(end)
            else:
                nodesDict[pre] = [end]
        return nodesDict


def main():
    """ Takes a sequence string and prints out component prefix/suffix kMers in lexicographic order """
    arguments = CommandLine()
    seqStr = sys.stdin.read().strip()
    graphKmers = DeBruijn(seqStr, arguments.args)
    graphList = graphKmers.getkMers()
    for kMer in graphList:
        print('{0} -> {1}\r'.format(kMer[0], ','.join(kMer[1])))

if __name__== "__main__":
    main()
