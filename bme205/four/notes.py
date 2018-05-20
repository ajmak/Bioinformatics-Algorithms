#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Finds all possible kmers of length k in provided sequence and returns them in lexicographic order
Input: Motif length as a positive integer, sequence as fasta file or plain seq file
Output: kMer composition of given string 
"""

from fastaReader import FastAreader
import argparse
import sys

class CommandLine():
    """ Parses command line argument for kmer length"""
    def __init__(self):
        #adapted from David Bernick's program.py skeleton code and http://docs.python.org/2/library/argparse.html
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-k', nargs='?', default=3, type=int, help='Set kmer length, default=3')
        self.args = self.parser.parse_args()
        

class node:
    def __init__(self, data, next=None, prev=None):
        self.data = data
        self.out = {}
        self.inc = {}
    def addEdge(self, node):
        self.out[node] = (names of the node)
        node.inc[self] = 
#        for each edge that you add must tell incoming node about previous node
        #each node knows about each of its edges 
        #if looking at some node a what are its edges a.out.keys() would give you all the nodes that that node points to 
        #for every key in list print node name and all of its edge node names
class LinkedList:
        def addNode():
            pass
        def removeNode():
            pass
        def pop():
            pass

def main():
    """Accepts file containing either sequence(s) alone and prints kmers in lexicographic order"""
    arguments = CommandLine()
    seqStr = sys.stdin.readlines()

    for seq in seqStr:
        kMerSearch = FindKmers(seq, arguments.args)
        sortedList = kMerSearch.getkMers()
        for kMer in sortedList:
            print(kMer)

if __name__== "__main__":
    main()
