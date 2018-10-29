#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Finds all possible kmers of length k in provided sequence and returns them in lexicographic order
Input: Motif length as a positive integer, sequence as fasta file or plain seq file
Output: kMer composition of given string 
"""

import argparse
import sys

class CommandLine():
    """ Parses command line argument for kmer length"""
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-k', nargs='?', default=3, type=int, help='Set kmer length, default=3')
        self.args = self.parser.parse_args()

class FindKmers:
    """ Finds kMers of length k in sequence and prints them in lexicographic order """
    def __init__(self, seqStr, args):
        self.seqStr = seqStr.strip()
        self.args = args
        self.k = args.k
        
    def getkMers(self):
        """ Finds and returns all kMers in a sequence as a sorted list """
        seqLength = len(self.seqStr)
        kList = [self.seqStr[x:self.k+x] for x in range(len(self.seqStr)-self.k+1)]
        #https://stackoverflow.com/questions/7371935/sort-a-string-in-lexicographic-order-python
        for line in kList: print(line)
        sortedkList = sorted(kList, key=str.upper)
        #return sortedkList


def main():
    """Accepts file containing a sequence and prints component kmers in lexicographic order"""
    arguments = CommandLine()
    seqStr = sys.stdin.readlines()

    for seq in seqStr:
        kMerSearch = FindKmers(seq, arguments.args)
        sortedList = kMerSearch.getkMers()
#    for kMer in sortedList:
#        print(kMer)

if __name__== "__main__":
    main()
