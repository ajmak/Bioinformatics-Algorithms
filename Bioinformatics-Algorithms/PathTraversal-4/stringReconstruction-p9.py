#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Recontructs a string from the overlap of its prefix/suffix k-mers.
Input: A prefix/suffix genome path given in order as a list of k-mers
Output: The string constructed from the given k-mers
"""

import sys

class ReconstructString:
    """ Reconstructs a string based on its component kMers """
    def __init__(self, kMerList):
        self.kMerList = kMerList
        
    def joinPath(self):
        """ Given a list of kMers, appends ordered genome graph into a full sequence  """
        sequence = ''
        # appends first letter of each kMer until last kMer and appends entire last kMer
        for node in self.kMerList:
            if node == self.kMerList[-1]:
                sequence += node
            else:
                chars = list(node)
                sequence += chars[0]
        
        return sequence


def main():
    """ Accepts file containing kMers ordered in the order of the sequence """
    kMers = sys.stdin.read().splitlines()
    kMerConstruct = ReconstructString(kMers)
    fullPath = kMerConstruct.joinPath()
    
    print(fullPath)
if __name__== "__main__":
    main()
