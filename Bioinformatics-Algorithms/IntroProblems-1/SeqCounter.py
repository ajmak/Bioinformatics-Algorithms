#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
This program counts the number of each nucleotide in a given string.
Input: A DNA string with max length of 1000
Output: The count of each nucleotide in that string in the order of 'A', 'C', 'G', 'T' separated by spaces
"""
import argparse
import sys

class nucleotideCounter (str):
    def __init__(self, s):
        self.s = s
    def count (self, s):
        """ counts number of each nucleotide in given DNA string """
        if len(s) <= 1000: pass
        else: 
            print("DNA string exceeds maximum of 1000 nt input") 
            sys.exit(1)
        """To ensure string is all uppercase """
        s = s.upper()
        A = s.count("A")
        C = s.count("C")
        G = s.count("G")
        T = s.count("T")
        print('{0} {1} {2} {3}'.format(A,C,G,T))


def main():
    """
    Allows functionality to be imported as a module. 
    Source: https://stackoverflow.com/questions/4041238/why-use-def-main
    """
    try:
        # Ensures that DNA string is provided as input
        dna = sys.argv[1]
    except:
        print ("DNA string input is required")
        sys.exit(1)
    x = nucleotideCounter(dna)
    x.count(dna)

if __name__== "__main__":
    main()
