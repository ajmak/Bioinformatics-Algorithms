#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Takes a file containing multiple lines of strings and returns the even lines of strings
Input: A file of strings with maximum of 1000 lines from stdin
Output: A file containing only the even numbered lines of the input file
"""
import sys
class evenRead():
    """ Extracts even lines from input file and prints to output file """
    def __init__(self, lines):
      
        self.lines = lines
        #removed argparse
        outfile = open('EvenReads.txt', 'w')
        self.outfile = outfile
    def ExtractLines(self):
        """ uses a for loop and interval of 2 to extract even lines and add to output file """
        i=0
        #removed exception handling
        for line in self.lines:
            if i%2 != 0:
                self.outfile.write(line)
                i+=1
            else: i+=1
        
def main():
    """Takes file from stdin """
    lines = sys.stdin
    Evens = evenRead(lines)

    Evens.ExtractLines()

if __name__== "__main__":
    main()
