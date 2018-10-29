
#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Given 4 positive integers defining two intervals and a string, this program will split the string and return the regions between both intervals
Input: 4 positive integers separated by spaces and a string
Output: 2 strings separated by a space defined by the provided intervals
"""

class StringParser():
    """ Contains command line parsing code and functions to sum the odd integers on the given interval """
    def __init__(self):
        """ Initializes argument parser """

        import sys
        import argparse
        #adapted from David Bernick's program.py skeleton code
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("sentence", type=str)
        self.parser.add_argument("a", type=int)
        self.parser.add_argument("b", type=int)
        self.parser.add_argument("c", type=int)
        self.parser.add_argument("d", type=int)
        self.args = self.parser.parse_args()
        
        if self.args.a > 0: pass
        else:
            print("Value of a must be a positive integer")
            sys.exit(1)
        if self.args.b > 0: pass
        else:
            print("Value of b must be a positive integer")
            sys.exit(1)

    def SplitStrings(self, interval):
#        print(self.args.sentence)
#        print(self.args.a)
        first = self.args.sentence[self.args.a: self.args.b + 1]
        
        second = self.args.sentence[self.args.c: self.args.d + 1]
        print('{0} {1}'.format(first, second))

def main():
    """ Retrieves command line arguments from class StringParser and passes them to the SplitStrings method """
    Odds = StringParser()
    Odds.SplitStrings(Odds.args)

if __name__== "__main__":
    main()
