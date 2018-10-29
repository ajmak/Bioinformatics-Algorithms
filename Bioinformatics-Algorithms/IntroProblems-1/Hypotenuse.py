#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Calculates the square of the hypotenuse of a right triangle given the length of the other two legs
Input: 2 positive integers separated by a space
Output: The square of the hypotenuse of a right triangle as defined by its two other leg lengths
"""

class Hypotenuse():
    """ Contains command line parsing code and a function to calculate the square of the sum of two squared numbers """
    def __init__(self):
        """ Initializes argument parser """
        #Code taken from my previous program, SumOddLoop.py
        import sys
        import argparse
        #adapted from David Bernick's program.py skeleton code
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("a", type=int)
        self.parser.add_argument("b", type=int)
        self.args = self.parser.parse_args()
        
        if self.args.a > 0: pass
        else:
            print("Value of a must be a positive integer")
            sys.exit(1)
        if self.args.b > 0: pass
        else:
            print("Value of b must be a positive integer")
            sys.exit(1)

    def SquareSumCalc(self, a, b):
        """ Calculates square and sums arguments a and b """
        sqsum = self.args.a**2 + self.args.b**2
        print(sqsum)
        
def main():
    """ Retrieves command line arguments from class Hypotenuse and passes them to the SquareSumCalc method """
    Leg = Hypotenuse()
    try:
        print (Leg.args)
    except:
        sys.exit(1)
    Leg.SquareSumCalc(Leg.args.a, Leg.args.b)

if __name__== "__main__":
    main()
