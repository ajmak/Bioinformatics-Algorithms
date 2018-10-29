#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""
Given 2 positive integers defining an interval, this program will return the sum of all of the odd integers within that interval, inclusive.
Input: 2 positive integers separated by a space
Output: Sum of the odd integers within that interval, inclusive
"""

class SumOddInts():
    """ Sums the odd integers on the given interval """
    def __init__(self, interval):
        import sys
        self.interval = interval
        # removed argparse and try/except for positive intervals
    def SumCalc(self):
        """ utilizes a while loop and interval of 2 to sum odd integers """
        a = int(self.interval.split()[0])
        b = int(self.interval.split()[1])
        #changed to list comprehension
        sumList = [x for x in range(a,b+1) if x%2!=0]
        total = sum(sumList)
        return total

def main():
    interval = input('Enter Interval: \n')
    Odds = SumOddInts(interval)
    getTotal=Odds.SumCalc()
    print(getTotal)
if __name__== "__main__":
    main()
