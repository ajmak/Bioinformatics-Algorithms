#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Lori Siao, Xian Chang, Ioannis Anastopulos

"""
Conducts a randomized motif search to identify most common motifs across sequences. Uses a profile matrix and calculates entropy on that matrix. The lowest entropy matrix found is used to determine the consensus motif.
Input: FASTA sequences in stdin file, k-mer length, # of iterations, pseudocounts, shuffling option, output motif option
Output: The highest probability motif across all sequences given the lowest entropy profile
"""
import sys
import os
import pdb
import argparse
import math
import string
import random
from random import shuffle
class FastAreader():
    """
    Takes input file from standard in and returns a generator of each sequence in the given FASTA file 
    Input: Fetches from stdin
    Output: Sequence as generator 

    Class written by David Bernick before minor adjustments to the yield statement 
    """
    

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        """Searches for file input in stdin and returns to readFasta method"""
        if self.fname is '':
            return sys.stdin
        else:
            return self.fname
 
    def readFasta(self):
        """Calls doOpen function to open FASTA file, parses it at each FASTA header, and returns each sequence as a generator to main"""
        header = ''
        sequence = ''

        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence


class CommandLine():
    """ Parses command line arguments using argparse, both required and optional, and returns them to main"""
    def __init__(self):
        """ Initializes argument parser for required and optional arguments and passes them to main"""

        #adapted from David Bernick's program.py skeleton code and http://docs.python.org/2/library/argparse.html
        
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-i', nargs='?', default = 100000, type = int, help = 'Set number of randomized iterations, default=100,000')
        self.parser.add_argument('-p', nargs='?', default = 1, type = int, help = 'Set number of pseudocounts, default=1')
        self.parser.add_argument('-k', nargs='?', default = 13, choices = range(0,51), type = int, help = 'Set k-mer length, default=13')
        self.parser.add_argument('-r', dest='r', action='store_true', help = 'Run on shuffled sequences, default = False')
        self.parser.add_argument('-m', dest='m', action='store_true', help = 'Print motifs that were used to generate the profile, default = False')
        self.parser.set_defaults(r=False)
        self.parser.set_defaults(m=False)
        self.args = self.parser.parse_args()



class MotifSearch:
    """Finds highest probability motif across multiple sequences using a profile of probabilities """
    def __init__(self, arguments, seqStr, headerList):
        """Instantiates necessary arguments, dictionaries, lists, and objects needed for MotifSearch class """
        self.arguments = arguments
        self.seqStr = seqStr
        self.headerList = headerList
        self.pseudo=self.arguments.p
        self.k=self.arguments.k 
        self.iterations=self.arguments.i
        self.r = self.arguments.r
        self.m = self.arguments.m
        self.nucs=['A', 'C', 'G', 'T']
        self.profileDict={}
        self.seqCount = 0
        self.seqCount = len(self.seqStr) +4*self.pseudo + self.k

    def genMotif(self):
        """Conducts a randomized motif search that restarts i times and returns the motif of highest probability and the entropy of the profile that generated the motif """

        if self.r == True: self.seqStr = [self.shuffleSeq(word) for word in self.seqStr]
        finalMotifList=[]

        for j in range(self.iterations):
            kmerList=[]
            for seq in self.seqStr:
                fastLength, kmer = self.pickRandom(seq)
                kmerList.append(kmer)
            
            prevProfile = self.makeProfile(kmerList, self.seqCount)
            firstEntropy = self.findE(prevProfile)
            prevEntropy = 2*self.k

            while True:
                highMotifs=[]
                for seq in self.seqStr:
                    highMotif, highProb = self.slideSeq(seq, prevProfile)
                    highMotifs.append(highMotif)

                curProfile = self.makeProfile(highMotifs, self.seqCount)
                curEntropy = self.findE(curProfile)
                # stop condition, stores profile and entropy in dict if entropy is no longer decreasing
                #CHANGED: now stores previous entropy instead of current
                if curEntropy >= prevEntropy:
#                    finalMotifList.append(highMotifs)
                    self.profileDict[prevEntropy] = prevProfile
                    break

                prevEntropy = curEntropy
                prevProfile = curProfile


            profileList = list(self.profileDict.keys())
            #This accounts for when only one profile is stored
            if len(profileList) == 1: 
                continue
            #removes profile and motif associated with the higher entropy
            #list size should only ever be 1 or 2
            else:
                if profileList[0] <= profileList[1]:
                    print(profileList[0])
                    self.profileDict.pop(profileList[1])
                    del profileList[1]
                    del finalMotifList[1]
                else:
                    self.profileDict.pop(profileList[0])
                    del profileList[0]
                    del finalMotifList[0]

        #sets back to regular list since no longer need list of lists
        print(finalMotifList)
        finalMotifList = finalMotifList[0]
        finalMotif = self.regenerateMotif(self.profileDict[profileList[0]])
        print('{0}\t{1}'.format(finalMotif, profileList[0]))

        if self.m == True:
            for head in range(len(self.headerList)-1):
                print('{0}\t{1}'.format(self.headerList[head].split()[0], finalMotifList[head]))

    def shuffleSeq(self, seq):
        """Shuffles a given sequence while maintaining composition"""
        #https://stackoverflow.com/questions/976882/shuffling-a-list-of-objects
        word = list(seq)
        word = shuffle(word)

        return ''.join(word)

    def pickRandom(self, seq):
        """Picks a kmer of length k at random from a provided sequence and returns the sequence length and the kmer """
        fastLength = len(seq)
        randStart = random.randint(0, len(seq) - self.k)
        kmer = seq[randStart: randStart + self.k]

        return fastLength, kmer
            
    def makeProfile(self, kmerList, seqCount):
        """Generates a probability profile given a list of kmers """

        profile = {key: [self.pseudo]*(self.k) for key in self.nucs}

        for kmer in kmerList:
            i=0
            for char in kmer:

                profile[char][i] += 1
                i+=1

        #newProfile accounts for pseudocounts and generates the probabilities for the profile
        newProfile = {key: [x/seqCount for x in profile[key]] for key in profile.keys()}
        return(newProfile)

    def slideSeq(self, seq, prevProfile):
        """Calculates the probability of every motif of length k in a fasta sequence given by a profile matrix and returns the motif with the highest probability """
        highestMotif=''
        highProb=0
        lenSeq = len(seq)
        for i in range(lenSeq-self.k+1):
            motif = seq[i:i+self.k]
            t=0
            tempProb=1

            for char in motif:
                prevProb = tempProb
                tempProb *= prevProfile[char][t]
                t+=1

            if tempProb > highProb:
                highProb = tempProb
                highestMotif = motif

        return highestMotif, highProb    

    def findE(self, profile):
        """Calculates the entropy of a provided profile of probabilities """        
        sumE=0
        Eprof={}
        #Accounts for if psuedocount=0 so log2(0) does not error out
        Eprof = {key: [-1*(x*math.log(x,2)) if x>0.0 else 0 for x in profile[key]] for key in profile.keys()}

        for vals in Eprof.values():
            sumE+=sum(vals)

        return sumE

        
    def regenerateMotif(self, profileDict):
        """Takes final profile and returns highest probability motif """
        A = list(profileDict['A'])
        C = list(profileDict['C'])
        G = list(profileDict['G'])
        T = list(profileDict['T'])
        indexDict = {0:'A', 1:'C', 2:'G', 3:'T'}
        finalMotif=''
        for i in range(self.k):
            findMax = max([A[i], C[i], G[i], T[i]])
            maxIndex = [A[i], C[i], G[i], T[i]].index(findMax)
            finalMotif += indexDict[maxIndex]
        return finalMotif

def main():
   
    seqList=[]
    headerList=[]
    fastLength=0
    seqCount=0

    arguments = CommandLine()
 
    if arguments.args.i <=0:
        print("Cannot iterate 0 times")
        sys.exit(0)
    readFast = FastAreader()
    seqStr = readFast.readFasta()


    for seq in seqStr: 
        seqList.append(seq[1])
        headerList.append(seq[0])
#    sys.stdin='/dev/tty'

    getMotif = MotifSearch(arguments.args, seqList, headerList)
    genProfile = getMotif.genMotif()


if __name__== "__main__":
    main()
