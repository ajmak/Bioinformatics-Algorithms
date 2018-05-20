#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Lori Siao, Xian Chang, Ioannis Anastopulos, Alexis Thornton, Balaji Sundararaman, Vinay Poodari

"""

Input: 
Output:
"""
import sys
import argparse
import math
import string
import random
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
                    yield sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield sequence


class CommandLine():
    """ Parses command line arguments using argparse, both required and optional, and returns them to main"""
    def __init__(self):
        """ Initializes argument parser for required and optional arguments and passes them to main"""

        #adapted from David Bernick's program.py skeleton code and http://docs.python.org/2/library/argparse.html
        
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-i', nargs='?', default = 100000, type = int, choices = range(3,7), help = 'Set number of randomized iterations, default=100,000')
        self.parser.add_argument('-p', nargs='?', default = 1, type = int, help = 'Set number of pseudocounts, default=1')
        self.parser.add_argument('-k', nargs='?', default = 13, type = int, help = 'Set k-mer length, default=13')

        self.args = self.parser.parse_args()



class MotifSearch:
    """ """
    def __init__(self, arguments, seqStr):
        """Instantiates necessary arguments, dictionaries, lists, and objects needed for MotifSearch class """
        self.arguments = arguments
        self.seqStr = seqStr
        self.pseudo=self.arguments.p
        self.k=self.arguments.k
        self.iterations=self.arguments.i
        self.nucs=['A', 'C', 'G', 'T']
        self.kmerList=[]
        self.seqCount = 0
        self.profile = {key: [self.pseudo]*(self.k+1) for key in self.nucs}
        print(self.profile)
        #change pseudocount to pseudo*4 + k
    def genMotif(self):
        """ """
        self.seqCount = len(self.seqStr) +4*self.pseudo + self.k
        for seq in self.seqStr:
            fastLength, kmer = self.pickRandom(seq)
            self.kmerList.append(kmer)
        #self.profile = {key: [self.pseudo]*(self.k+1) for key in self.nucs}
        prevProfile = self.makeProfile(self.kmerList, self.seqCount, self.profile)
#        prevEntropy = self.findE(profile)

#        self.newProfile = {key: [x/self.seqCount for x in self.profile[key]] for key in self.profile.keys()}
        
        for seq in self.seqStr:
            self.slideSeq(seq)

#        profile = self.makeProfile(seqCount, kList)

    def pickRandom(self, seq):
        """ """
        fastLength = len(seq)
        randStart = random.randint(0, len(seq) - self.k + 1)
        print(randStart)
        kmer = seq[randStart: randStart + self.k]
        return fastLength, kmer
            
    def makeProfile(self, kmerList, seqCount, profile):
        """ """
        i=0
        profile = {key: [self.pseudo]*(self.k+1) for key in self.nucs}
        for kmer in kmerList:
            for char in kmer:
                profile[char][i] += 1
                i+=1

        self.newProfile = {key: [x/seqCount for x in profile[key]] for key in self.profile.keys()}
        
        print(kmer)
        print(self.newprofile)

    def slideSeq(self, seq):
        """ """
        
        lenSeq = len(seq)
        for i in range(lenSeq-self.k+1):
            motif = seq[i:i+self.k]
            t=0
            highProb=0
            tempProb=1
            for char in motif:
                #cange back to probability
                prevProb = tempProb
                tempProb *= self.newProfile[char][t]
                t+=1
            #    tempE += self.findE(tempProb)
            
            if tempProb > highProb:
                highProb = tempProb
                highestMotif = motif
            
        return highestMotif, highestProb    
#        print(highestMotif)
#        print(highProb)
        #don't recalculate the score - store the score --> feed this back as dictionary with score and prob, then generate profile and compare to previous profile 
    def findE(self, prob):
        """ """
        return -1*(prob*math.log(prob,2))
        
    def outToFile(self, outTuples):
        """ """
        pass

def main():
    """ Retrieve command line arguments from class CommandLine and passes them to MotifSearch method, passes each sequence from FastAReader to MotifSearch """
    seqList=[]
    fastLength=0
    seqCount=0
    arguments = CommandLine()
    readFast = FastAreader()
    seqStr = readFast.readFasta()
    for seq in seqStr: seqList.append(seq)
    getMotif = MotifSearch(arguments.args, seqList)
    genProfile = getMotif.genMotif()


    
#        fastLength = fastLength + tempLength

#    getCounts = MotifSearch(arguments.args, seqs, motifDict, fastLength)
    
#    doCalcs = getCounts.computeSeq(fastLength, motifDict)
#    sortAndPrint = getCounts.outToFile(doCalcs)

if __name__== "__main__":
    main()
