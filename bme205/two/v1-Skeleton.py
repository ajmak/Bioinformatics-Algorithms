#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""

TODO: write usage error for if k-mer length < 2
Provides comparison of expected counts for all possible motifs of given length (default 3 to 8) to actual counts in genome. This information can be used to help determine underrepresented sequences in the given genome.
Input: FASTA file for a single genome 
Output: All possible motifs with corresponding counts from the genome, the theoretical expected counts for that motif, and the Z-scores for comparison with the null model
"""
import sys
import argparse
import math

class FastAreader():
    """
    Takes input file parsed by the CommandLine class and returns a generator of each sequence in the given FASTA file 

    Class written by David Bernick before minor adjustments to the return statement 
"""
    

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        """Searches for file input in stdin if not provided as instance to the constructor and returns or returns file provided to the constructor"""
        if self.fname is '':
            return sys.stdin
        else:
            return self.fname
            #return open(self.fname)
 
    def readFasta(self):
        """Calls doOpen function to open FASTA file, parses it at each FASTA header, and returns each sequence as a generator to main """
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



        #adapted from David Bernick's program.py skeleton code
        
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('--minMotif', nargs='?', default = 3, type = int, help = 'Set minimum motif length, default=3')
        self.parser.add_argument('--maxMotif', nargs='?', default = 8, type = int, help = 'Set maximum motif length, default=8')
        #Change this later!!!
        self.parser.add_argument('--cutoff', nargs='?', default = 0, type = int, help = 'Set minimum z cutoff score, default=0')
        self.parser.add_argument("infile", type=argparse.FileType('r'))
#        self.parser.add_argument("outfile", type=argparse.FileType('w'))
        self.args = self.parser.parse_args()
        #http://docs.python.org/2/library/argparse.html#type


class StoreCounts:
    """Counts all possible motifs in provided range,  """
    def __init__(self, arguments, seq, motifDict, fastLength):
        self.arguments = arguments
        self.seq = seq
        self.motifDict = motifDict
#        self.motifDict={}
        self.mini=self.arguments.minMotif
        self.maxi=self.arguments.maxMotif
        z=self.arguments.cutoff
        self.midmin = self.mini - 2
        self.outTuple = ()        
        self.indexDict={'A':0,'C':1,'G':2,'T':3}
        self.intDict= {0:'A', 1:'C', 2:'G', 3:'T'}
        self.stop = 0
        

    def slideSeq(self):
        print(self.seq)
        allowed=set('ACGT')
        fastLength = len(self.seq)
        j=0
        
        for i in self.seq[j:].rstrip():

            for k in range(self.midmin, self.maxi+1):
                
                if j >= fastLength-k+1:
                    break
                self.motif = self.seq[j:j+k]
                index = self.motif
                #index = self.seqToIndex(self.motif)
                if set(self.motif) <= allowed:

                    if index not in self.motifDict.keys():
                        self.motifDict[index] = 1 
                    else: self.motifDict[index] += 1
                else:
                    break
            j+=1
        self.motifDict['None'] = 1
        
#        print(self.motifDict)

        return fastLength, self.motifDict
    
    def computeSeq(self, fastLengths, motifDict):
        nextSwitch=0
        maxIndex = (4**self.maxi + 4**(self.maxi -1) )
        if self.mini == 2:minIndex = 4**(self.mini - 1)
        else: minIndex = 4**(self.mini - 1) + 4**(self.mini -2)
#        print(minIndex)
        word = self.mini
        print(str(word) + 'word')
        prevSwitch = minIndex
        for i in range(1,self.mini + 1): nextSwitch += 4**i
        #nextSwitch = (4**word + 4**(word -1) )
        print(str(nextSwitch) + ' nextSwitch')
        for i in range(minIndex, maxIndex):
            print(i)
            if i == nextSwitch: 
                word+=1
                
                print(str(i) + ' switch')
                prevSwitch = nextSwitch
                print(str(prevSwitch) + ' prevSwitch')
                nextSwitch = nextSwitch + 4**word #4**word + 4**(word -1) 
                #prevSwitch = 4**(word - 1) + 4**(word -2)
            seqInt = i - prevSwitch
            motifAtIndex = self.indexToSeq(seqInt, word)
            print(motifAtIndex + str(seqInt))
        
    def indexToSeq(self, seqInt, k):
        if k==1:
            #returns first letter of string
            return self.intDict[seqInt]
        prefixIndex = seqInt//4
        remainder = seqInt%4
#        print((self.indexToSeq(prefixIndex, k-1)) + self.intDict[remainder])
        return (self.indexToSeq(prefixIndex, k-1)) + self.intDict[remainder]

    def expectedVal(self, motif):
        if motif[:-1] not in self.motifDict.keys(): left='None'
        else: left=motif[:-1]
        if motif[1:] not in self.motifDict.keys(): right='None'
        else: right=motif[1:]
        if motif[1:-1] not in self.motifDict.keys():
            expected = 'None'
            pass
        else: 
            mid=motif[1:-1]

#            print(motif + ' motif')
#            print(left + ' left')
#            print(right + ' right')
#            print(mid + ' mid')
#            print(self.motifDict[motif])
#            print(self.motifDict[left])
#            print(self.motifDict[right])
#            print(self.motifDict[mid])

            expected = (self.motifDict[left] * self.motifDict[right]) / self.motifDict[mid]
        return expected

    def zCalc(self, fastLength):
        for key in self.motifDict.keys():
            exp = self.expectedVal(key)
            count = self.motifDict[key]

            # Discards motifs if expected value could not be calculated
            if exp == 'None': continue
            p = (self.motifDict[key] / fastLength)
#            print(str(p) + ' prob')
            mean = (fastLength * p)
#            print(str(mean) + ' mean')

            sd = math.sqrt((mean*(1-p)))
            z = (self.motifDict[key] - mean) / sd
            self.outTuple = self.outTuple + (len(key), key, count, exp, z)
            print(self.outTuple)
            
        #send to sort z then length now

    
    def seqToIndex(self, motif):
#        print(self.motif)
        if self.stop >= len(motif):
            return 0
#        self.stop += 1
#        print(self.stop)
        symbol = motif[self.stop:len(self.motif)].rstrip()
        if len(symbol) > 1: symbol = 4^(len(symbol)-1) + self.indexDict[symbol]
        print(symbol)
        print(self.indexDict[symbol])
        prefix = motif[:self.stop]

        print(prefix)
        self.stop += 1
        return 4*self.seqToIndex(prefix) + self.indexDict[symbol]
        
    def outToFile(self, outTuples):
        pass
        #don't forget to write if in range of min-max+1
    

def main():
    """ Retrieve command line arguments from class CommandLine and passes them to the StoreCounts method """
    motifDict={}
    fastLength=0
    arguments = CommandLine()
    readFast = FastAreader(arguments.args.infile)
    
    seqStr = readFast.readFasta()
    for seqs in seqStr:
        tempLength = fastLength
        getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
        
        fastLength, motifDict = getCounts.slideSeq()
        fastLength = fastLength + tempLength
        print(str(fastLength) + 'fastaLength')
    getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
    
    doCalcs = getCounts.computeSeq(fastLength, motifDict)
    
#    passToFastA = FastAreader()
#    passToFastA.FastAreader(arguments)
#    arguments.slideSeq(arguments.args)

if __name__== "__main__":
    main()
