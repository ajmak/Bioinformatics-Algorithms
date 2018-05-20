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
        self.parser.add_argument('--minMotif', nargs='?', default = 3, type = int, help = 'Set minimum motif length, default=3')
        self.parser.add_argument('--maxMotif', nargs='?', default = 8, type = int, help = 'Set maximum motif length, default=8')
        #Change this later!!!
        self.parser.add_argument('--cutoff', nargs='?', default = 0, type = int, help = 'Set minimum z cutoff score, default=0')

        self.args = self.parser.parse_args()



class StoreCounts:
    """Counts all motifs for word lengths between --minMotif and --maxMotif in sequence provided from main and returns observed motifs as a dictionary of counts"""
    def __init__(self, arguments, seq, motifDict, fastLength):
        """Instantiates necessary arguments, dictionaries, lists, and objects needed for StoreCounts class """
        self.arguments = arguments
        self.seq = seq
        self.motifDict = motifDict
        self.mini=self.arguments.minMotif
        self.maxi=self.arguments.maxMotif
        self.cutoff=self.arguments.cutoff
        self.midmin = self.mini - 2
        
        self.outTuple = ()
        self.tupleList=[]
        self.intDict= {0:'A', 1:'C', 2:'G', 3:'T'}
        self.stop = 0
        

    def slideSeq(self):
        
        allowed=set('ACGT')
        fastLength = len(self.seq)
        print(fastLength)
        j=0
        
        for k in range(self.midmin, self.maxi+1):
            for i in range(fastLength-k+1):
                self.motif = self.seq[i:i+k]
                if set(self.motif) <= allowed:
                    if self.motif not in self.motifDict.keys():
                        self.motifDict[self.motif] = 1 
                    else: self.motifDict[self.motif] += 1
                else:
                    print(self.motif)
                    pass
            j+=1

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
#        print(str(nextSwitch) + ' nextSwitch')
        for i in range(minIndex, maxIndex):
#            print(i)
            if i == nextSwitch: 
                word+=1
                
#                print(str(i) + ' switch')
                prevSwitch = nextSwitch
#                print(str(prevSwitch) + ' prevSwitch')
                nextSwitch = nextSwitch + 4**word #4**word + 4**(word -1) 
                #prevSwitch = 4**(word - 1) + 4**(word -2)
            seqInt = i - prevSwitch
            motifAtIndex = self.indexToSeq(seqInt, word)
            self.outTuple = self.zCalc(fastLengths, motifAtIndex)
            self.tupleList.append(self.outTuple)
        return self.tupleList

        
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
        if motif[1:-1] not in self.motifDict.keys() or left=='None' or right=='None':
            expected = 'None'
            return
        else: 
            mid=motif[1:-1]
            expected = (self.motifDict[left] * self.motifDict[right]) / self.motifDict[mid]
        return expected

    def zCalc(self, fastLength, motifKey):
        if motifKey in self.motifDict.keys():
#            print(self.motifDict.keys())
            exp = self.expectedVal(motifKey)
            count = self.motifDict[motifKey]

            # Discards motifs if expected value could not be calculated
            if exp == 'None': return
            else:
                mean = exp
                p = (mean / fastLength)
#                print(str(p) + ' prob' + motifKey)
                
                #(float(fastLength) * float(p))

#                print(str(mean) + ' mean')

                sd = math.sqrt((mean*(1-p)))
#                print(str(sd) + ' sd')
                z = ((count - mean) / sd)
#                print(z)
        else:
            exp = self.expectedVal(motifKey)
            if exp=='None': return
            count = 0
            z = 0
        self.outTuple = (motifKey, count, exp, z) 
        return self.outTuple

 

        
    def outToFile(self, outTuples):
        #https://stackoverflow.com/questions/19643099/sorting-a-list-of-tuples-with-multiple-conditions
        
        sortedTuple = sorted(outTuples, key=lambda x: ( (len(x[0])), -x[3] ), reverse = True)
#        print(sortedTuple)
#        print(outTuples)
        for item in sortedTuple:
#            print(item[3])
            if type(item[2]) != float or item[3] >= self.cutoff: continue

            print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(item[0], item[1], item[2], item[3]))
        #don't forget to write if in range of min-max+1
    

def main():
    """ Retrieve command line arguments from class CommandLine and passes them to the StoreCounts method """
    motifDict={}
    fastLength=0
    arguments = CommandLine()
    readFast = FastAreader()
    
    seqStr = readFast.readFasta()
    for seqs in seqStr:
        tempLength = fastLength
        getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
        fastLength, motifDict = getCounts.slideSeq()
        fastLength = fastLength + tempLength

    getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
    
    doCalcs = getCounts.computeSeq(fastLength, motifDict)
    sortAndPrint = getCounts.outToFile(doCalcs)

if __name__== "__main__":
    main()
