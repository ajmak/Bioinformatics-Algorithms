#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: Lori Siao, Xian Chang, Ioannis Anastopulos, Alexis Thornton, Balaji Sundararaman, Vinay Poodari

"""

Provides comparison of expected counts for all possible motifs of given length (default 3 to 8) to actual counts in genome. This information can be used to help determine underrepresented sequences in the given genome.
Input: FASTA file for a single genome 
Output: All possible motifs with corresponding counts from the genome, the theoretical expected counts for that motif, and the Z-scores for comparison with the null model
"""
import sys
import argparse
import math
import string
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
        # CHANGED: ranges are now correct 
        self.parser.add_argument('--minMotif', nargs='?', default = 3, type = int, choices = range(3,8), help = 'Set minimum motif length, default=3')
        self.parser.add_argument('--maxMotif', nargs='?', default = 8, type = int, choices = range(4,9), help = 'Set maximum motif length, default=8')
        self.parser.add_argument('--cutoff', nargs='?', default = 0, type = int, help = 'Set minimum z cutoff score, default=0')

        self.args = self.parser.parse_args()



class StoreCounts:
    """Counts all motifs for word lengths between --minMotif and --maxMotif in sequence provided from main and returns observed motifs as a dictionary of counts"""
    def __init__(self, mini, maxi, cutoff):
        """Instantiates necessary arguments, dictionaries, lists, and objects needed for StoreCounts class """


        self.mini = mini
        self.maxi= maxi
        self.cutoff= cutoff
        self.midmin = self.mini - 2
        
        self.intDict= {0:'A', 1:'C', 2:'G', 3:'T'}
        self.stop = 0
        

    def slideSeq(self, seq, motifDict):
        """Finds and stores motifs and counts of motifs for word lengths between --minMotif and --maxMotif by sliding a window of size word across sequence """
        allowed=set('ACGT')
        fastLength = len(seq)
        
        for k in range(self.midmin, self.maxi+1):
            for i in range(fastLength-k+1):
                self.motif = seq[i:i+k]
                # https://stackoverflow.com/questions/20726010/how-to-check-if-a-string-contains-only-characters-from-a-given-set-in-python
                if set(self.motif) <= allowed:
                    
                    if self.motif not in motifDict.keys():
                        motifDict[self.motif] = 1 
                    else: motifDict[self.motif] += 1
                else:
                    continue

        return fastLength, motifDict
    
    def computeSeq(self, fastLengths, motifDict):
        """
        Utilizes a integer to string formula for letters A, C, G, T to determine all possible motifs for word lengths from --minMotif to --maxMotif

        For wordlength k, the formula
        4^k + 4^k-1 + ... + 4^1 
        stores the length of the word (k)

        The formula
        4*indexToSeq(prefix) + a
        where a is the value of the last letter of the motif as defined by A:0, C:1, G:2, T:3
        stores each nucleotide of the motif
        """
        tupleList = []
        nextSwitch=0
        ###########################################################
        # By determining the minimum and maximum integers as a    #
        # function of k length, we find an integer representation #
        # of all possible motifs of lengths between --minMotif and# 
        # --maxMotif                                              #
        maxIndex = minIndex = 0
        for i in range(1, self.maxi): maxIndex += 4**i
        for i in range(1, self.mini-1): minIndex += 4**i
#        if self.mini == 4:
#            minIndex = 4**(self.mini - 1) + 4**(self.mini -2) + 4
#            maxIndex = 4**self.maxi + 4**(self.maxi -1) + 4**self.maxi + 4*self.maxi 
#        else:
#            minIndex = 4**(self.mini - 1) + 4**(self.mini -2) 
#            maxIndex = 4**self.maxi + 4**(self.maxi -1) + 4*self.maxi + 4
        word = self.mini
        # Each 'Switch' defines a point where word length increases by 1
        prevSwitch = minIndex
        nextSwitch = prevSwitch
#        for i in range(1,self.mini + 1):
#            prevSwitch = nextSwitch
#            nextSwitch += 4**i
        for i in range(minIndex, maxIndex):
            if i == nextSwitch: 
                word+=1
                prevSwitch = nextSwitch
                nextSwitch = nextSwitch + 4**word 
            seqInt = i - prevSwitch
            print(str(i) + ' i')
            print(str(prevSwitch) + ' prevswitch')
            print(str(nextSwitch) + ' nextswitch')
            print(str(seqInt) + ' seqInt')
            print(str(word) + ' minimum')
            print(self.maxi)
            motifAtIndex = self.indexToSeq(seqInt, word)
            print(motifAtIndex)
            # performs z-score calculations
            outTuple = self.zCalc(fastLengths, motifAtIndex, motifDict)
            # CHANGED only appends tuple if is not None
            if outTuple != None: tupleList.append(outTuple)
        return tupleList

        
    def indexToSeq(self, seqInt, k):
        """Works with computeSeq by processing the portion of the integer that stores the sequence of the motif and returns the sequence """
        if k==1:
            #returns first letter of string
#            print(self.intDict[seqInt])
            return self.intDict[seqInt]
        
        prefixIndex = seqInt//4
        remainder = seqInt%4
        return (self.indexToSeq(prefixIndex, k-1)) + self.intDict[remainder]


    def expectedVal(self, motif, motifDict):
        """Computes the expected value for a given motif """
        if motif[:-1] not in motifDict.keys(): left='None'
        else: left=motif[:-1]
        if motif[1:] not in motifDict.keys(): right='None'
        else: right=motif[1:]
        if motif[1:-1] not in motifDict.keys() or left=='None' or right=='None':
            #used as flag so that motifs that have uncomputable expected values are not printed
            expected = None
            return
        else: 
            mid=motif[1:-1]
            expected = (motifDict[left] * motifDict[right]) / motifDict[mid]
        return expected

    def zCalc(self, fastLength, motifKey, motifDict):
        """Calculates Z-score from expected value, mean, standard deviation, and genome length and returns the motif, observed count of motif, expected count of motif, and z-score """
        outTuple = ()
        if motifKey in motifDict.keys():
#            print(motifDict.keys())
            exp = self.expectedVal(motifKey, motifDict)
            count = motifDict[motifKey]

            # Discards motifs if expected value could not be calculated
            if exp == None: return
            else:
                mean = exp
                p = (mean / fastLength)
                sd = math.sqrt((mean*(1-p)))

                z = ((count - mean) / sd)
                
            

        else:
            exp = self.expectedVal(motifKey, motifDict)
            if exp == None: return
            count = 0
            z = 0
        outTuple = (motifKey, count, exp, z) 
        return outTuple

        
    def outToFile(self, outTuples):
        """ Prints Motif, Actual Count, Expected Count, and Z-score to stdout """
        #https://stackoverflow.com/questions/19643099/sorting-a-list-of-tuples-with-multiple-conditions
        
        sortedTuple = sorted(outTuples, key=lambda x: ( (len(x[0])), -x[3] ), reverse = True)
        print('{0}\t{1}\t{2}\t{3}'.format('Motif', 'Actual Count', 'Expected Count', 'Z-score'))
        for item in sortedTuple:
            # does not print if expected value could not be calculated or if the z-score is above --cutoff
            if item[3] <= self.cutoff:
                # print statement written by David Bernick
                print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(item[0], item[1], item[2], item[3]))
      

def main():

    motifDict={}
    fastLength=0
    arguments = CommandLine()
    # CHANGED pass all arguments parsed in main instead of parsed in __init__
    mini = arguments.args.minMotif
    maxi = arguments.args.maxMotif
    cutoff = arguments.args.cutoff
    getCounts = StoreCounts(mini, maxi, cutoff)

    readFast = FastAreader()
    seqStr = readFast.readFasta()
    for seqs in seqStr:
        tempLength = fastLength
#        getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
        fastLength, motifDict = getCounts.slideSeq(seqs, motifDict)
        fastLength = fastLength + tempLength
#        print(seqs)
#        print(motifDict)
#        print(fastLength)
#    print(len(motifDict))
#    getCounts = StoreCounts(arguments.args, seqs, motifDict, fastLength)
    
    doCalcs = getCounts.computeSeq(fastLength, motifDict)
    sortAndPrint = getCounts.outToFile(doCalcs)

if __name__== "__main__":
    main()
