#!/usr/bin/env python3
#Name: Allysia Mak (amak)
#Group Members: none

"""

Input: 
Output:
"""
import sys


class FastAreader():
    
    def __init__ (self, infile):
        '''contructor: saves attribute fname '''
        self.fname = infile
        #print(infile)    
    def doOpen (self, fname):
        if self.fname is '':
            return sys.stdin
        else:
            return self.fname
            #return open(self.fname)
 
    def readFasta(self, infile):
    
        header = ''
        sequence = ''


        with self.doOpen(self.fname) as fileH:

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
    """ Parses command line arguments, both required and optional and returns them to main"""
    def __init__(self):
        """ Initializes argument parser """


        import argparse
        #adapted from David Bernick's program.py skeleton code
        
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('--minMotif', nargs='?', default = 3, type = int, help = 'Set minimum motif length, default=3')
        self.parser.add_argument('--maxMotif', nargs='?', default = 8, type = int, help = 'Set maximum motif length, default=8')
        #Change this later!!!
        self.parser.add_argument('--cutoff', nargs='?', default = 0, type = int, help = 'Set minimum z cutoff score, default=0')
        self.parser.add_argument("infile", type=argparse.FileType('r'))
        self.parser.add_argument("outfile", type=argparse.FileType('w'))
        self.args = self.parser.parse_args()
        #http://docs.python.org/2/library/argparse.html#type


class StoreCounts:
    """ """
    def __init__(self, arguments, seq):
        self.arguments = arguments
        self.seq = seq
        self.motifDict={}
        self.mini=self.arguments.minMotif
        self.maxi=self.arguments.maxMotif
        z=self.arguments.cutoff
        self.midmin = self.mini - 2
        
        self.indexDict={'A':0,'C':1,'G':2,'T':3}
        self.stop = 0
        

    def slideSeq(self, args, seq):
        print(self.maxi)
        for each in self.seq:
            print(each)
            for k in range(self.midmin, self.maxi+1):
                for i in range(len(each)-k+1):

                    self.motif = each[i:i+k]
                    #print(k)
                    print('sending')
                    index = self.motif
                    #index = self.seqToIndex(self.motif)
                    if index not in self.motifDict.keys():
                        self.motifDict[index] = 1 
                    else: self.motifDict[index] += 1
        print(self.motifDict)
    def expectedVal(self, motif):
        if motif[:-1] not in self.motifDict.keys(): left=1
        else: left=motif[:-1]
        if motif[0:] not in self.motifDict.keys(): right=1
        else: right=motif[0:]
        if motif[1:-1] not in self.motifDict.keys(): mid=1
        else: left=motif[1:-1]
        expected = (self.motifDict[left] * self.motifDict[right]) / self.motifDict[mid]


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
        
    def indexToSeq(self, seqInt, k):
        if k==1:
            #returns first letter of string
            return self.intDict[seqInt]
        prefixIndex = seqInt//4
        remainder = seqInt%4
#        print((self.indexToSeq(prefixIndex, k-1)) + self.intDict[remainder])
        return (self.indexToSeq(prefixIndex, k-1)) + self.intDict[remainder]
    def outToFile(self, motifDict):
        pass
        #don't forget to write if in range of min-max+1
    

def main():
    """ Retrieve command line arguments from class CommandLine and passes them to the StoreCounts method """
    arguments = CommandLine()
    readFast = FastAreader(arguments.args.infile)

    seqStr = readFast.readFasta(arguments.args.infile)
    doCalcs = StoreCounts(arguments.args, seqStr)
    doCalcs.slideSeq(arguments.args, seqStr)
#    passToFastA = FastAreader()
#    passToFastA.FastAreader(arguments)
#    arguments.slideSeq(arguments.args)

if __name__== "__main__":
    main()
