import array
import sys
import numpy as np
import pysam
import argparse
import logging 
import time
logger = logging.getLogger()

"""See the comments below to see the code you need to complete.
"""

class MinimizerIndexer(object):
    """ Simple minimizer based substring-indexer. 
    
    Please read: https://doi.org/10.1093/bioinformatics/bth408
    
    Related to idea of min-hash index and other "sketch" methods.
    """
    def __init__(self, targetString, w, k, t):
        """ The target string is a string/array of form "[ACGT]*".
        
        Stores the lexicographically smallest k-mer in each window of length w, such that w >= k positions. This
        smallest k-mer is termed a minmer.
        
        If a minmer occurs more than t times then it is omitted from the index.
        """
        
        self.targetString = targetString
        self.w = int(w)
        self.k = int(k)
        self.t = int(t)
        # If a minmer occurs more than t times then its entry is removed from the index
        # This is a heuristic to remove repetitive minmers that would create many spurious alignments between
        # repeats
        
        # Hash of minmers to query locations, stored as a map whose keys
        # are minmers and whose values are lists of the start indexes of
        # occurrences of the corresponding minmer in the targetString, 
        # sorted in ascending order of index in the targetString.
        #
        # For example if k = 2 and w = 4 and targetString = "GATTACATTT"
        #
        # GATTACATTT
        # GATT (AT)
        #  ATTA (AT)
        #   TTAC (AC)
        #    TACA (AC)
        #     ACAT (AC)
        #      CATT (AT)
        #       ATTT (AT)
        #
        # then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        self.minimizerMap = {}
        kCount = {}
        allowed=set('ACGT')
        fastLength = len(targetString)

#        print(self.w)
        for i in range(fastLength-self.w+1):
            window = targetString[i:i+self.w]
            kList = []
            if set(window) <= allowed:
                for kmer in range(self.w - self.k+1):
                    kmer = window[kmer:kmer+self.k]
                    #add to list in lexicographic order
                    kList.append(kmer)
                #change later to if < kList[0] then insert
                sortedkList = np.sort(kList)
                minmer = sortedkList[0]
                if minmer not in self.minimizerMap.keys():
                    self.minimizerMap[minmer] = (i+kList.index(minmer),)
                elif i+kList.index(minmer) not in self.minimizerMap[minmer]:
                    self.minimizerMap[minmer] += (i+kList.index(minmer),)
                                                                
                #add to counts dict
                if minmer not in kCount.keys():
                    kCount[minmer] = 1
                else:
                    kCount[minmer] +=1
            else:
                continue
        print(self.minimizerMap)
        #remove minmers with more than t occurences 
        for minCount in kCount.keys():
                kCount.pop(minCount)
        for tup in self.getMatches('GATTTA'):
            print(tup)
        
    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index searchString and y is an occurrence in targetString.
        
        For example if k = 2 and w = 4 and targetString = "GATTACATTT" and searchString = "GATTTA"
        then self.minimizerMap = { "AT":(1,6), "TA":(3,), "AC":(4,), "CA":(5,), }
        and getMatches will yield the following sequence:
        (1, (1,6)), (4, (3,))
        
        You will need to use the "yield" keyword
        """
        allowed=set('ACGT')
        fastLength = len(searchString)
        matchTuple=()
        print(searchString)
        for i in range(fastLength-self.w+1):
            window = searchString[i:i+self.w]
            kList = []
            if set(window) <= allowed:
                for kmer in range(self.w - self.k+1):
                    kmer = window[kmer:kmer+self.k]
                    #add to list in lexicographic order
                    kList.append(kmer)
                    #change later to if < kList[0] then insert
                sortedkList = np.sort(kList)
                minmer = sortedkList[0]


                if minmer in self.minimizerMap.keys():
                    if i+kList.index(minmer) not in matchTuple:
                       
                        matchTuple += (i+kList.index(minmer), self.minimizerMap[minmer],)
            else:
                continue

        yield(matchTuple)


class SeedCluster:
    """ Represents a set of seeds between two strings.
    """
    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed 
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = map(lambda(x, y) : y, seeds)
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.
        
        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of shared k-mer in both strings, termed a *seed*, such that the k-mer 
        occurrence starts at position x in the first string and starts at position y_i in the second string.
        
        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.   
        
        Consider a graph in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).
        
        The return value is a Python set of components, each component is a SeedCluster object.
        
        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use.  
        """
        edgeDict = {}
        for i in range(len(seeds)):
            x1, y1s = seeds[i]

            for y1 in y1s:
                for y1b in y1s:
                    if abs(y1-y1b) <=l:
                        if (x1,y1) not in edgeDict.keys():
                            edgeDict[(x1,y1)] =[(x1,y1b)]
                        else:
                            edgeDict[(x1,y1)].append((x1,y1b))
                            
            for j in range(i+1, len(seeds)):
                x2, y2s = seeds[j]
                if (abs(x1-x2) <= l):
                    for y1 in y1s:
                        for y2 in y2s:
                            print(x1,y1)
                            print(x2,y2)
                            if (abs(y1-y2)) <= l:
                                if (x1,y1) not in edgeDict.keys():
                                    edgeDict[(x1,y1)] = [(x2,y2),(x1,y1)]
                                else:
                                    edgeDict[(x1,y1)].append((x2,y2))
                                if (x2,y2) not in edgeDict.keys():
                                    edgeDict[(x2,y2)] = [(x1,y1),(x2,y2)]
                                else:
                                    edgeDict[(x2,y2)].append((x1,y1))
        print('my edgeDict')
        print(edgeDict)
        print()
        #this is the ugliest code I've ever written

        
        finalPath = []
        # while stack is not empty or dictionary of edges is not empty
        while len(self.stack) != 0 or len(edgeDict[curNode]) != 0:
            if curNode not in edgeDict.keys() or not edgeDict[curNode]:
                finalPath.append(curNode)
                curNode = self.stack.pop()
                continue
            self.stack.append(curNode)
            prevNode = curNode
            curNode = edgeDict[curNode][0]
            del edgeDict[prevNode][0]
        finalPath.append(curNode)
        return finalPath
        
        

class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.
        
        Implements the Smith-Waterman algorithm: 
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        
        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to 
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?
        """
        # Code to complete to compute the edit matrix
        self.editMatrix = np.zeros(shape=(len(string1)+1, len(string2)+1))
        self.traceBack_i = np.zeros(shape=(len(string1)+1, len(string2)+1))
        self.traceBack_j = np.zeros(shape=(len(string1)+1, len(string2)+1))
        print(string1)
        print(string2)
        self.maxScore = 0

        for i in range(1,len(string1)+1):
            for j in range(1,len(string2)+1):
                parent1 = self.editMatrix[i-1,j] - 2
                parent2 = self.editMatrix[i, j-1] - 2
                if string1[i-1] == string2[j-1]:
#                    parent1 = self.editMatrix[i-1,j] 
#                    parent2 = self.editMatrix[i, j-1] 
                    parent3 = self.editMatrix[i-1, j-1] + 3
                else:
#                    parent1 = self.editMatrix[i-1,j] - 2
#                    parent2 = self.editMatrix[i, j-1] - 2
                    parent3 = self.editMatrix[i-1, j-1] - 3
                if parent1 > parent2:
                    bestParent = parent1
                    bestIndex = (int(i-1),j)
                # if parent2 is greater or equal to parent1
                else:
                    bestParent = parent2
                    bestIndex = (i,int(j-1))
                # choose diagonal if parent 3 greater or equal to bestParent
                if bestParent <= parent3:
                    bestParent = parent3
                    bestIndex = (int(i-1),int(j-1))
                if bestParent < 0: bestParent = 0
                self.editMatrix[i,j] = bestParent
                self.traceBack_i[i,j] = int(bestIndex[0])
                self.traceBack_j[i,j] = int(bestIndex[1])
                if bestParent > self.maxScore:
                    self.maxScore = bestParent
                    self.maxIndex = bestIndex
                    
        print(self.editMatrix)
        print self.traceBack_i
        print self.traceBack_j        

    def getAlignment(self):
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.
        
        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ] 
        """
        self.matchIndexes = []
        currentIndex = self.maxIndex
        print "currentIndex:" + str(currentIndex)
        print "[0]:" + str(currentIndex[0])
        print "[1]:" + str(currentIndex[1])
        parentIndex = ((int(self.traceBack_i[currentIndex[0],currentIndex[1]])), int(self.traceBack_j[currentIndex[0],currentIndex[1]]))
        matchList = []

        print parentIndex
        print currentIndex
        while currentIndex[0] != 0 and currentIndex[1] != 0:
            curMatch = int(self.editMatrix[parentIndex[0], parentIndex[1]])
            if self.editMatrix[currentIndex[0],currentIndex[1]] == ((self.editMatrix[parentIndex[0],parentIndex[1]]) + 3):
                self.matchIndexes.append((currentIndex[1]-1, currentIndex[0]-1))
            # increment and loop
            currentIndex = parentIndex
            print('traceback')
            print(self.traceBack_i[currentIndex[0], currentIndex[1]])
            parentIndex = ( int(self.traceBack_i[currentIndex[0], currentIndex[1]]), int(self.traceBack_j[currentIndex[0], currentIndex[1]]) )
        print(self.maxScore)
        print(self.matchIndexes)
        print(reversed(self.matchIndexes))
    
    def maxScore(self):
        """ Returns the maximum alignment score
        """
        return self.maxScore
    
def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.
    
    Maps the string in both its forward and reverse complement orientations.
    
    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it, 
    and note any drawbacks this might have.
    """
    bestAlignment = [None]
    
    def mapForwards(queryString):
        """ Maps the query string forwards
        """
        # Find seed matches, aka "aligned kmers"
        seeds = list(minimizerIndex.getMatches(queryString))
        
        # For each cluster of seeds
        for seedCluster in SeedCluster.clusterSeeds(list(seeds), l=config.l):
            
            # Get substring of query and target to align
            queryStringStart = max(0, seedCluster.minX - config.c) # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.c) # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]
            
            targetStringStart = max(0, seedCluster.minY - config.c) # Inclusive coordinate
            targetStringEnd = min(len(targetString), seedCluster.maxY + config.k + config.c) # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]
            
            #print "target_aligning", targetStringStart, targetStringEnd, targetSubstring
            #print "query_aligning", queryStringStart, queryStringEnd, querySubstring
            
            # Align the genome and read substring
            alignment = SmithWaterman(targetSubstring, querySubstring, 
                                      gapScore=config.gapScore, 
                                      matchScore=config.matchScore,
                                      mismatchScore=config.mismatchScore)
            
            # Update best alignment if needed
            if bestAlignment[0] == None or alignment.getMaxAlignmentScore() > bestAlignment[0].getMaxAlignmentScore():
                bestAlignment[0] = alignment
        
        return bestAlignment
    
    def reverseComplement(string):
        """Computes the reverse complement of a string
        """
        rMap = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
        return "".join(rMap[i] for i in string[::-1])
                
    # Run mapping forwards and reverse
    mapForwards(queryString)
    mapForwards(reverseComplement(queryString))
    
    return bestAlignment[0]

        
def main():
    w = 4
    k = 2
    t = 1000
    queryString = "TGTTACGG"
    targetString = "GGTTGACTA"
#    queryString = "CGTTGCCATCTCCGTATGA"
#    targetString = "ATTCTGC"
#    targetString = "GATTACATTT"
#    queryString = "GATTTA"


    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w, k, t)

    matches = minimizerIndex.getMatches(queryString)
    alignmentScores = [] # Array storing the alignment scores found
    #class SmithWaterman(object):def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3)
    alignment = SmithWaterman(targetString, queryString)
#    getEditMatrix = 
    alignment.getAlignment()
if __name__ == '__main__':
    main()
