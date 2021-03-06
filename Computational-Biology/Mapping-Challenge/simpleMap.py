"""
Allysia Mak
CURRENT VERSION
"""


import array
import sys
import numpy as np
import pysam
import argparse
import logging 
import time
from collections import defaultdict
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
        self.w = w
        self.k = k
        self.t = t # If a minmer occurs more than t times then its entry is removed from the index
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
        self.minimizerMap = self.getMinimizerMap()
        
    def getMinimizerMap(self):
        self.minimizerMap = defaultdict(set)
        kCount = defaultdict(int)
        allowed=set('ACGT')
        fastLength = len(self.targetString)

        for i in range(fastLength-self.w+1):
 
            window = self.targetString[i:i+self.w]
            kList = []
            minMinmer = None
            minDex = None
            if set(window) <= allowed:
                #for kmer in window
                for j in range(self.w - self.k+1):
                    kmer = window[j:j+self.k]
                    #keep track of minmer
                    if minMinmer == None or kmer < minMinmer:
                        minMinmer = kmer
                        minDex = i+j

                #add to minimizerMap 
                kCount[minMinmer] += 1
                if kCount[minMinmer] <= self.t:
                    self.minimizerMap[minMinmer].add(minDex)
                    

            else:
                continue


        return self.minimizerMap
    
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
        #rewrote minimizer code so could compare to target while getting minmers in searchString 
        fastLength = len(searchString)
        matchTuple=[]
        sortedkList = []
        #Finds minmers in search string and finds matches to minmers in target string
        for i in range(fastLength-self.w+1):
            window = searchString[i:i+self.w]
            minMinmer = None
            minDex = None

            for j in range(self.w - self.k+1):
                kmer = window[j:j+self.k]
                if minMinmer == None or kmer < minMinmer:
                        minMinmer = kmer
                        minDex = i+j
            #compare minmer in target to minmers in query
            if minMinmer in self.minimizerMap:
                if minDex not in matchTuple:

                    yield(minDex, self.minimizerMap[minMinmer])



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
        ys = map(lambda (x, y): y, seeds)
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
        (ANSWER 1): You could build the graph of connected seeds completely, then build connected components by running depth or breadth first search starting with random nodes. A DFS or BFS will touch each node of a connected component of a graph, so you just need to pick a random node that is not yet included in a connected component as the starting node for the search, and the search will build a new connected component.
        My method for seed clustering compares each seed to seeds already added into components. This is able to take care of all seed clusters in the x range and most of the clusters in the y range since the seeds are given in order. There is a case that a seed with the same x but much higher y should be in that cluster but was not added. To take care of this, every time the x value is incremented the first value of each list is compared to the last value of every list preceding it. If the x's and y's are in range then the two components are merged.

        """
        components = []
        edgeDict = {}
        #for all x's

        if len(seeds) == 0 or seeds == [[]]: return 
        for i in range(len(seeds)):

            x1, y1s = seeds[i]
            #for al y's of that x

            for y1 in y1s:
                #add very first seed to its own component list
                if len(components) == 0:
                    components.append([(x1,y1)])
                    continue
                added = False
                #compare to all seeds already put into components
                for comp in range(len(components)):
                    for tup in components[comp]:
                        #if you found the seed in a cluster then it was already added
                        if tup == (x1,y1): added = True
                        if (abs(x1 - tup[0]) <= l) and (abs(tup[1] - y1) <= l) and tup != (x1,y1):
               
                            components[comp].append((x1,y1))
                
                            added = True
                            break
                #if seed was not added to a component list then put it in its own component list
                if added == False: components.append([(x1,y1)])
                added = False
            #when x increments check if ends of each component are within l and merge if so
            
            a=b=0
            if i != 0:
                for comp1 in components:
                    b = 0
                    for comp2 in components:
                        if comp1 == comp2:
                            b+=1
                            continue
                        if abs(comp1[-1][0] - comp2[0][0]) <= l and abs(comp1[-1][1] - comp2[0][1]) <= l:

                            comp1 = comp1 + comp2

                            components[a] = comp1
                            _ = components.pop(b)
                        b+=1
                    a+=1
                

        for component in components:
            yield SeedCluster(set(component))
        #this is the ugliest code I've ever written


    def seeds(self):
        return self.clusterSeeds()
        
                

class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.
        
        Implements the Smith-Waterman algorithm: 
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        
        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to 
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?
        (ANSWER 2): K-tuple methods (such as BLAST) improve on the speed of the alignment but compromise on accuracy due to the heuristic nature of the algorithm.

        """
        
        self.editMatrix = np.zeros(shape=(len(string1)+1, len(string2)+1))
        self.traceBack_i = np.zeros(shape=(len(string1)+1, len(string2)+1))
        self.traceBack_j = np.zeros(shape=(len(string1)+1, len(string2)+1))
        self.maxScore = 0
        
        for i in range(1,len(string1)+1):
            for j in range(1,len(string2)+1):
                parent1 = self.editMatrix[i-1,j] + gapScore
                parent2 = self.editMatrix[i, j-1] + gapScore
                                    
                if string1[i-1] == string2[j-1]:
                    parent3 = self.editMatrix[i-1, j-1] + matchScore
                else:
                    parent3 = self.editMatrix[i-1, j-1] + mismatchScore
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
                #store maxScore and maxIndex
                if bestParent > self.maxScore:
                    self.maxScore = bestParent
                    self.maxIndex = bestIndex

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
        #start at the index of where the maxScore occurred
        currentIndex = self.maxIndex
        #get backpointer to best parent that produced currentIndex
        parentIndex = ((int(self.traceBack_i[currentIndex[0],currentIndex[1]])), int(self.traceBack_j[currentIndex[0],currentIndex[1]]))
        self.matchIndexes.append((currentIndex[0], currentIndex[1]))
        
        while currentIndex[0] != 0 and currentIndex[1] != 0:
            curMatch = int(self.editMatrix[parentIndex[0], parentIndex[1]])
            #if the currentIndex + 3 to the best parent indicates it was a match then append the currentIndex to the final path
            if self.editMatrix[currentIndex[0],currentIndex[1]] == ((self.editMatrix[parentIndex[0],parentIndex[1]]) + 3):
                self.matchIndexes.append((currentIndex[0]-1, currentIndex[1]-1))
            # increment and loop
            currentIndex = parentIndex
            parentIndex = ( int(self.traceBack_i[currentIndex[0], currentIndex[1]]), int(self.traceBack_j[currentIndex[0], currentIndex[1]]) )
        #reverse indices list
        revMatches = reversed(self.matchIndexes)

        return(list(revMatches))
    def maxScore(self):
        """ Returns the maximum alignment score
        """
        return(maxScore)
    
def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.
    
    Maps the string in both its forward and reverse complement orientations.
    
    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it, 
    and note any drawbacks this might have.
    (ANSWER 3): Filtering out non-(ACGT) sequences would likely help speed up the downstream methods. seeds is cast into a list twice which is confusing. If you want this code to be faster you should probably not use Python. Drawbacks of this are maintainability.
    
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
            if bestAlignment[0] == None or alignment.maxScore > bestAlignment[0].maxScore:
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

class Config():
    """ Minimal configuration class for handing around parameters
    """
    def __init__(self):
        self.w = 30
        self.k = 20
        self.t = 10
        self.l = 30
        self.c = 100
        self.gapScore=-2
        self.matchScore=3
        self.mismatchScore=-3
        self.logLevel = "INFO"
        
def main():
    # Read parameters
    config = Config()
    
    #Parse the inputs args/options
    parser = argparse.ArgumentParser(usage="target_fasta query_fastq [options]", version="%prog 0.1")

    parser.add_argument("target_fasta", type=str,
                        help="The target genome fasta file.")
    parser.add_argument("query_fastq", type=str,
                        help="The query sequences.")
    
    parser.add_argument("--w", dest="w", help="Length of minimizer window. Default=%s" % config.w, default=config.w)
    parser.add_argument("--k", dest="k", help="Length of k-mer. Default=%s" % config.k, default=config.k)
    parser.add_argument("--t", dest="t", help="Discard minmers that occur more frequently " 
                                            "in the target than t. Default=%s" % config.w, default=config.w)
    parser.add_argument("--l", dest="l", help="Cluster two minmers into the same cluster if within l bases of"
                                            " each other in both target and query. Default=%s" % config.l, default=config.l)
    parser.add_argument("--c", dest="c", help="Add this many bases to the prefix and suffix of a seed cluster in the"
                                            " target and query sequence. Default=%s" % config.c, default=config.c)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" % 
                      config.mismatchScore, default=config.mismatchScore)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % 
                      config.logLevel, default=config.logLevel)
    
    options = parser.parse_args()
    
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")
    
    startTime = time.time()
    
    # Parse the target sequence and read the first sequence
    with pysam.FastaFile(options.target_fasta) as targetFasta:
        targetString = targetFasta.fetch(targetFasta.references[0])
    logger.info("Parsed target string. Length: %s" % len(targetString))
    
    #TAKE THIS OUT LATER
#    targetString = targetString[80000:80500]

    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w=options.w, k=options.k, t=options.t)

    minmerInstances = sum(map(len, minimizerIndex.minimizerMap.values()))
    logger.info("Built minimizer index in %s seconds. #minmers: %s, #minmer instances: %s" %
                 ((time.time()-startTime), len(minimizerIndex.minimizerMap), minmerInstances))
    
    # Open the query files
    alignmentScores = [] # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, xrange(sys.maxint)):

            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.maxScore
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" % 
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore)) 
            # Comment this out to test on a subset
#            if queryIndex > 10:
#                break
    
    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" % 
                    (time.time()-startTime, float(sum(alignmentScores))/len(alignmentScores)))
    
if __name__ == '__main__':
    main()
