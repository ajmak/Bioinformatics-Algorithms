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
        self.minimizerMap = {}
        kCount = {}

        allowed=set('ACGT')
        fastLength = len(targetString)

        for i in range(fastLength-self.w+1):
            window = targetString[i:i+self.w]
            kList = []
            sortedkList = []
            if set(window) <= allowed:
                for kmer in range(self.w - self.k+1):
                    kmer = window[kmer:kmer+self.k]
                    kList.append(kmer)
                    #add to sorted list in lexicographic order while maintaining a list ordered by its occurrence in string
                    if len(sortedkList) != 0 and kmer < sortedkList[0]:
                        sortedkList.insert(0,kmer)
                    else:
                        sortedkList.append(kmer)

                minmer = sortedkList[0]
                #add to minimizerMap 
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

        #remove minmers with more than t occurences 
        for minCount in kCount.keys():
            if kCount[minCount] > self.t:
                kCount.pop(minCount)
        
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
        matchTuple=[]
        sortedkList = []
        #Finds minmers in search string and finds matches to minmers in target string
        for i in range(fastLength-self.w+1):
            window = searchString[i:i+self.w]
            kList = []
            #if all nucleotides in the window are A,C,G, or T 
            if set(window) <= allowed:
                for kmer in range(self.w - self.k+1):
                    kmer = window[kmer:kmer+self.k]
                    kList.append(kmer)
                    #add to sorted list in lexicographic order while maintaining a list ordered by its occurrence in string
                    if len(sortedkList) != 0 and kmer < sortedkList[0]:
                        sortedkList.insert(0,kmer)
                    else:
                        sortedkList.append(kmer)

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
        components = []
        edgeDict = {}
        #for all x's
        for i in range(len(seeds)):
            print('seeds ' + str(seeds[i]))
            x1, y1s = seeds[i]
            #fix later (should only be comparing y's after cur y)
            #for al y's of that x
            print(i)
            for y1 in y1s:
                print('y ' + str(y1))
                if len(components) == 0:
                    print([(x1,y1)])
                    components.append([(x1,y1)])
                    continue
                added = False
                for comp in range(len(components)):
                    for tup in components[comp]:
                        print('tup')
                        print(tup)
                        print('comp' + str(components))
                        if tup == (x1,y1): added = True
                        if (abs(tup[1] - y1) <= l) and tup != (x1,y1):
                            print('adding ' + str((x1,y1)))
                            components[comp].append((x1,y1))
                
                            added = True
                            break
                if added == False: components.append([(x1,y1)])
            #when x increments check if ends of each component are within l and merge if so
            
            a=b=0
            if i != 0:
                for comp1 in components:
                    b = 0
                    for comp2 in components:
                        if comp1 == comp2:
                            b+=1
                            continue
                        if abs(comp1[-1][0] - comp2[0][0]) <= l and abs(comp1[-1][1] - comp2[-1][0]) <= l:
                            print('joining ' + str(comp1) + str(comp2))
                            comp1 = comp1 + comp2
                            print comp1
                            components[a] = comp1
                            _ = components.pop(b)
                        b+=1
                    a+=1
                
            
        print('my components')
        print(components)
        print()
        return(components)
        #this is the ugliest code I've ever written

        """ Takes the dictionary of available edges and the defined starting node and uses a stack to store traversed nodes while removing e
dges from the dictionary of available edges. A node is appended to a connected component when it has no outgoing edges, A new connected component is created when the stack is empty"""
        """
        allComponents = [[]]
        stack = []
        curNode = edgeDict.keys()[0]
        prevNode = ''
        # while stack is not empty or dictionary of edges is not empty
        while len(edgeDict.keys()) != 0:
            if edgeDict[curNode]:
                for child in edgeDict.curNode:
                    stack.append(child)
            if len(stack) == 0: allComponents.append([])
            allComponents.append[-1](curNode)
            if not edgeDict[curNode]:
                del edgeDict[curNode]
            else:
                
                                                   



            print(curNode)
            print
            print('stack')
            print(stack)
#            if len(stack) == 0:
#                allComponents.append([])
            if curNode in edgeDict.keys() and not edgeDict[curNode]:
                print('appending ' + str(curNode))
                allComponents[-1].append(curNode)
                del edgeDict[curNode]
                curNode = stack.pop()
#                if not edgeDict[curNode]: break
                continue
            stack.append(curNode)
            prevNode = curNode
            print(allComponents)
            if curNode not in edgeDict.keys():
                curNode = edgeDict.keys()[0]
                allComponents.append([])
                continue
            curNode = edgeDict[curNode][0]
            print('deleting ' + str(edgeDict[prevNode][0]))
            del edgeDict[prevNode][0]
            print(edgeDict)
            
        allComponents[-1].append(curNode)
        print('these are the components')
        print(allComponents)
        return allComponents
        """
        

class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.
        
        Implements the Smith-Waterman algorithm: 
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        
        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to 
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?
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
        currentIndex = self.maxIndex
        parentIndex = ((int(self.traceBack_i[currentIndex[0],currentIndex[1]])), int(self.traceBack_j[currentIndex[0],currentIndex[1]]))
        self.matchIndexes.append((currentIndex[0], currentIndex[1]))

        while currentIndex[0] != 0 and currentIndex[1] != 0:
            curMatch = int(self.editMatrix[parentIndex[0], parentIndex[1]])
            if self.editMatrix[currentIndex[0],currentIndex[1]] == ((self.editMatrix[parentIndex[0],parentIndex[1]]) + 3):
                self.matchIndexes.append((currentIndex[0]-1, currentIndex[1]-1))
            # increment and loop
            currentIndex = parentIndex
            parentIndex = ( int(self.traceBack_i[currentIndex[0], currentIndex[1]]), int(self.traceBack_j[currentIndex[0], currentIndex[1]]) )
        revMatches = reversed(self.matchIndexes)

        return(list(revMatches))
    def maxScore(self):
        """ Returns the maximum alignment score
        """
        return(self.maxScore)
    
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
            if bestAlignment[0] == None or alignment.maxScore() > bestAlignment[0].maxScore():
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
            print queryIndex
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.maxScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" % 
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore)) 
            # Comment this out to test on a subset
            #if queryIndex > 100:
            #    break
    
    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" % 
                    (time.time()-startTime, float(sum(alignmentScores))/len(alignmentScores)))
    
if __name__ == '__main__':
    main()
