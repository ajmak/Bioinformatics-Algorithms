import array
import sys
import numpy as np
import pysam
import argparse
import logging
import time
from collections import defaultdict
from collections import Counter

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

        If a minmer occurs in the target sequence more than t times as a minmer then it is omitted from the index, i.e. if the given minmer (kmer) is a minmer
        in more than t different locations in the target string. Note, a minmer may be the minmer for more than t distinct windows
        and not be pruned, we remove minmers only if they have more than t distinct occurrences as minmers in the sequence.
        """

        self.targetString = targetString
        self.w = w
        self.k = k
        self.t = t
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
        self.minimizerMap = self.generateMinimizerDict()

    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index in searchString and y is an occurrence in targetString.

        For example if k = 2 and w = 4 and targetString = "GATTACATTT" and searchString = "GATTTAC"
        then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        and getMatches will yield the following sequence:
        (1, (1,6)), (5, (4,))

        You will need to use the "yield" keyword
        """
        # Counter dictionary for keeping count of minmers
        minmer_count_dict = Counter()

        # Iterate through string at given window size
        for w in range(len(searchString) - self.w + 1):
            # Splice
            seq_splice = searchString[w:w + self.w]
            temp_minmer = [0]

            for k in range(len(seq_splice) - self.k + 1):
                minmer = seq_splice[k:k + self.k]
                index = w + k


                if temp_minmer[-1] == 0 or (minmer,index) < temp_minmer:
                    temp_minmer = (minmer,index)

            minmer_count_dict[temp_minmer[0]] +=1

            print("yieldVal: ("+str( temp_minmer[1])+","+str(self.minimizerMap[temp_minmer[0]])+")")
            if minmer_count_dict[temp_minmer[0]] <= self.t:
                yield temp_minmer[1], self.minimizerMap[temp_minmer[0]]

    def generateMinimizerDict(self):
        """

        :return:
        """
        minmer_count_dict = Counter()
        self.minimizerMap = defaultdict(set)
        for w in range(len(self.targetString) - self.w + 1):
            seq_splice = self.targetString[w:w + self.w]
            temp_minmer = [0]

            for k in range(len(seq_splice) - self.k + 1):
                minmer = seq_splice[k:k + self.k]
                index = w + k

                if temp_minmer[-1] == 0 or (minmer, index) < temp_minmer:
                    temp_minmer = (minmer, index)

            minmer_count_dict[temp_minmer[0]] += 1

            if minmer_count_dict[temp_minmer[0]] <= self.t:
                self.minimizerMap[temp_minmer[0]].add(temp_minmer[1])

        return self.minimizerMap


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
    def generateGraph(seeds, l):
        """

        :param seeds:
        :param l:
        :return:
        """

        # flatting out seeds to individual (x,y) values
        seeds = [(x[0], y) for x in seeds for y in x[1] if
                 len(x[1]) > 0]

        graph = {}  # building graph with immediate neighbors for each seed
        for seed1 in seeds:
            # comparing each seed to every other seed
            for seed2 in seeds:
                # checking if condition of distance is met
                if abs(seed1[0] - seed2[0]) <= l and abs(seed1[1] - seed2[1]) <= l:
                    # adding immediate neighbors to each seed
                    if seed1 in graph:
                        graph[seed1].add(seed2)
                    else:
                        graph[seed1] = {seed2}
        return graph, seeds

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.

        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of a shared k-mer in both strings, termed a *seed*, such that the k-mer
        occurrence starts at position x in the first string and starts at position y_i in the second string.
        The input seeds list contains no duplicates and is sorted in ascending order,
        first by x coordinate (so each successive tuple will have a greater
        x coordinate), and then each in tuple the y coordinates are sorted in ascending order.

        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.

        Consider a *seed graph* in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).

        The return value is a Python set of SeedCluster object, each representing a connected component of seeds in the
        seed graph.

        Acknowledgements:
        http://code.activestate.com/recipes/578507-strongly-connected-components-of-a-directed-graph/
        Yianni


        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use.

        Alternative strategies by which the seeds could be clos
        """

        # Code to complete - you are free to define other functions as you like

        graph, seeds = SeedCluster.generateGraph(seeds,l)

        visited = set()  # set of already visited seeds
        stack = []  # list of seeds that have already been visited, used to make connected components
        index = {}  # {seed: index of seeds in stack}
        boundaries = []  # list of indeces of seeds
        #
        for seed in seeds:
            if seed not in index:
                # three separate tags are used:
                # visit tag are for seed that need to be visited
                # visited tag are for seeds that have already been visited
                # neighbors tag are for the immediate neighbors of a seed from graph
                nodes = [('visit', seed)]  # (type=to be visited, seed)
                while nodes:  # when nodes is empty
                    node, seed = nodes.pop()
                    if node == 'visit':  # checking seeds labeled with  "to be visited"
                        index[seed] = len(stack)
                        stack.append(seed)
                        boundaries.append(index[seed])
                        nodes.append(('visited', seed))  # seed has been visited, so type changes
                        nodes.extend([('neighbors', neigh) for neigh in
                                      graph[seed]])  # adding immediate neighbors to nodes list with "neighbors" as type
                    elif node == 'neighbors':  # checking seeds labeled with 'neighbors'
                        if seed not in index:
                            nodes.append(('visit',
                                          seed))  # neighbors have now to be visited, so added to nodes with label "visit"
                        elif seed not in visited:
                            while index[seed] < boundaries[-1]:
                                boundaries.pop()
                    else:  # this case is if seeds have already been visited
                        if boundaries[-1] == index[seed]:
                            boundaries.pop()
                            connected_component = set(stack[index[seed]:])  # creating connected components
                            del stack[index[seed]:]
                            visited.update(connected_component)  # connected components have now been visited
                            # print(list(SeedCluster(connected_component)))
                            yield (SeedCluster(connected_component))

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

        Alternative implementations to accelerate the finding of reasonable local alignments may be to
        """

        self.edit_matrix = np.zeros((len(string2) + 1, len(string1) + 1))

        for row in range(1, len(string2) + 1):
            for column in range(1, len(string1) + 1):
                if column == 1:
                    if string2[row - 1] == string1[column - 1]:
                        self.edit_matrix[row][column] = matchScore
                    else:
                        if self.edit_matrix[row - 1][column] + gapScore > self.edit_matrix[row - 1][
                            column] + mismatchScore and self.edit_matrix[row - 1][column] + gapScore > 0:
                            self.edit_matrix[row][column] = self.edit_matrix[row - 1][column] + gapScore
                        if self.edit_matrix[row - 1][column] + mismatchScore > self.edit_matrix[row - 1][
                            column] + gapScore and self.edit_matrix[row - 1][column] + mismatchScore > 0:
                            self.edit_matrix[row][column] = self.edit_matrix[row - 1][column] + gapScore
                if row == 1:
                    if string2[row - 1] == string1[column - 1]:
                        self.edit_matrix[row][column] = matchScore
                    else:
                        if self.edit_matrix[row][column - 1] + gapScore > self.edit_matrix[row][
                                    column - 1] + mismatchScore and self.edit_matrix[row][column - 1] + gapScore > 0:
                            self.edit_matrix[row][column] = self.edit_matrix[row][column - 1] + gapScore
                        if self.edit_matrix[row][column - 1] + mismatchScore > self.edit_matrix[row][
                                    column - 1] + gapScore and self.edit_matrix[row][column - 1] + mismatchScore > 0:
                            self.edit_matrix[row][column] = self.edit_matrix[row][column - 1] + mismatchScore
                if column != 1 and row != 1:
                    if string2[row - 1] == string1[column - 1]:
                        self.edit_matrix[row][column] = self.edit_matrix[row - 1][column - 1] + matchScore
                    else:
                        temp = []
                        top = self.edit_matrix[row - 1][column] + gapScore
                        diagnol = self.edit_matrix[row - 1][column - 1] + mismatchScore
                        left = self.edit_matrix[row][column - 1] + gapScore
                        temp.append((top, left, diagnol))
                        if max(temp[0]) < 0:
                            self.edit_matrix[row][column] = 0
                        else:
                            self.edit_matrix[row][column] = max(temp[0])
                        temp = []
        self.string1 = string1
        self.string2 = string2
        self.gapScore = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore

    def getAlignment(self):
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.

        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ]

        In cases where there is a tie between optimal sub-alignments use the following rule:
        Let (i, j) be a point in the edit matrix, if there is a tie between possible sub-alignments
        (e.g. you could choose equally between different possibilities), choose the (i, j) to (i-1, j-1)
        (match) in preference, then the (i, j) to (i-1, j) (insert in string1) in preference and
        then (i, j) to (i, j-1) (insert in string2).
        """
        # Code to complete - generated by traceback through matrix to generate aligned pairs

        alignment = []
        max_index = np.unravel_index(np.argmax(self.edit_matrix, axis=None), self.edit_matrix.shape)
        row = max_index[0]
        column = max_index[1]
        while row > 0 and column > 0:
            # print(row, column)
            if self.string2[row - 1] == self.string1[column - 1]:
                alignment.append((row - 1, column - 1))
                row -= 1
                column -= 1
            # Check if it is a mismatch or gap
            else:

                # Check gap

                # Mismatch - Diagonal
                if self.edit_matrix[row - 1][column - 1] + self.mismatchScore == self.edit_matrix[row][column]:
                    row -= 1
                    column -= 1

                # Gap- Side
                if self.edit_matrix[row][column - 1] + self.gapScore == self.edit_matrix[row][column]:
                    column -= 1

                # Gap- Up
                if self.edit_matrix[row - 1][column] + self.gapScore == self.edit_matrix[row][column]:
                    row -= 1

        # Shitty reversal because I want to pass the unit test
        pre = list(reversed(alignment))
        alignment = [tuple(reversed(x)) for x in pre]
        return alignment

    def getMaxAlignmentScore(self):
        """ Returns the maximum alignment score
        """
        # Code to complete
        return np.amax(self.edit_matrix)


def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.

    Maps the string in both its forward and reverse complement orientations.

    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it,
    and note any drawbacks this might have.

    Potential ways to accelerate this code may be to
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
            queryStringStart = max(0, seedCluster.minX - config.column)  # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.column)  # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]

            targetStringStart = max(0, seedCluster.minY - config.column)  # Inclusive coordinate
            targetStringEnd = min(len(targetString),
                                  seedCluster.maxY + config.k + config.column)  # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]

            # print "target_aligning", targetStringStart, targetStringEnd, targetSubstring
            # print "query_aligning", queryStringStart, queryStringEnd, querySubstring

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
        rMap = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
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
        self.column = 100
        self.gapScore = -2
        self.matchScore = 3
        self.mismatchScore = -3
        self.logLevel = "INFO"


def main():
    # Read parameters
    config = Config()

    # Parse the inputs args/options
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
                                              " each other in both target and query. Default=%s" % config.l,
                        default=config.l)
    parser.add_argument("--column", dest="column",
                        help="Add this many bases to the prefix and suffix of a seed cluster in the"
                             " target and query sequence. Default=%s" % config.column, default=config.column)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" %
                                                            config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" %
                                                                config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" %
                                                                      config.mismatchScore,
                        default=config.mismatchScore)
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
                ((time.time() - startTime), len(minimizerIndex.minimizerMap), minmerInstances))

    # Open the query files
    alignmentScores = []  # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, xrange(sys.maxint)):
            print
            queryIndex
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.getMaxAlignmentScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" %
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore))
            # Comment this out to test on a subset
            # if queryIndex > 100:
            #    break

    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" %
                    (time.time() - startTime, float(sum(alignmentScores)) / len(alignmentScores)))


if __name__ == '__main__':
    main()
