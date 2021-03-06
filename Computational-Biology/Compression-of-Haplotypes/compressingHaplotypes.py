import sys
import numpy as np
import challenge1 as cA


"""
Name: Allysia Mak
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
N. 

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort 
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the 
index in X of the ith string ordered by jth reverse prefix. To break ties 
(equal prefixes) the ordering of the strings in X is used. 

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings. 
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> 

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the 
symbol in each sequence at index j-1.  

Question 1: In terms of M and N what is the asymptotic cost of your algorithm?
    Answer: O(NM)
"""
def constructReversePrefixSortMatrix(X):
    #Creates the Mx(N+1) matrix

    k0 = [x for x in range(0,len(X))]
    A = np.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1], dtype=int)
    #Tests for empty array
    if len(X) == 0:
        A[:,:] = 0
        return A
    #Set first column to 0-len(X)
    A[:,0] = k0
    #Tests for empty strings 
    if len(X[0]) == 0: return A    
    A = np.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1], dtype=int) 
    A[:,0] = k0

    for k in range(len(X[0])):
        zeros = []
        ones = []
        if k == 0:
            ak = range(len(X))
        else:
            ak = column
        for i in ak:
            i = int(i)
            #if value is 0, append index
            if int(X[i][k]) == 0:
                zeros.append(i)
            #if value is 1, append index
            if int(X[i][k]) == 1:
                ones.append(i)
        column = zeros + ones
        A[:,k+1] = column    
    return A

"""Problem 2: 

Following on from the previous problem, let Y be the MxN matrix such that for 
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X. 

Hint: You can either use your solution to constructReversePrefixSortMatrix() 
or adapt the code from that algorithm to create Y without using 
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm?
    Answer: O(NM)
"""
def constructYFromX(X):
    #Creates the MxN matrix
    #Test for empty array
    if len(X) == 0: return np.zeros(shape=[len(X)])
    #Test for empty strings
    if len(X[0]) == 0: return ['' for i in range(len(X))]
    Y = np.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0]) ], dtype=int)
    #Gets array A (reverse prefix matrix)
    A = constructReversePrefixSortMatrix(X)

    for j in range(len(X[0])):
        for i in range(len(X)):
            Y[i,j] = X[A[i,j]][j]

    return Y

"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y, 
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?
     Answer: O(NM)
Question 3b: What could you use the transformation of Y for? 
Hint: consider the BWT.
     Answer: Y is a more efficiently compressable than X and can be transformed back into the original data (X) 
Question 3c: Can you come up with a more efficient data structure for storing Y?
     Answer: Perhaps an adaptation of a compressed suffix array for the sorted prefix array Y, using its neighbors 
"""
def constructXFromY(Y):
    #Creates the MxN matrix
    #Test for empty array
    if len(Y) == 0: return []
    #Test for empty strings 
    if len(Y[0]) == 0: return ['' for i in range(len(Y))]
    

    X = np.empty(shape=[len(Y), 0 if len(Y) == 0 else len(Y[0]) ], dtype=int)
    A = np.empty(shape=[len(Y), 1 if len(Y) == 0 else len(Y[0]) ], dtype=int)
    #set first column of X to first column of Y
    X[:,0] = Y[:,0]
    #starting after first column of X since input straight from Y
    for k in range(len(X[0])-1):
        zeros = []
        ones = []
        if k == 0:
            Aatk = range(len(Y))
        #rebuilds A from X
        for i in Aatk:
            i = int(i)
            if int(X[i][k]) == 0:
                zeros.append(i)
            if int(X[i][k]) == 1:
                ones.append(i)
        Aatk = zeros + ones
        A[:,k] = Aatk
        #rebuilds X from A and Y
        for y in range(len(Y[:,k])):
            xi = Aatk[y]
            X[xi,k+1] = Y[y,k+1]
        
    return map(lambda i : "".join(map(str, i)), X) #Convert back to a list of strings

"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared 
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because 
both end with "10" but not both "110" or both "010". 

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N, 
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j] 
and X[A[i-1,j]][:j].  

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints: 

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1 
and thesymbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix 
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is 
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm? 
    Answer: O(NM)

"""
def constructCommonSuffixMatrix(A, X):


    D = np.zeros(shape=A.shape, dtype=int) #Creates the Mx(N+1) D matrix
    #Test for empty array
    if len(X) == 0:
        D[:,:] = 0
        return D
    D[:,0] = [0 for x in range(0,len(X))]
    #Test for empty strings, return column of zeros
    if len(X[0]) == 0: return D
                                    

    for k in range(len(X[0])):
        p = k + 1
        q = k + 1
        j = k + 1
        a = []
        b = []
        d = []
        e = []
        #to store j-p, j-q values for D
        dsub = []
        esub = []
        if k == 0:
            ak = range(len(X))
            dk = [0 for i in range(len(D))]
        for i in range(len(X)):
            ak_i = ak[i]
            dk_i = dk[i]            
            if dk_i > p:        
                p = dk_i                

            if dk_i > q:
                q = dk_i
            # if X at index given by list of A indicies ak_i is 0
            if int(X[ak_i][k]) == 0:
                a.append(ak_i)
                d.append(p)
                dsub.append(j-p)                
                p = 0
            else:
                b.append(ak_i)
                e.append(q)
                esub.append(j-q)
                q = 0
        ak = a + b
        dk = d + e
        #concat subtracted values to put into divergence array
        dksub = dsub + esub
        D[:,k+1] = dksub
    
    return D

    
"""Problem 5.
    
For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.
    
The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.
    
Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?
     Answer: O(MN) if # long matches < MxN, else O(# long matches)
 
Question 5b: Can you see any major time efficiencies that could be gained by
refactoring?
     Answer: Arrays A and D could be built simultaneously
    
Question 5c: Can you see any major space efficiencies that could be gained by
refactoring?
    Answer: You may be able to take in compressed Y as an argument instead of X and uncompress it as needed
            Another thing you could do, since A and D can be constructed simultaneously, is compute a column of A and D at a time so that you don't need to store the entire arrays 

Question 5d: Can you imagine alternative algorithms to compute such matches?,
if so, what would be the asymptotic cost and space usage?
     Answer: You could use a seed algorithm for indexing sequences, not sure how fast it would be but it would require more space

"""
def getLongMatches(X, minLength):
    assert minLength > 0
    
    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)
    
    #For each column, in ascending order of column index
    for j in xrange(1, 0 if len(X) == 0 else len(X[0])):
        #Working arrays used to store indices of strings containing long matches
        #b is an array of strings that have a '0' at position j
        #c is an array of strings that have a '1' at position j
        #When reporting long matches we'll report all pairs of indices in b X c,
        #as these are the long matches that end at j.
        b, c = [], []
        
        #Iterate over the aligned symbols in column j in reverse prefix order
        for i in xrange(len(X)):
            #For each string in the order check if there is a long match between
            #it and the previous string.
            #If there isn't a long match then this implies that there can
            #be no long matches ending at j between sequences indices in A[:i,j]
            #and sequence indices in A[i:,j], thus we report all long matches
            #found so far and empty the arrays storing long matches.
            if D[i,j] < minLength:
                for x in b:
                    for y in c:
                        #The yield keyword converts the function into a
                        #generator - alternatively we could just to append to
                        #a list and return the list
                        
                        #We return the match as tuple of two sequence
                        #indices (ordered by order in X) and coordinate at which
                        #the match ends
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []
            
            #Partition the sequences by if they have '0' or '1' at position j.
            if X[A[i,j]][j] == '0':
                b.append(A[i,j])
            else:
                c.append(A[i,j])
        
        #Report any leftover long matches for the column
        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)
