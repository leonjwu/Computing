# [Project 2: Searching Algorithms and Graph Path-Finding](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/)

## Overview
- Developed efficient DNA matching algorithms using rolling hash functions (Rabin-Karp) and tables including runtime analysis
- Designed route planning algorithms for special directed and weighted graphs using BFS, DFS, Djikstras.




A brief outline of the algorithm is as follows, with L the list of M N -length
sequences, and pairs the list of J K-length pairs:

- Convert all sequences and pairs into base 4 by mapping [A,C,G,T] onto
[0,1,2,3]. Running time: O(M N + JK) since the dictionary lookup
is of constant time to convert a letter into a digit.

- Create a dictionary, where the keys are the hashes of the first of each pair,
and the values are dictionaries consisting of {second of the pairs : list of
pair indices}. In each case, the hashes are checked in the dictionary and
are created/added to. This hash function converts base 4 integers into
base 10 integers. There is no need to mod the hashes since there in no
performance penalty associated with large integers, and also this means
there are no hash collisions possible since each K-length string uniquely
maps to its hashed representation. This part is O(JK), since again
the checking of the hashed pair and adding to the dictionary is
done in constant time, and since the hashing function runs in
O(K) and is repeated J times.

- The hashes for first K-mers in all sequences are calculated. O(M K).
Then, Iterate down each sequence, calculating the rolling hash as in the
Rabin-Karp method. For each sequence, check whether the hash exists
in the pairs dictionary. If so, check all corresponding pairs in the next
sequence, and if there is a match, append to the output. There is no need
to check for hash collisions since K is constant, as discussed earlier. The
code is arranged so that the K-mers in the sequences are only hashed once
to avoid any duplicate calculations. Roughly, this is O((N − K)M J) since
the rolling hash is calculated in constant time, the matches are checked in
a dictionary so is also constant time. The dependence on J here is actually
eliminated in the case where there are no exact repeated pairs, such that
index list is of length 1, which will be assumed.

- Overall, we have a worst-case run time of ((N −K)M +M N +JK +M K).
Asymptotically as N and M get large, and since N >> K and M >> K,
we have an asymptotic worst-case run time of O(N M + J). A naive
approach would loop through each sequence, and
check one character at a time for matches with all the pairs. This would
be O((N M ∗ K ∗ J)), as discussed in lectures for the worst case with many
near misses, but in this case M sequences must be checked for J patterns.
Asymptotically as N and M get large, and since N >> K and M >> K,
we have an asymptotic worst-case run time of O(N M J).

- In summary, naive search runs in O(M N J) whereas my algorithm runs
in O(M N + J), with the smaller terms omitted. Clearly the latter is
preferred since the J term is additive rather than multiplicative. Also, my
algorithm has no N K term whereas the naive search does, which could
add significant time to the naive search algorithm compared to mine.
