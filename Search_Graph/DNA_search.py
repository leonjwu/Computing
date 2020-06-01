# DNA Matching Algorithm

def char2base4(S):
    """ Convert gene sequence to base 4 string
    """
    c2b = {}
    c2b['A'] = '0'
    c2b['C'] = '1'
    c2b['G'] = '2'
    c2b['T'] = '3'

    L = ''
    for s in S:
        L += c2b[s]
    return L


def hash10(S, base):
    """Convert list S to base-10 number where
    base specifies the base of S
    """
    f = 0
    for s in S[:-1]:
        f = base*(int(s)+f)
    f += int(S[-1])
    return f


def pairSearch(L,pairs):
    """Find locations within adjacent strings (contained in input list,L)
    that match k-mer pairs found in input list pairs. Each element of pairs
    is a 2-element tuple containing k-mer strings
    """

    # Convert L to base 4 (Note that the input is being modified)
    n = len(L)
    for i in range(n):
        L[i] = char2base4(L[i])

    # Convert each pairs to base 4
    pairs4 = []
    n = len(pairs)
    for i in range(n):
        pairs4.append((char2base4(pairs[i][0]), char2base4(pairs[i][1])))

    # Create a dictionary, with keys as hashes of the first of the pairs,
    # and values as lists of tuples, with the second of the pairs with its index in pairs
    pairs_hash_dict = {}
    for index, pair in enumerate(pairs4):
        # Calculate hash of the base 4 k-mers in pairs
        hashed0 = hash10(pair[0], 4)
        hashed1 = hash10(pair[1], 4)
        if hashed0 in pairs_hash_dict:
            if hashed1 in pairs_hash_dict[hashed0]:
                pairs_hash_dict[hashed0][hashed1].append(index)
            else:
                pairs_hash_dict[hashed0].update({hashed1: [index]})
        else:
            pairs_hash_dict[hashed0] = {hashed1: [index]}


    # Initialize
    k = len(pairs[0][0])
    N = len(L[0])  # Assumes all l in L are of length N
    locations = []
    first_digit = [int(l[0]) for l in L]

    # Compute first k-mer hash for all sequences
    hashes = [hash10(l[:k], 4) for l in L]

    # Main loop across base 4 strings in l for loop below
    for j in range(k-1, N):

        # Loop across l in L
        match = False
        for i, l in enumerate(L):

            # Separate case for first step
            if j != k-1:
                # Complete the rolling hash calculation
                hashes[i] = 4 * (hashes[i] - first_digit[i]*4**(k-1)) + int(l[j])

                # Store first base 4 digit to subtract at next iteration
                first_digit[i] = int(l[j-k+1])

            # If previous sequence matched with first of a pair
            if match:
                # Check for a match in this sequence
                if hashes[i] in pairs_hash_dict[match]:
                    index_list = pairs_hash_dict[match][hashes[i]]
                    # Below loop accounts for EXACT duplicates in pairs
                    for index in index_list:
                        locations.append([j-k+1, i-1, index])

            # Check if k-mer exists in first of pairs in pairs
            if hashes[i] in pairs_hash_dict:
                # Check pair2 in next sequence
                match = hashes[i]
            else:
                match = False

    return locations


if __name__ == "__main__":
    pass
