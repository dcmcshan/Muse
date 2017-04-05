# DNA is comprised of A, T, G, and C letters. Sometimes it is useful to
# compress DNA sequence information for simplicity. The conversion dictionary
# provided here shows how DNA letters can be compressed. For example, the
# sequences AAA, TAA, GAA, and CAA can all be captured using the expression
# NAA.

conversiondict = {'N': '(A|T|G|C)',
                  'S': '(C|G)',
                  'W': '(A|T)',
                  'R': '(A|G)',
                  'Y': '(C|T)',
                  'M': '(A|C)',
                  'K': '(G|T)',
                  'B': '(C|G|T)',
                  'D': '(A|G|T)',
                  'H': '(A|C|T)',
                  'V': '(A|C|G)',
                  'A': 'A',
                  'C': 'C',
                  'G': 'G',
                  'T': 'T'}

# In nature, DNA exists in a 'double stranded' form, such that two sequences
# are paired using the complementarity rules (A:T and G:C)
# Below is a representation of a double stranded DNA sequence:
#
#           5' - ATTGCGCAATAGCCGAATGCAGCA - 3'    <--- "Forward Strand"
#                ||||||||||||||||||||||||
#           3' - TAACGCGTTATCGGCTTACGTCGT - 5'    <--- "Reverse Strand"
#
# As indicated in the diagram above, each of the two strands is labeled as
# "forward" or "reverse" depending on the orientation of it's "5' end". For
# this exercise, the meaning of "5' end" and "3' end" is not important. What
# matters is that you can use the fuctions below to convert an input sequence
# into it's "reverse complement" if you find this functionality necessary
# to solve the below problem.

def revcomp(s):
    """Return the reverse complement of the given DNA or RNA string

    Parameters
    ----------
    s : str
        DNA or RNA string

    Returns
    -------
    str
        Rever complemented DNA or RNA string

    Notes
    -----
    Degenerate bases (W, Y, N, etc) are allowed

    Raises
    ------
    ValueError
        Unknown character(s) passed
    """
    return reverse(complement(s))


def reverse(s):
    """Reverses a string

    Parameters
    ----------
    s : str
        String to reverse

    Returns
    -------
    str
        Reversed string
    """
    return s[::-1]

def complement(s):
    """Return the complement of the given DNA or RNA string

    Parameters
    ----------
    s : str
        DNA or RNA string

    Returns
    -------
    str
        Complemented DNA or RNA string

    Notes
    -----
    Degenerate bases (W, Y, N, etc) are allowed

    Raises
    ------
    ValueError
        Unknown character(s) passed
    """
    known_bases = set('ATUCGNatucgnKRSBDMYWVHkm')
    unknown = set(s) - known_bases
    if len(unknown) > 0:
        raise ValueError('Unknown bases passed: %s' % ', '.join(unknown))

    complement = str.maketrans('ATUCGNatucgnKRSBDMYWVHkm',
                               'TAAGCNtaagcnMYSVHKRWBDmk')
    return s.translate(complement)


# The function below has a docstring, and you are tasked with filling in
# the body of the function.

def Site_finder(full_seq, subseq):
    """ Gets a list of information on all sub-sequences in a given DNA sequence

    Parameters
    ----------
    full_seq : str
        Input DNA sequence
    subseq : str
        target sub-sequence (i.e. NRG)

    Returns
    -------
    list of tuple of objects
        This is a list of tuples with the following format:
            [target1, target2, target3, ...., targetN]

        where targetX has the format:
        (positon from start, exact sub-sequence,
         strand bool(True=Fwd, False=Rev))

        List is sorted by the position from start ON FWD STRAND,
        least to greatest.

    Raises
    ------
    ValueError
        input sequence must contain only A, T, G, C letters
        subseq is greater than 6 bases
        subseq contains non-nucleotide characters

    Notes
    -----
    The input sequence is the Fwd strand. It's reverse complement is the Rev
    strand.
    The exact sub-sequence returned can be comprised of only ATGC letters.
    """
    pass
