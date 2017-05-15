from string import maketrans

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

    reverse = maketrans('[]', #modified to allow reversing regexp 
                            '][') #modified to allow reversing regexp 
    s.translate(reverse)

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
    known_bases = set('ATUCGNatucgnKRSBDMYWVHkm[]')  #modified to allow complementing regexp 
    unknown = set(s) - known_bases
    if len(unknown) > 0:
        raise ValueError('Unknown bases passed: %s' % ', '.join(unknown))

    complement = maketrans('ATUCGNatucgnKRSBDMYWVHkm',  
                               'TAAGCNtaagcnMYSVHKRWBDmk') 
    return s.translate(complement)

def degenerate(seq):
    """ Takes a sequence and returns the a regular expression """
    

    reconversiondict = {'N': '[ATGC]',
                        'S': '[CG]',
                      'W': '[AT]',
                      'R': '[AG]',
                      'Y': '[CT]',
                      'M': '[AC]',
                      'K': '[GT]',
                      'B': '[CGT]',
                      'D': '[AGT]',
                      'H': '[ACT]',
                      'V': '[ACG]',
                      'A': 'A',
                      'C': 'C',
                      'G': 'G',
                      'T': 'T'}


    import re
    #convert degeneracies
    for r1, r2 in reconversiondict.items():
        seq = re.sub(r1,r2,seq)
    return seq

import unittest

class TestDegenerate(unittest.TestCase):
    def test_degenerate(self):
        self.assertEqual(degenerate("GNRA"),"G[ATGC][AG]A")
        
suite = unittest.TestLoader().loadTestsFromTestCase(TestDegenerate)
unittest.TextTestRunner(verbosity=2).run(suite)

# The function below has a docstring, and you are tasked with filling in
# the body of the function.

def Site_finder_fd(full_seq, subseq):
    """ Gets a list of information on all sub-sequences in a given DNA sequence
    This one just finds the subseq in the forward strand.
    Also, it outputs a dict rather than a list.
    
    Below is the requested Site_finder

    Parameters
    ----------
    full_seq : str
        Input DNA sequence
    subseq : str
        target sub-sequence (i.e. NRG)

    Returns
    -------
    list of tuple of objects
        This is a list of dicts with the following format:
            [target1, target2, target3, ...., targetN]

        where targetX has the format:
        {"start":start, "seq":exact sub_sequence, "forward": 
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
    
    subseq = degenerate(subseq);
    
    import re
    #print (full_seq)
    fmatches = [{"start":a.start(), "seq":full_seq[a.start() : a.end()], "forward":True} for a in list(re.finditer(subseq, full_seq))]
    #print("F=",fmatches)
    return fmatches
    
    #matches_list = list(re.finditer(subseq, full_seq))
    #rmatches = [(len(full_seq)-a.start(), full_seq[a.start() : a.end()], True) for a in reversed(matches_list)]
    #print("R=",rmatches)
    
import unittest
    
class TestSite_finder_df(unittest.TestCase):
    def test_find_a_simple_substring(self):
        m = Site_finder_fd("ATTCGGTAACGCAT","GTAA")
        self.assertEqual(m,[{'start': 5, 'seq': 'GTAA', 'forward': True}])
    def test_find_a_degenerate_substring(self):
        m = Site_finder_fd("ATTCGGTAACGCAT","GNRA")
        self.assertEqual(m,[{'start': 5, 'seq': 'GTAA', 'forward': True}])
            
suite = unittest.TestLoader().loadTestsFromTestCase(TestSite_finder_df)
unittest.TextTestRunner(verbosity=2).run(suite)    




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
    sites_forward = Site_finder_fd(full_seq,subseq)
    sites_reverse = Site_finder_fd(complement(full_seq), (reverse(subseq)))
    sites = [(a['start'],a['seq'],True) for a in sites_forward] + [(a['start'],a['seq'],False) for a in sites_reverse] 
    return sorted(sites)
    

class TestSite_finder(unittest.TestCase):
    def test_find_subseq(self):
        m = Site_finder("ATTCGGTAACGCATTTTTATGCGTTACCGAATATTCGGTAACGCATTTTTATGCGTTACCGAAT","GNRA")
        self.assertEqual(m,[(5, 'GTAA', True), (23, 'AATG', False), (37, 'GTAA', True), (55, 'AATG', False)]) #sorted and a list
            
suite = unittest.TestLoader().loadTestsFromTestCase(TestSite_finder)
unittest.TextTestRunner(verbosity=2).run(suite)    


def generate_matchstring(m,l):
    ''' generate a matchstring from a list of matches and a provided string length '''
    s = list(" "*l)
    for a in m:
        #print(a[0])
        #print(a[1])
        s[a['start']]=a['seq'][0]
        s[a['start']+1]=a['seq'][1]
        s[a['start']+2]=a['seq'][2]
        s[a['start']+3]=a['seq'][3]
    return "".join(s)

def FindHairpins (sequence, forward):
    ''' Gets a list of hairpins in a given DNA sequence

    Parameters
    ----------
    full_seq : str
        Input DNA sequence
    forward : bool
        Forward/Reverse

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
    Rule 1: substring length >= 12 and <=100
    Rule 2: A ‘GNRA’ motif is at the center of the subsequence 
    Rule 3: The subsequence contains at least 4 pairs of complementary nucleotides.
    Rule 4: Each half of each complementary pair must be located equidistant from the center of the ‘GNRA’ motif as shown in the diagram of the hairpin structure above.
    Rule 5: Both nucleotides of one of the complementary pairs must be a distance of <= 2 from either end of the ‘GNRA’ motif
    Rule 6: At least 70% of nucleotides in a complementarity region must be participating in a pairing interaction.

    '''
    
    if forward:
        #print("Forward")
        subseq = "GNRA"
    else:
        #print("Reverse Complement")
        subseq = reverse("GNRA")
        sequence = complement(sequence)
    
   
    #print(subseq)
    
    # Start with Rule 2: find all the GNRA motifs
    matches = Site_finder_fd (sequence,subseq)
    
    # Convert to lists so we can index easily
    l = list(sequence)
    
    hairpins = []
    
    # For each match in the string
    for m in matches:
        
        #reload matchstring each time
        matchstring = generate_matchstring(matches,len(sequence))
        lm = list(matchstring)
    
        #print(m['start'])
       
        # Rule 1 suggests that we have max string length of 100, so it's 50 comparisons
        r1 = -1
        r2 = -1
        
        found_i = 0
        found_complements = 0
        found_complement_ratio = 0
        found_r1 = 0
        found_r2 = 0
        
        complements = 0
        complement_ratio = 0
        best_possible = 1
        
        for i in range (1,48):  #from Rule 1
                
            r1 = m['start']-i    #residue 1
            r2 = m['start']+3+i  #residue 2
            
            # GNRA motif terminates search
            #if (lm[r1] != " "): break
            #if (lm[r2] != " "): break
              
            if (l[r1][0] == complement(l[r2])):
                complements += 1
                lm[r1] = "-"
                lm[r2] = "-"
            else:
                lm[r1] = "."
                lm[r2] = "."
                
            if (i==2):
                if ((lm[m['start']-1] != "-") & (lm[m['start']-2] != "-")): break
                    
            complement_ratio = (1.0 * complements) / i 
            
            #print ("%4d %4d %1s %1s %1d %100s %2d/%2d %2.2f %2.2f %4d %4d"%(r1,r2,l[r1],l[r2],l[r1][0] == complement(l[r2]),"".join(lm[r1:r2+1]),complements,i,complement_ratio,best_possible, found_r1, found_r2))
            
            if ((i>=4) & (complement_ratio > 0.7) & (complements > 4)):
                found_r1 = r1
                found_r2 = r2
                found_i  = i
                found_complements = complements
            
            if ((i>4) & (best_possible<0.7) & (complements > 4)):
                r1 = found_r1
                r2 = found_r2
                complement_ratio = found_complements/found_i
                break
               
            best_possible = (complements + (49-i)) / 50.0
            if (best_possible < 0.70):
                if (found_i):
                    r1 = found_r1
                    r2 = found_r2
                    complement_ratio = found_complements/found_i
                break
            
            # string ends terminate search
            if (r1==0): break
            if (r2==len(lm)-1): break
            
            
                
        hstring = "".join(lm[r1:r2+1])
        
        #hairpins should always be symmetric
        #assert (len(hstring))%2 == 0
        
        #print ('%4d %4d %4d %-100s %2d %2.2f' % (m['start'],r1,r2,hstring,len(hstring),complement_ratio))
        
        if (complement_ratio >= 0.70):
           if ((lm[m['start']-1] == "-") | (lm[m['start']-2] == "-")):
                hairpins.append({"start":r1, "motif":m['start'], "hairpin":hstring})
                
                ## in-line verification of rules
                
                ##Rule 1: substring length >= 12 and <=100
                assert len(hstring) >= 12, "Rule 1: substring length >= 12 and <=100"
                assert len(hstring) < 100, "Rule 1: substring length >= 12 and <=100"
                
                ## Rule 2: A ‘GNRA’ motif is at the center of the subsequence 
                assert len(Site_finder_fd(hstring,subseq)) > 0, "Rule 2: A ‘GNRA’ motif is at the center of the subsequence "
            
                ## Rule 3: The subsequence contains at least 4 pairs of complementary nucleotides.             
                complements = hstring.count("-")
                assert complements >= 8, "Rule 3: The subsequence contains at least 4 pairs of complementary nucleotides."
                
                ## Rule 4: Each half of each complementary pair must be located equidistant from the center of the ‘GNRA’ motif as shown in the diagram of the hairpin structure above. 
                startMotif = (Site_finder_fd(hstring,subseq)[0]['start'])
                assert (startMotif == (len(hstring)-4)/2),  "Rule 4: Each half of each complementary pair must be located equidistant from the center of the ‘GNRA’ motif as shown in the diagram of the hairpin structure above."
                
                ## Rule 5: Both nucleotides of one of the complementary pairs must be a distance of <= 2 from either end of the ‘GNRA’ motif
                assert (lm[m['start']-1] == "-") | (lm[m['start']-2] == "-"), "Rule 5: Both nucleotides of one of the complementary pairs must be a distance of <= 2 from either end of the ‘GNRA’ motif"
                
                ## Rule 6: At least 70% of nucleotides in a complementarity region must be participating in a pairing interaction.  
                complement_ratio = (1.0 * complements) / (len(hstring)-4) 
                assert (complement_ratio >= 0.7), 'Rule 6: At least 70% of nucleotides in a complementarity region must be participating in a pairing interaction. %2.2f' % (complement_ratio)
        
        
    return hairpins 
    return "".join(lm)

def generate_hairpinstring (hairpins, length):
    s = list(" "*length)
    for h in hairpins:
        #print(h['start'])
        #print(h['hairpin'])
        for i in range (0, len(h['hairpin'])):
            s[h['start']+i]=(h['hairpin'][i])
    return "".join(s)
    


given_seq = "ATTCGGTAACGCAT"   
example_seq =  given_seq + "TTTT" + revcomp(given_seq) 

class TestFindHairpins(unittest.TestCase):
    
    
    def test_find_forward_hairpin(self):
        h = FindHairpins(example_seq, True)
        self.assertGreater(len(h), 0)
    def test_find_reverse_hairpin(self):
        h = FindHairpins(example_seq, False)
        self.assertGreater(len(h), 0)
    def test_eyeball_strings(self):
        print ("\n")
        hf = FindHairpins(example_seq, True)
        hfs = generate_hairpinstring (FindHairpins(example_seq, True), len(example_seq))
        
        hc = FindHairpins(example_seq, False)
        hcs = generate_hairpinstring (FindHairpins(example_seq, False), len(example_seq))
        
        print (hfs + "\n" + example_seq + "\n" + hcs)        
            
suite = unittest.TestLoader().loadTestsFromTestCase(TestFindHairpins)
unittest.TextTestRunner(verbosity=2).run(suite)    





print ("\n\nFind all the hairpins in the test sequence")


test_sequence = "TTCATGGTATACTAGTCACGGTGCGCGCATAATGAAGGACTGTTCGCGTCCGTCACGTCTTAGGCCTGGAAAGGCCTAAGACGTGACGGACGCGAACAGTCCTTCATTATGCGCGCACCGTGAAAGTACCATATGATCAGGCACGAGATGGAAATCGGGAGATTTCCATCTCGTGCCTCCGAAGTAGCCACCGGCGCGCGGCTTCATCGGTGGCCCGATTATGAGGGTTTCACACGGGAGTGTGGCTAATACTCCCAAAGCCAGTTTAGTGAAACGCATAGAACAAGGCGATGACAGTGACCCATAACATGGGAGGGTTATTTAGGGATAAATAACCCTCCCATGTTATGGGTCACTGTCATCGCCTTGTTCTATGCGTTTCACTAAACTGGGCCCTGCACTGCGCGGTTCAAGCAGTAATGCTTGAACCGCGCAGTGCAGGGCTACCTAGCGGTCTCACTATTAAACAATGGGGGACCATTGTTTAATAGTGAGACCGCTAGGTATCATAGGTAGTATCGTCCTAGAGTAATCTAGGACGATACTACCTATGATCCAGCCTTTCAGACTACTAGCTAGGTCGGAAAGTCTGATGACAAACATTTGCGTCTAAACAGGGGACTGTTTAGACGCAAATGTTTGAATAACGCTGGTCGGACGCAGTATCGGTGGAGGATATCAACCCCGTGCGTATGGTCCCCTATCGTCCGGACGATAGGGGACCATACGCACGGGGTTGATATCCTCCACCGATACTGCTTATTGCGACCAGCCTGCAACTAATCTCGCGTGTTGCTCTCGAGAGCAACACGCGAGATTAGTTGGACGTCAATACGCAGAAATCTAGATTTCTGCGTATTGACGTCATCAGTGCCTGTCCGGAACACTCGCGAGAGTTAGTCACGGACAGGCCTTGGTACTCGACAAAGAAAGTCTGTTTCCATATATGGAAACACATGAGCTGTTTCTTTCAGTTTCCAGACCAATGTAGGTAAGCGTCGGGGACGACGCTTACCTACAAAGGTCTGGTTAGAACAAGGCACCCCTTAGGATCCTAAGGGGTGCCTTGTTCATTTCCTGGAGTGCCATAAGAGCTCACCGGTGGACTTTATTGAGCCACTAAGTACATGAGTAAGCGGACTGACCAAGATTGTTAAGATGCAAGCTTGCATCTTAACAATCTTGGTCAGTCCGCTTACTCATGTACTTAGTGGCTCAATAAAGTCCACCGGTGAGCTCTAAAGGACCTCACGGTATTCGGTGAGGCGCGTGGTCATAAAGGATTCGAGAGAATCCTTTATGACCACGCGCCTCACCGTCAGTTTCTGATATTGGGCACTGACCCAGACCTACCAAAGGGCGTCGAGTAGATCGCCGCGCTCAATTGAGCACGTTCGACGTCGAACGTGCTCAATTGAGCGCGGCGATCTACTCGACGCCCTTTGGTAGGTCTGGGTCAGTGCCCAGTCAAAGACTATAAAATCCCATGCCCCGTCCTGCAGGACTTAGGGTACGGGGGCCGGGTAGTAGCCCGCGCGCGGCGAGAGCCGCCGGCCCATCATCGGGCGCGCTGGGGAAACTTGTGTCTCTAGAGACACAAGTTTCCCCAGTTGGTATCTTGGCGTGACCCTCACGGGGGCGACCCCGTGAGGGTCACGCCAAGATACCAATAACGCTGATGTTTGAGTTGTTGAGTGAAAACTCAACAACTCAAACATCAGCGTTAGCGGACACACTGTTCCGGCCGCGCCTGTGTGACAAGTTATTAGCCTATCCAGGGATGAGCCAAGAGTCCCACCCGGTGTGTAAGTTATCCTGATTTAGTGCAATGAGGGATCATTGCACTAAATCAGGATAACTTACACACCGGGTGGGACTCTTGGCTCATCAATAATCGGATAGGTCCGTAAGCCCGGACTGATCTCTGGCGGGGACGCCAGAGATCAGTCCGGGCTTACGATAGGCCTCTATATACCTCCCCGCGAGAGCGGGGCTATCCGGAGATATATGGACGCGTGTTGCCAGTACGTTGCAAAACGTACTGGCAACACGCGACCATAGGCTTCTTGCGCTGATCTGGTAACATTTAAACGTGAGACTTCTAACCCTCGTCAGAGATGACGAGGGTTAGAAGTCTCACGTTTAAATGTTACCAGATCAGCGCAAGAAGCCTATGGTCGCCTAAGACAGGTGAATCGCTCGAATGAAACAGAAATGTTTCATTCGAGCGCGGATTCTGTCCACTTAGTTTACCTCGCTTGCGACAGGGCAACCTAAATGGAGCGAACGCTGTTGCCTCGAGTGTTCATGATATCTACTTTATCCGTACGATTACCATATCAGAGCCCCGTGTCAACAAGTTTGAGATGGCTTACCGCAGCATTGGAAAATGCTGCGGTAAGCCATCTCAAACTTGTTGACACGGGGCTCTGATATGGTAATCGTACGGATAAAGTAGATATCATGAACACTCGAGGCAAGCCCATCGCCCAGAGCTACTCCAAAACAACGTCGCGGAGTAGTCCCCCAAAGCGCGCTCTCGCGCGAGAGCGCGCTTTGGGGGACTACTCCGCGACGTTGTTTTGGAGTAGCGGGTAGCGGGTCTCTATTCGCTTATTGGCGATTCGAATCGCCAATAAGCGAATAGAGGTTCTATGCCTGCATAACCTTCGGAGGAAACTCCGAAGGTTATGCAGGCATAGAACCTCGCTTTGTCCTCCAGGCCGGGTAGGGAACTACCCGGCCTGGAGGACAAAGCACGAGGTGGAGACTGCTACTAGGGATAGTAGCAGTCTCCACCTCGTCGTCTACAATCGTACTTGCACTACGAGTCATGACTATTTACGAAACGACGTGAGTCGTTTCGTAAATAGTCATGACTCGTAGTGCAAGTACGATTGTAGACGGGGCGACTCGGGTAGCTGGCTGAAGGGATTCAGCCAGCTACCCGAGTCGCCCTGGACTAATTGATATTAGCTAATATCAATTAGTCCA"

hf = FindHairpins(test_sequence, True)
hfs = generate_hairpinstring (FindHairpins(test_sequence, True), len(test_sequence))

hc = FindHairpins(test_sequence, False)
hcs = generate_hairpinstring (FindHairpins(test_sequence, False), len(test_sequence))

print (hfs + "\n" + test_sequence + "\n" + hcs)
