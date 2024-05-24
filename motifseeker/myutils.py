"""
Utilities for mypileup
"""
import sys

def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write("{msg}\n".format(msg=msg) )
	sys.exit(1)

def GetPFM(sequences):
    """ Compute the PFM for a set of sequences
    
    Parameters
    ----------
    sequences : list of str
        List of sequences (e.g. binding_sites)
    
    Returns
    -------
        pfm : 2d np.array
        
    Assumes all sequences have the same length
    """
    pfm = np.zeros((4, len(sequences[0])))
    for seq in sequences:
        for j, char in enumerate(seq):
            pfm[nucs[char], j] += 1

    return pfm

def GetPWM(binding_sites, background_freqs=[0.25, 0.25, 0.25, 0.25]):
    """ Compute the PWM for a set of binding sites
    
    Parameters
    ----------
    binding_sites : list of str
        List of sequences 
    background_freqs: list of float
        Background frequency of A, C, G, T
    
    Returns
    -------
        pwm : 2d np.array
        
    Assumes all sequences have the same length
    """
    pwm = np.zeros((4, len(binding_sites[0])))
    pfm = GetPFM(binding_sites)
    pfm = pfm + 0.01 

    for j in range(pfm.shape[1]):
        col_sum = np.sum(pfm[:, j])
        for i in range(4):
            p_ij = pfm[i, j] / col_sum
            pwm[i, j] = math.log2(p_ij / background_freqs[i])
            
    return pwm

def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    for i, nucleotide in enumerate(sequence):
        row_index = nucs[nucleotide]
        score += pwm[row_index, i]

    return score

def RandomSequence(n):
    """ Generate a random string of nucleotides of length n
    
    Parameters
    ----------
    n : int
       Length of random string to generate
       
    Returns
    -------
    seq : str
       Random nucleotide string
    """
    seq = ""
    for i in range(n):
        seq += ["A","C","G","T"][random.randint(4)]
    return seq

def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 
    thresh = np.percentile(null_dist, 100 * (1 - pval))
    return thresh

def ReverseComplement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    revcomp = ''.join(complement[nt] for nt in reversed(sequence))
    return revcomp

def FindMaxScore(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
    max_score = -1*np.inf
    n = pwm.shape[1]
    max_score = -np.inf

    def scan_seq(seq):
        max_seq_score = -np.inf
        for i in range(len(seq) - n + 1):
            subseq = seq[i:i+n]
            score = ScoreSeq(pwm, subseq)
            max_seq_score = max(max_seq_score, score)
        return max_seq_score

    max_score = max(scan_seq(sequence), scan_seq(ReverseComplement(sequence)))

    return max_score

def ComputeNucFreqs(sequences):
    """ Compute nucleotide frequencies of a list of sequences
    
    Parameters
    ----------
    sequences : list of str
       List of sequences
       
    Returns
    -------
    freqs : list of float
       Frequencies of A, C, G, T in the sequences
    """
    freqs = [0.25, 0.25, 0.25, 0.25]
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    # Count each nucleotide
    for seq in sequences:
        for nucleotide in seq:
            if nucleotide in count:
                count[nucleotide] += 1
                total += 1

    # Calculate frequencies
    if total > 0:
        freqs = [count[n] / total for n in ['A', 'C', 'G', 'T']]
    else:
        freqs = [0.25, 0.25, 0.25, 0.25]  # Default to equal distribution if no valid data

    return freqs

def RandomSequence(n, freqs):
    """ Generate a random string of nucleotides of length n
    
    Use the given nucleotide frequences
    
    Parameters
    ----------
    n : int
       Length of random string to generate
    freqs : list of float
       List of frequencies of A, C, G, T
       
    Returns
    -------
    seq : str
       random sequence of length n with the specified allele frequencies
    """
    seq = "A"*n
    nucs = ['A', 'C', 'G', 'T']
    seq = ''.join(random.choices(nucs, weights=freqs, k=n))
    return seq

def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       pval% of null_dist should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0
    null_dist.sort()
    index = int((1 - pval) * len(null_dist))
    thresh = null_dist[index - 1] if index > 0 else null_dist[0]
    return thresh

def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    table = [[peak_motif, peak_total - peak_motif],
             [bg_motif, bg_total - bg_motif]]

    _, pval = scipy.stats.fisher_exact(table)
    return pval

# for i in range(len(PWMList)):
    # pwm = PWMList[i]
    # thresh = pwm_thresholds[i]
    # num_peak_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
    # num_bg_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
    # pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
