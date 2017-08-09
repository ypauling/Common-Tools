#!/usr/local/bin/env python

##### This is a script modified from Abby #####
##### Author: Bing Yang                   #####
##### Date: 2017-07-27                    #####

from __future__ import print_function, division
import sys
import re
import regex

###Function reverse_complement(seq):
###    Description: The function to obtain the reverse complement
###      sequence of the seq.
###    Parameters:
###      seq: a string representing the sequence
###    Return:
###      A string representing the reverse complement


def reverse_complement(seq: str) -> str:
    reversed = seq.upper()[::-1]
    base_paring = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    ret = ''.join([base_paring[i] for i in reversed])
    return ret


###Function get_targets(fasta):
###    Description: The function to get all the sites
###      that Cas9 will recognize and interrogate for base-pairing
###      with the crRNA. It searches for both sense and anti-sense
###      target site, which are defined as 21bp followed by a 3 prime
###      'GG' dinucleotide. If there are 'N's in the fasta sequence
###      (or any other degenerate nucleotides), they will not be
###      identified as target sites, since all of the testing
###      criteria used in this program use only the 4 regular bases.
###    Parameters:
###      fasta: a string representing the search sequence
###    Return:
###      A list containing all the matched target sites in the fasta
###      sequence

def get_targets(fasta: str) -> list:

    targets = regex.findall('([ACGT]{21}GG|CC[ACGT]{21})',
            fasta, overlapped=True)
    return(targets)


###Functions below:
###    Descriptions: All functions below test a 23bp CRISPR target
###      site for characteristics that are either enriched in
###      high-functioning or low-functioning crRNA sequences
###      as determined by three studies
###      (Ren et al, Cell reports 2014;
###      Doench et al, Nature biotechnology 2014;
###      Wong et al, Genome Biology 2015).
###      Each test function returns a score of either 0 or 1.
###      Score 1 is returned either when a target site is enriched
###      in high-functioning crRNA sequences or depeleted in
###      low-funcitoning crRNA sequences. Score 0 is returned if
###      vice versa.
###    Parameters:
###      target: a string representing the target site sequence.
###    Return:
###      A score either 0 or 1.

def PAM_proximal_GC(target: str) -> int:

    PAM_proximal = target[15:20]
    gc_count = sum([1 if i in set(['C','G']) else 0 for i in PAM_proximal])
    if gc_count > 3:
        return 1
    return 0

def overall_GC(target: str) -> int:

    crRNA_seq = target[0:20]
    gc_count = sum([1 if i in set(['G','C']) else 0 for i in crRNA_seq])
    gc_content = (gc_count/len(crRNA_seq))*100
    if gc_content < 80 and gc_content > 30:
        return 1
    return 0

def seed_TTT(target: str) -> int:

    seed = target[6:20]
    match = re.search('TTT', seed)
    if match == None:
        return 1
    return 0

def repeats(target: str) -> int:

    crRNA_seq = target[0:20]
    match = re.search('[A]{5}|[C]{5}|[G]{4}|[T]{4}', crRNA_seq)
    if match == None:
        return 1
    return 0

def base_19_check(target: str) -> int:

    if target[18] == 'T':
        return 0
    return 1

def undesirable_base_20(target: str) -> int:

    if target[19] in set(['C','T']):
        return 0
    return 1

def desirable_base_20(target: str) -> int:

    if target[19] == 'G':
        return 1
    return 0

def undesirable_base_16(target: str) -> int:

    if target[15] == 'G':
        return 0
    return 1

def desirable_base_16(target: str) -> int:
    if target[15] == 'C':
        return 1
    return 0
