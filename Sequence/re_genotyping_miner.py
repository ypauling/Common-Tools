#! /usr/local/bin/env python

### A script to search for restriction enzyme cut sites and
### compare the cutting sites among several sequences

### Author: BING YANG
### Date: 2017-08-09
### Usage:
###     python3 re_genotyping_miner.py <INPUT.fasta>

from __future__ import print_function,division
import re
import sys
from Bio import SeqIO
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

def main(argv: list) -> None:
    if len(argv) < 2:
        print('Please specify the input fasta file.')
        exit(1)

    ### Prepare for the analysis
    amb = IUPACAmbiguousDNA()
    seq_ids = []
    sites_results = []
    output = dict()

    ### Read the fasta file
    ### Each fasta sequence should have a fasta description
    for seq_item in SeqIO.parse(argv[1], 'fasta', alphabet=amb):
        seq_ids.append(seq_item.id)
        ana = Restriction.Analysis(
        Restriction.AllEnzymes, seq_item.seq, linear=True
        )
        sites_results.append(ana.full())

    for enzyme in Restriction.AllEnzymes:
        sites = [r[enzyme] for r in sites_results]
        nub_sites = [len(s) for s in sites]
        ### Check if the number of sites are the same
        if nub_sites.count(nub_sites[0]) != len(nub_sites):
            output[str(enzyme)] = nub_sites

    print_fmt = '{:>20}'*(len(seq_ids)+1)
    print(print_fmt.format('Enzyme Name', *seq_ids))
    for k, v in sorted(output.items()):
        print(print_fmt.format(k, *v))

if __name__ == '__main__':
    main(sys.argv)
