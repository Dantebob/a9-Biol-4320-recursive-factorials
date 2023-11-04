#! /bin/python3

import sys
import re
from Bio.seq import seq

class Fasta:
    def __init__(self, sample_id, date, phenotype, nuc_seq):
        self.sample_id = sample_id
        self.date = date
        self.phenotype = phenotype
        self.nuc_seq = nuc_seq
    def __str__(self):
        return f"Sample_ID: {self.sample_id} \nDate: {self.date} \nPhenotype: {self.phenotype} \nNucleotide Sequence: {self.nuc_seq}"

def get_fasta_obj(fasta_file):
    while open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                header = re.split(r'[_\s]', line)
                sample_id.append = header[0]
                date = header[1]
            else:
                nuc_seq = seq(line)
                nuc_seq = nuc_seq.transcribe()
            aa_seq = nuc_seq.translate()
            if aa_seq[3] = "S":
                phenotype = "blue"
            elif aa_seq[3] = "R":
                phenotype = "orange"
    fasta_obj = Fasta(sample_id, date, phenotype, nuc_seq)
    return fasta_obj

def get_probability()

def print_output_message(p, n, k):
    result = 0
    output_message = 
        (f'Results\n\np (the frequency of "orange" in the population) = {p}\n'
        f'n (the number of sampled individuals) = {n}\n'
        f'k (the number of "orange" individuals in the sample set) = {k}\n\n'
        f'Probability of collecting {n} individuals with {k} being "orange" (given a population frequency of {p}) = {result}')
    return output_message
def get_factorial(n):
    if n == 1:
        return 1
    else:
        return n * get_factorial(n-1)

def __name__ == '__main__':
    fasta_file = sys.argv[0]
    trait_freq = sys.argv[1]
    results_file = sys.argv[2]


