#! /bin/python3

import sys
import re
from Bio.Seq import Seq

class Fasta:
    def __init__(self, sample_id, date, phenotype, nuc_seq, n, k):
        self.sample_id = sample_id
        self.date = date
        self.phenotype = phenotype
        self.nuc_seq = nuc_seq
        # k = the number of occurences of a trait, n = the sample size,
        self.n = n
        self.k = k
    def __str__(self):
        return f"Sample_ID: {self.sample_id} \nDate: {self.date} \nPhenotype: {self.phenotype} \nNucleotide Sequence: {self.nuc_seq}"

    def get_factorial(self, n):
        if n == 1:
            return 1
        else:
            return n * self.get_factorial(n-1)

    def get_probability(self, p):  # p = the frequency of a given trait in a population
        p = float(p)
        n = self.n
        k = self.k
        n_fac = self.get_factorial(n)
        k_fac = self.get_factorial(k)
        nmk_fac = self.get_factorial(n - k)
        q = 1.0 - p
        probability = (n_fac / (nmk_fac * k_fac))*(p**k)*(q**(n-k))
        return probability

    def write_output_message(self, p, probability, output_file):
        n = self.n
        k = self.k
        output_message = (f'Results\n\np (the frequency of "orange" in the population) = {p}\n'
            f'n (the number of sampled individuals) = {n}\n'
            f'k (the number of "orange" individuals in the sample set) = {k}\n\n'
            f'Probability of collecting {n} individuals with {k} being "orange" (given a population frequency of {p}) = {probability}')
        with open(output_file, "w") as out_file:
            out_file.write(output_message)


def get_fasta_obj(fasta_file):
    sample_id = []
    date = []
    phenotype = []
    with open(fasta_file, "r") as fasta:
        n = 0
        k = 0
        for line in fasta:
            if line.startswith(">"):
                n += 1
                header = re.split(r'[_\s]', line)
                sample_id.append(header[0])
                date.append(header[1])
            else:
                nuc_seq = Seq(line)
                nuc_seq_t = nuc_seq.transcribe()
                aa_seq = nuc_seq_t.translate()
                if aa_seq[3] == "S":
                    phenotype.append("blue")
                elif aa_seq[3] == "R":
                    phenotype.append("orange")
                    k += 1
    fasta_obj = Fasta(sample_id, date, phenotype, nuc_seq, n, k)
    return fasta_obj


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    trait_freq = sys.argv[2]
    results_file = sys.argv[3]
    fasta_obj = get_fasta_obj(fasta_file)
    probability = fasta_obj.get_probability(trait_freq)
    fasta_obj.write_output_message(trait_freq, probability, results_file)
    print(f"results have been written to {results_file}")
