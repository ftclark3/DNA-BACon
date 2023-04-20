from random import choice
import numpy as np
## Program to generate acceptable DNA sequences for use in readin.py

##------------------------------------------------------------------------------------------------------------
## function to generate random DNA sequence
def gen_DNA(length):
    DNA = ''
    for count in range(length):
        DNA += choice('CGTA')
    return DNA

## function to check if two strings are complementary
def comp_check(first, second):
    comp_count = 0
    for i in range(len(first)):
        if first[i] == 'A' and second[i] == 'T':
            comp_count += 1
        if first[i] == 'T' and second[i] == 'A':
            comp_count += 1
        if first[i] == 'G' and second[i] == 'C':
            comp_count += 1
        if first[i] == 'C' and second[i] == 'G':
            comp_count += 1
        else:
            comp_count += 0
    if comp_count/len(first) >= 0.5:
        return True

## function to check for triple repeats and high GC content
def base_check(sequence):
    if 'AAA' in sequence or 'GGG' in sequence or 'CCC' in sequence or 'TTT' in sequence:
        return True
    gc_count = 0
    for i in sequence:
        if i == 'G' or i == 'C':
            gc_count += 1
    if gc_count//len(sequence) >= 0.6:
        return True

## function to convert base to complementary nucleotide
def conv_base(seq):
    comp_seq = ''
    for base in seq:
        if base == 'A':
            comp_seq += 'T'
        if base == 'T':
            comp_seq += 'A'
        if base == 'G':
            comp_seq += 'C'
        if base == 'C':
            comp_seq += 'G'
    return comp_seq
