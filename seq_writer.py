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


def cagewriter():
    ##------------------------------------------------------------------------------------------------------------
    ## settings
    DNA = input("Enter your formatted DNA sequence now or press (enter) to generate one: ")
    filename = input('File name: ')
    linker = input("Linker filename: ")
    mg = input('Add Mg atoms? (y/n): ')
    min = input('Minimize? (y/n): ')
    inert_length = 1

    if DNA == 'y':
        ## shape input
        vertices = int(input('Vertices: '))
        bp_len = int(input('Middle sequence base pairs: '))

        ##------------------------------------------------------------------------------------------------------------
        ## generate inital middle lines
        ## array to hold all sequences as seperate rows
        st = np.empty((vertices, 5), dtype='object')

        ## iterate through array to gen sequences
        regen = True
        for row in st:
            row[2] = gen_DNA(bp_len)

        ## check initial middle sequences - gen new if necessary
        regen = True
        while regen: # this is garbage - fix it
            regen = False
            for row in st:
                for i in range(vertices):
                    regen = comp_check(row[2], st[i, 2])
                    regen = base_check(row[2])
                    if regen:
                        print('regenerating middle..')
                        row[2] = gen_DNA(bp_len)

        ## generate end sequences
        for i in range(vertices):
            if i == 0:
                st[i, 0] = conv_base(st[vertices-1, 2][:(bp_len//2 - 1)])
                st[i, 4] = conv_base(st[vertices-1, 2][(bp_len//2):bp_len])
            if i == vertices - 1:
                st[i, 0] = conv_base(st[0, 2][:(bp_len//2 - 1)])
                st[i, 4] = conv_base(st[0, 2][(bp_len//2):bp_len])
            else:
                st[i, 0] = conv_base(st[i + 1, 2][:(bp_len//2 - 1)])
                st[i, 4] = conv_base(st[i + 1, 2][(bp_len//2):bp_len])

        ##------------------------------------------------------------------------------------------------------------
        ## generate inert sequences
        for row in st:
            row[1] = gen_DNA(bp_len)
            row[3] = gen_DNA(bp_len)

        ## check inert against middle sequences
        regen = True
        while regen:
            regen = False
            for row in st:
                for i in range(vertices):
                    if comp_check(row[1], st[i, 2]) or base_check(row[1]):
                        print('regenerating left inert end..')
                        row[1] = gen_DNA(bp_len)
                    regen = False
                    if comp_check(row[3], st[i, 2]) or base_check(row[3]):
                        print('regenerating right inert end..')
                        row[3] = gen_DNA(bp_len)

        ## check inert against end sequences
        while regen:
            regen = False
            for row in st:
                for i in range(vertices):
                    if comp_check(row[1], st[i, 0]) or comp_check(row[1], st[i, 4]) or base_check(row[1]):
                        print('regenerating left inert end..')
                        row[1] = gen_DNA(bp_len)
                    regen = False
                    if comp_check(row[3], st[i, 2]) or comp_check(row[3], st[i, 4]) or base_check(row[1]):
                        print('regenerating right inert end..')
                        row[3] = gen_DNA(bp_len)
    else:
        st = DNA
    #print(st)
    ##------------------------------------------------------------------------------------------------------------
    ## create file with sequence
    filename = input('File name: ')
    f = open('input/' + filename, 'w')

    ## write settings to file
    f.write('# !Name ' + filename + '\n')
    if mg == 'y':
        f.write('# !Ion Mg' + '\n')
    if min == 'y':
        f.write('# !Min True' + '\n')
    f.write('# !Linker ' + linker + '\n')

    ## write sequence to file
    f.write('>\n')
    for row in st:
        for seq in row:
            f.write(seq + ' ')
        f.write('\n')
    f.write('<')
    f.close()

    ## return filename for use in main
    return filename

## method to generate junction sequences
def junctionwriter():
    ##------------------------------------------------------------------------------------------------------------
    ## shape input
    filename = input('File name: ')
    bp_len = int(input('Middle sequence base pairs: '))
    mg = input('Add Mg atoms? (y/n): ')
    min = input('Minimize? (y/n): ')
    vertices = input('Vertices: ')

    ##------------------------------------------------------------------------------------------------------------
    ## generate outer sequences
    st = np.empty((3, 1), dtype='object')  ## will probably want to change this to accept vertices at some point

    ## iterate through array to gen sequences
    regen = True
    for row in st:
        row[0] = gen_DNA(bp_len)

    ## check initial middle sequences - gen new if necessary
    regen = True
    while regen: # this is garbage - fix it
        regen = False
        for row in st:
            for i in range(3):
                regen = comp_check(row[0][0], st[i][0])
                regen = base_check(row[0][0])
                if regen:
                    print('regenerating middle..')
                    row[0][0] = gen_DNA(bp_len)

    ##------------------------------------------------------------------------------------------------------------
    ## write to file
    f = open('input/' + filename, 'w')

    f.write('# !Name ' + filename + '\n')
    if mg == 'y':
        f.write('# !Ion Mg' + '\n')
    if min == 'y':
        f.write('# !Min True' + '\n')

    f.write('>\n')
    for row in st:
        for seq in row:
            f.write(seq + ' ')
        f.write('\n')
    f.write('<')
    f.close()
