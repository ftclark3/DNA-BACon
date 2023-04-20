from writer_dev import gen_DNA, comp_check, conv_base
import numpy as np
## program to create threeway junction sequences
## if youre reading this, ik my code is disgusting (ill fix it later)

##------------------------------------------------------------------------------------------------------------
## create array to store first pillar sequence
st1 = np.array([['                                    '], # junction backbone
                ['       '], # first u seq
                ['       '], # second u seq
                ['       '], # third u seq
                ['       '], # first addn (closest to end of pillar)
                ['       ']]) # second addn
## create array to store second pillar sequence
st2 = np.array([['                                    '], # junction backbone
                ['       '], # first u seq
                ['       '], # second u seq
                ['       '], # third u seq
                ['       '], # first addn (closest to end of pillar)
                ['       ']]) # second addn
## create array to store third pillar sequence
st3 = np.array([['                                    '], # junction backbone
                ['       '], # first u seq
                ['       '], # second u seq
                ['       '], # third u seq
                ['       '], # first addn (closest to end of pillar)
                ['       ']]) # second addn
print('arrays generated')

##------------------------------------------------------------------------------------------------------------
## generate pillar sequences
st1[0] = gen_DNA(36)
gen = True
## generate second pillar and check if if it is comp
st2[0] = gen_DNA(36)
while gen:
    gen = False
    if comp_check(st1[0], st2[0]):
        gen = True
        st2[0] = gen_DNA(36)
## generate third pillar and check if it is comp
gen = True
st3[0] = gen_DNA(36)
while gen:
    gen = False
    if comp_check(st2[0], st3[0]) or comp_check(st1[0], st3[0]):
        gen = True
        st3[0] = gen_DNA(36)
print('pillars generated')

##------------------------------------------------------------------------------------------------------------
## generate first u sequence - comp to its own pillar
st1[1] = conv_base(st1[0][0][-1:-8:-1])
st2[1] = conv_base(st2[0][0][-1:-8:-1])
st3[1] = conv_base(st3[0][0][-1:-8:-1])
print('first u generated')

## generate third u sequence - comp to the next pillar
st1[3] = conv_base(st2[0][0][1:8])
st2[3] = conv_base(st3[0][0][1:8])
st3[3] = conv_base(st1[0][0][1:8])
print('second u generated')

## generate second/middle u sequence - comp to nothing
## for st1
gen = True
st1[2] = gen_DNA(7)
while gen:
    gen = False
    for row in st1:
        if comp_check(row, st1[2][0]) and row != st1[2][0]:
            gen = True
    for row in st2:
        if comp_check(row, st1[2][0]):
            gen = True
    for row in st3:
        if comp_check(row, st1[2][0]):
            gen = True
    if gen:
        st1[2] = gen_DNA(7)
## gen and check for st2
gen = True
st2[2] = gen_DNA(7)
while gen:
    gen = False
    for row in st1:
        if comp_check(row, st2[2][0]):
            gen = True
    for row in st2:
        if comp_check(row, st2[2][0]) and row != st2[2][0]:
            gen = True
    for row in st3:
        if comp_check(row, st2[2][0]):
            gen = True
    if gen:
        st2[2] = gen_DNA(7)
## gen and check for st3
gen = True
st3[2] = gen_DNA(7)
while gen:
    gen = False
    for row in st1:
        if comp_check(row, st3[2][0]):
            gen = True
    for row in st2:
        if comp_check(row, st3[2][0]):
            gen = True
    for row in st3:
        if comp_check(row, st3[2][0]) and row != st2[2][0]:
            gen = True
    if gen:
        st3[2] = gen_DNA(7)
print('u-sequence generated')

##------------------------------------------------------------------------------------------------------------
## gen two additional comp seq with jb comp
st1[4] = conv_base(st1[0][0][19:26])
st1[5] = conv_base(st1[0][0][12:19])

st2[4] = conv_base(st2[0][0][19:26])
st2[5] = conv_base(st2[0][0][12:19])

st3[4] = conv_base(st3[0][0][19:26])
st3[5] = conv_base(st3[0][0][12:19])
print('additional sequences generated')

##------------------------------------------------------------------------------------------------------------
## write to file
filename = input('Filename: ')
f = open(filename, 'w')
f.write('# !Name ' + filename + '\n')

## write sequence to file
f.write('>\n')
for row in st1:
    f.write(row[0] + ' ')
f.write('\n')
for row in st2:
    f.write(row[0] + ' ')
f.write('\n')
for row in st3:
    f.write(row[0] + ' ')
f.write('\n')
f.write('<')
f.close()
