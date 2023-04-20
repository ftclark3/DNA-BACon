from writer_dev import gen_DNA, comp_check, conv_base, base_check
import numpy as np

## program to create threeway junction sequences
##------------------------------------------------------------------------------------------------------------
## shape input
bp_len = int(input('Number of bp for middle sequence: '))

##------------------------------------------------------------------------------------------------------------
## generate outer sequences
st = np.empty((3, 1), dtype='object')

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
filename = input('File name: ')
f = open(filename, 'w')
mg = input('Add Mg atoms? (y/n): ')
min = input('Minimize? (y/n): ')

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
