## program to check sequence compatability
import matplotlib.pyplot as plt
import numpy as np

## import sequence
sides = [['TCGCTGAGTA', 'TCAACTGCTCTCAACTGCTC', 'GCAAGTGTGGGCACGCACAC', 'CACAAATCTG'], ['TCGCTGAGTA', 'TCAACTGCTCTCAACTGCTC', 'GCAAGTGTGGGCACGCACAC', 'CACAAATCTG']]

## create figure
bbox = {'fc': '1.5', 'pad': 0} # the text bounding box
fig, axs = plt.subplots(1, 1)
props = {'ha': 'center', 'va': 'center', 'bbox': bbox}

## iterate through sides and place text
addendum = 1/40*len(sides[0][0])
for i in range(len(sides)):
    axs.text(addendum*2, 0.6, sides[i][2], props, rotation = 90)
    axs.text(addendum, 0.38, sides[i][1], props)
    axs.text(addendum, 0.82, sides[i][1], props)
    addendum += 1/40*len(sides[i][0])




#addtext(axs, {'ha': 'center', 'va': 'center', 'bbox': bbox})
plt.show()
