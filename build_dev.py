from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from math import tan
import build_dna2b
import sys, numpy
import random
import copy

''' Functions to build a DNA cage from DNA walls '''

##------------------------------------------------------------------------------------------------------------
## function to bond walls of cage
def cage_builder(wall_list, edge, chain_length, linker_start, linker_end, resnum_start, resnum_end):

    # calculate geometric measurements
    iterator = len(wall_list)
    angle = 2*numpy.pi/iterator
    radius = edge/(2*numpy.tan(angle/2))

    # iterate through wall_list and place each wall
    # also mark the DNA atoms needed to connect the walls
    connectme = []
    for i in range(iterator):

        it = wall_list[i]

        # get new posn
        newpos = (radius, 0.0, 0.0, 0.0)
        (x, y, z, q) = newpos - transform.get_centroid(it)

        # move wall into place
        transform.translate_structure(it, 0, 0, radius)
        transform.rotate_structure(it, 0.0, angle*i, 0.0)

        # if not first structure, merge with previous walls
        if i == 0:
            st = it
            # find previously marked atoms for bonding
            for atom in st.atom:
                if atom.element == "Xe": # atoms were changed to Xe in wallDev so they could be identified here
                    connectme.append(atom.index)
                    atom.element = "O" # change it back to O now that index is recorded
        if i > 0:
            st = st.merge(it)
            # find previously marked atoms for bonding
            for atom in st.atom:
                if atom.element == "Xe": # atoms were changed to Xe in wallDev so they could be identified here
                    connectme.append(atom.index)
                    atom.element = "O" # change it back to O now that inex is recorded

    # mark the linker atoms needed to connect the walls
    linkme1=[]
    linkme2=[]
    for res in st.residue:
        # put some of these atoms in linkme1
        if res.resnum <= (1000 + iterator + (resnum_end * 10)) and res.resnum >= 1000 + (resnum_end * 10):
            enum = res.resnum - 1000 - (resnum_end * 10)
            for atom in res.atom:
                if atom.pdbname.strip() == linker_end:
                    linkme1.append(atom.index)
        # put the rest of these atoms in linkme2
        if res.resnum < (4000 + iterator + (resnum_start * 10)) and res.resnum >= 4000 + (resnum_start * 10):
            enum = res.resnum - 4000 - (resnum_start * 10)
            for atom in res.atom:
                if atom.pdbname.strip() == linker_start:
                    linkme2.append(atom.index)

    # combine the items of linkme1 and linkme2 in a new list called linkme, such that
    # the atom at each index of linkme is the linker atom that should be bonded
    # to the DNA atom at the same index of connectme
    if len(wall_list) < 3:
        print(f"You have requested {len(wall_list)} walls, but that doesn't really make sense,\
                so this program doesn't know how to handle it. Exiting now")
        raise SystemExit()
    elif len(wall_list) > 8:
        print(f"You have requested {len(wall_list)} walls, but this program has only been tested\
                on cages of up to and including 8 walls. Depending on the structure of the individual walls,\
                the assembly area may become increasingly crowded with increasing wall size, leading to errors.")
    linkme = []
    wall_list_length = len(wall_list)
    for i in range(2,wall_list_length+1):
        linkme.append(linkme1[wall_list_length-i])
        linkme.append(linkme2[wall_list_length-i])
    linkme.append(linkme1[wall_list_length-1])
    linkme.append(linkme2[wall_list_length-1])
    '''
    if len(wall_list) == 4:
        linkme = [linkme1[2],linkme2[2],linkme1[1],linkme2[1],linkme1[0],linkme2[0],linkme1[3],linkme2[3]]
    elif len(wall_list) == 3:
        linkme = [linkme1[1],linkme2[1],linkme1[0],linkme2[0],linkme1[2],linkme2[2]]
    elif len(wall_list) == 5:
        linkme = [linkme1[3],linkme2[3],linkme1[2],linkme2[2],linkme1[1],linkme2[1],linkme1[0],linkme2[0],linkme1[4],linkme2[4]]
    elif len(wall_list) == 6:
        linkme = [linkme1[4],linkme2[4],linkme1[3],linkme2[3],linkme1[2],linkme2[2],linkme1[1],linkme2[1],linkme1[0],linkme2[0],linkme1[5],linkme2[5]]
    elif len(wall_list) == 7:
        linkme = [linkme1[5],linkme2[5],linkme1[4],linkme2[4],linkme1[3],linkme2[3],linkme1[2],linkme2[2],linkme1[1],linkme2[1],linkme1[0],linkme2[0],linkme1[6],linkme2[6]]
    else: # invalid number of walls
        print(f"This program has only been tested on cages with 3-8 walls, but you have requested {len(wall_list)} walls\n\
                If you really want more walls, you can go to the point of this message in the code and straightforwardly
                add an additional elif statement following the pattern ")
        raise SystemExit()
    '''

    # connect the walls by bonding the linker atoms to the DNA atoms
    # the atoms to bond will depend on the number of walls in the cage
    st.addBonds([(connectme[i], linkme[i], 1) for i in range(iterator*2)])

    # remove and add hydrogens in order to get correct formal charges
    delete_hydrogens(st)
    add_hydrogens(st)

    return st

##------------------------------------------------------------------------------------------------------------
# function to add coordinates according to magnitude
def add_coords(coord1, coord2):
    if coord1 > 1 and coord2 > 1:
        return coord1 + coord2
    elif coord1 > 1 and coord2 < 1:
        return coord1 - coord2
    elif coord1 < 1 and coord2 > 1:
        return coord1 - coord2
    elif coord1 < 1 and coord2 < 1:
        return coord1 + coord2

##------------------------------------------------------------------------------------------------------------
## function to bond junction blocks
def junction_builder(block_list, edge, DNA_seq):
    ## constants
    iterator = len(block_list)
    sclr = 1 + 0.1*iterator
    angle = 2*numpy.pi/iterator
    radius = len(DNA_seq[0][0]) + (2 * edge)/tan(angle/2)

    ## move blocks into position
    for i in range(iterator):
        it = block_list[i]
        #newpos = (radius, 0.0, 0.0, 0.0)
        #(x, y, z, q) = newpos - transform.get_centroid(it)

        transform.translate_structure(it, radius * sclr, 0, 0)
        transform.rotate_structure(it, 0.0, 0.0, angle*i)

        if i == 0:
            st = it
        if i > 0:
            st = st.merge(it)

    ## bond blocks together
    connectme_top = []
    connectme_bottom = []
    deleteme = []

    for i in range(iterator):
        ires = i * 100
        for res in st.residue:
            if res.resnum == 1000 + ires + 1:
                for atom in res.atom:
                    if atom.pdbname.strip() == "O5'":
                        connectme_top.append(atom.index)
            if res.resnum == 2000 + ires + len(DNA_seq[0][0])//2+1:
                for atom in res.atom:
                    if atom.pdbname.strip() == 'P':
                        connectme_bottom.append(atom.index)

    # print(connectme_top)
    # print(connectme_bottom) # empty
    # print(deleteme) # empty

    # print('iterator' + str(range(iterator)))
    for i in range(iterator):
        if i < (iterator-1):
            st.addBond(connectme_top[i], connectme_bottom[i + 1], 1)
        else:
            st.addBond(connectme_top[i], connectme_bottom[0], 1)

    st.deleteAtoms(deleteme)
    add_hydrogens(st)

    return st
