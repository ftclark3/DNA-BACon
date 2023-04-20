import random
import sys, numpy
import build_dna2b
from math import tan
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from sklearn.externals import joblib
# functions to build cage and its' additions

##------------------------------------------------------------------------------------------------------------
## function to bond parts of cage
def cage_builder(wall_list, edge, chain_length, linker_start, linker_end, resnum_start, resnum_end):
    a = 0
    b = 0

    iterator = len(wall_list)
    angle = 2*numpy.pi/iterator
    radius = edge/(2*numpy.tan(angle/2))
    for i in range(iterator):

        it = wall_list[i]

        newpos = (radius, 0.0, 0.0, 0.0)
        (x, y, z, q) = newpos - transform.get_centroid(it)

        transform.translate_structure(it, 0, 0, radius)
        transform.rotate_structure(it, 0.0, angle*i, 0.0)

        if i == 0:
            st = it
        if i > 0:
            st = st.merge(it)

    i=0
    linkme1=[]
    linkme2=[]
    #linkme=iterator*2*[None]
    connectme=iterator*2*[None]
    #connectme1 = []
    #connectme2 = []
    deleteme=[]
    for res in st.residue:
        if res.chain < chr(65+iterator) and res.chain >= chr(64):
            enum = ord(res.chain) - 65
            if res.resnum == chain_length[1]/2:
                for atom in res.atom:
                    if atom.pdbname.strip() == "O3'":
                        #print('o3 firing')
                        connectme[enum] = atom.index
                        #connectme1.append(atom.index)
            if res.resnum == chain_length[1]*2 + chain_length[2] + chain_length[0] - chain_length[1]/2 + 1:
                for atom in res.atom:
                    if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                        deleteme.append(atom.index)
                    if atom.pdbname.strip()=="O5'":
                        #print('o5 firing')
                        connectme[enum + iterator] = atom.index
                        #connectme2.append(atom.index)

        if res.resnum <= (1000 + iterator + (resnum_end * 10)) and res.resnum >= 1000 + (resnum_end * 10):
            enum = res.resnum - 1000 - (resnum_end * 10)
            for atom in res.atom:
                if atom.pdbname.strip() == linker_end:
                    linkme1.append(atom.index)

        if res.resnum < (4000 + iterator + (resnum_start * 10)) and res.resnum >= 4000 + (resnum_start * 10):
            enum = res.resnum - 4000 - (resnum_start * 10)
            for atom in res.atom:
                if atom.pdbname.strip() == linker_start:
                    linkme2.append(atom.index)



    linkme = linkme1 + linkme2

    st.addBonds([(connectme[i], linkme[i], 1) for i in range(iterator*2)])
#    for i in range(iterator):
#	st.addBond(connectme[i],linkme[i],1)
#	print(i)

    st.deleteAtoms(deleteme)
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
# function to add functional additions to cage
def add_func(st, struct, npillar):
    keepme = []
    keepme_resnum = []
    coords_final = []
    linkchain = []
    # initially mark all comp chains for deletion
    deleteme = []
    for atom in st.atom:
        if atom.chain >= chr(65 + npillar) and atom.chain < chr(97): # see if this includes linkers
            deleteme.append(atom.index)
    i = 0
    # iterate through chains in order to find complimentary residue to keep
    for chain in st.chain:
        # check if chain is a complimentary chain
        if chain.name > chr(65 + npillar - 1) and chain.name < chr(97):
            chain_seq = ""
            for residue in chain.residue:
                keepme_resnum.append(residue.resnum)
                # add nucleotide to chain's seq
                if "D" in residue.pdbres:
                    chain_seq += residue.pdbres.lstrip(" D").strip(" ")
            # check if chain_seq matches any additions
            for addn in struct.addSeqs:
                if addn in chain_seq:
                    # unmark residue for deletion
                    start_resnum = keepme_resnum[chain_seq.find(addn)]
                    end_resnum = keepme_resnum[len(addn) + start_resnum]
                    for residue in chain.residue:
                        if (residue.resnum >= start_resnum) & (residue.resnum <= end_resnum):
                            for atom in residue.atom:
                                keepme.append(atom.index)
                        # remember chain name for later
                        if residue.resnum == start_resnum - 1:
                            for atom in residue.atom:
                                if atom.pdbname.strip() == "O3'":
                                    coords_final.append(atom.xyz)
                                    deleteme.append(atom.index)

                        if residue.resnum == start_resnum + 1:
                            linkchain.append(chain.name)
    # iterate thru all additions
    i = 0
    bindme = []
    for decoration in struct.addStructs:
        # create 5T linker and decoration as structures
        to_add = next(structure.StructureReader(decoration))

        for atom in to_add.atom:
            if atom.index == 1:
                coords_initial = atom.xyz


        # get updated coordinates
        transform.translate_to_origin(to_add)
        #(a, b, c) = add_coords(coords_final[i][0], coords_initial[0]), add_coords(coords_final[i][1], coords_initial[1]), add_coords(coords_final[i][2], coords_initial[2])

        translation_matrix = transform.get_translation_matrix(coords_final[0])
        #print(translation_matrix)
        transform.transform_structure(to_add, translation_matrix)

        # rename chain
        for chain in to_add.chain:
            chain.name = "z"

        st = st.merge(to_add, copy_props=True)

        # bind addition to complimentary strand
        bindme_temp = []
        for chain in st.chain:
            if chain.name == "z":
                for atom in chain.atom:
                    bindme_temp.append(atom.index)
        bindme.append(bindme_temp[0])

        i += 1


        deleteme = [x for x in deleteme if x not in keepme]

    # bind decoration to cube
    linkme = [] # cage
    decoration_start = 1
    # for atom in st.atom:
    #     print(atom.index)
    for chain in st.chain:
        if chain.name in linkchain:
            for residue in chain.residue:
                if residue.resnum == start_resnum:
                    for atom in residue.atom:
                        if atom.pdbname.strip() == "P":
                            linkme.append(atom.index)



    for i in range(len(bindme)):
        st.addBond(bindme[i], linkme[i], 1)
    st.deleteAtoms(deleteme)

    return st
