# import sys
import numpy as np
import build_dna2b
from math import pi
from schrodinger import structure
from schrodinger.structutils import transform, analyze, minimize
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens

# script to create DNA rings
# ------------------------------------------------------------------------------
# function to trim DNA and prepare it for bonding
def trim_DNA(wall_length, st):
    to_delete = []
    for chain in st.chain:
        if chain.name == "A":
            for res in chain.residue:
                # print("res")
                # print(res.resnum)
                if res.resnum == 1: # this residue appears TWICE! Second occurance is the weird HXL one
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3T" or atom.pdbname.strip() == "H3T":
                            to_delete.append(atom.index)
                if res.resnum == (wall_length + 2):
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O5T" or atom.pdbname.strip() == "H5T":
                            to_delete.append(atom.index)
                        if atom.pdbname.strip() == "P" or atom.pdbname.strip() == "O1P":
                            to_delete.append(atom.index)
                        if atom.pdbname.strip() == "O2P" or atom.pdbname.strip() == "O3T":
                            to_delete.append(atom.index)
        if chain.name == "B":
            for res in chain.residue:
                if res.resnum == 1: # this residue appears TWICE! Second occurance is the weird HXL one
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3T" or atom.pdbname.strip() == "H3T":
                            to_delete.append(atom.index)
                if res.resnum == (wall_length + 2):
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O5T" or atom.pdbname.strip() == "H5T":
                            to_delete.append(atom.index)
                        if atom.pdbname.strip() == "P" or atom.pdbname.strip() == "O1P":
                            to_delete.append(atom.index)
                        if atom.pdbname.strip() == "O2P" or atom.pdbname.strip() == "O3T":
                            to_delete.append(atom.index)
    # print(to_delete)
    st.deleteAtoms(to_delete)
    return st

# ------------------------------------------------------------------------------
# function to bond first and last walls
def final_bond(res_length, st):
    # initialize lists
    bindme = []
    linkme = []
    to_delete = []
    # parse through structure and find atoms to bind/delete
    for chain in st.chain:
        if chain.name == "A":
            for res in chain.residue:
                if res.resnum == res_length:
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3'":
                            bindme.append(atom.index)
                if res.resnum == 1:
                    for atom in res.atom:
                        if atom.pdbname.strip() == "P":
                            linkme.append(atom.index)
        if chain.name == "B":
            for res in chain.residue:
                if res.resnum == res_length:
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3'":
                            bindme.append(atom.index)
                if res.resnum == 1:
                    for atom in res.atom:
                        if atom.pdbname.strip() == "P":
                            linkme.append(atom.index)

    # merge structures and bind and delete atoms
    for i in range(len(bindme)):
        st.addBond(bindme[i], linkme[i], 1)
    st.deleteAtoms(to_delete)
    return st

# ------------------------------------------------------------------------------
# function to bond wall to previous wall
def bond_walls(res_num, st):
    # initialize lists
    bindme = []
    linkme = []
    # parse through structure and find atoms to bind/delete
    for chain in st.chain:
        if chain.name == "A":
            for res in chain.residue:
                if res.resnum == (res_num + 1):
                    for atom in res.atom:
                        if atom.pdbname.strip() == "P":
                            bindme.append(atom.index)
                if res.resnum == (res_num):
                    # print("end" + str(res.resnum))
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3'":
                            linkme.append(atom.index)
                            # print("hit4")
        if chain.name == "B":
            for res in chain.residue:
                if res.resnum == (res_num + 1):
                    for atom in res.atom:
                        if atom.pdbname.strip() == "P":
                            bindme.append(atom.index)
                if res.resnum == (res_num):
                    # print("end" + str(res.resnum))
                    for atom in res.atom:
                        if atom.pdbname.strip() == "O3'":
                            # print("hit3")
                            linkme.append(atom.index)
    # print(bindme)
    # print(linkme)

    # bind and return structure
    for i in range(len(bindme)):
        st.addBond(bindme[i], linkme[i], 1)
    # st.deleteAtoms(to_delete)
    return st

# ------------------------------------------------------------------------------
# DEPRECIATED function to rename resids of a structure:
def renumber(start, struct, total_resids):
    # print("start")
    # print(start)
    for chain in struct.chain:
        # name resids of chain A backwards
        if chain.name == "A":
            i = total_resids - start
            print(i)
            for res in chain.residue:
                print(res.inscode)
                res.resnum = i
                i += 1
                print(i)
        else:
            i = start
            for res in chain.residue:
                res.resnum = i
                i += 1
    return [struct, i-1]

# ------------------------------------------------------------------------------
# function used to move and bind alkyl additions to DNA chains
def add_addns(seq, file, indices, st):
    # add one to each index
    # print(indices)
    for i in range(len(indices)):
        indices[i] += 1
    # print(indices)
    linkme = []
    bindme = []
    deleteme = []
    BP_length = 2.5

    # iterate through each residue
    for res in st.residue:
        if res.resnum in indices and res.chain == "A": # check that we have the correct index and chain
            # load in structure to add
            addn = next(structure.StructureReader(file))
            for chain in addn.chain:
                chain.name = "C"

            # delete extra from carbon we will bond to (CAN DELETE IF YOU BUILD WITHOUT THIS H)
            for atom in addn.atom:
                if atom.index == 22: #THIS IS DECYL SPECIFIC SORRY
                    # print(atom.index)
                    delete = [atom.index]
            addn.deleteAtoms(delete)
            # addn.write("test_delete.pdb")

            # move addition along z axis
            addn_mvmt = res.resnum*BP_length
            transform.translate_structure(addn, x=20, z=addn_mvmt) # move out MUST CHANGE X MANUALLY RIGHT NOW



            # merge addn and original structure
            st = st.merge(addn, copy_props=True)


            # for atom in addn.atom:
            #     # print(atom.chain)
            #     # print(atom.index)
            #     # print("\n")
            #     if atom.index == 1 and atom.chain == "C":
            #         linkme.append(atom.index)
            #         # this seems redundent but we'll need it when we merge structures
            #     if atom.index == 22 and atom.chain == "C": # delete extra hydrogen
            #         deleteme.append(atom.index)
            # print(linkme)
    # find oxygen and replace with sulfur
    for res in st.residue:
        if res.resnum in indices and res.chain == "A": # check that we have the correct index and chain
            for atom in res.atom:
                if atom.pdbname.strip() == "OP2":
                    # print("hit")
                    atom.element = "S"
                    atom.pdbname = "S"
                    atom.formal_charge = 0
                    # print(atom.atom_type)
                    bindme.append(atom.index)

    #
    # # add sulfur to first carbon
    # for res in st.residue:
        if res.chain == "C":
            to_add_index = 100000000
            for atom in res.atom:
                if atom.index < to_add_index and atom.element == "C":
                    to_add_index = atom.index
                    linkme.append(atom.index)

            # bind addn to structure
    # print(linkme)
    # print(bindme)
    for i, bind in enumerate(bindme):
        # print(linkme[i], bind)
        st.addBond(linkme[i], bind, 1)

    # st.addBond(linkme, bindme, 1)
    st.deleteAtoms(deleteme)

            # move addn into place
    minimize.minimize_structure(st)


    return st

# ------------------------------------------------------------------------------
# initialize struct
filename = "DEB_ring.pdb"
total_walls = 7
res_num = 0
# init_seq = "TTTTTCACACTTTTTCACACT" # nanopore
init_seq = "TTTTTCACACTTTTTCACACT" # DEB
addn_file = "addn.mae"
# addn_place = [2, 4, 13, 15] # nanopore
addn_place = [6, 8, 17, 19] # DEB
BP_length = 3.7
BP_per_twist = 5
circ = (len(init_seq)+1)*total_walls*BP_length
radius = circ/(2*pi)

# print(len(init_seq))
total_resids = total_walls*(len(init_seq))
# print(total_resids)
# true_total_resids = total_walls*(len(init_seq)+1)
init_struct = build_dna2b.process_sequence(init_seq)
init_struct = add_addns(init_seq, addn_file, addn_place, init_struct)
init_struct = trim_DNA(len(init_seq), init_struct)

# add_hydrogens(init_struct)

# renumber initial struct
for chain in init_struct.chain:
    # CHAIN A IS WHAT WE BIND TO
    if chain.name == "A":
        i = 1
        for res in chain.residue:
            res.resnum = i
            i += 1
    if chain.name == "B":
        i = 1
        for res in chain.residue:
            res.resnum = i
            i += 1

# translate centroid to origin
transform.translate_centroid_to_origin(init_struct)
transform.rotate_structure(init_struct, z_angle=pi) # rotate along axis of DNA


# add addns to initial structure
# init_struct.write("init.pdb")
# exit()

# move and renumber init struct
init_struct.name = "A"
final_struct = init_struct.copy()
transform.translate_structure(final_struct, radius, y=0, z=0) # move out

# ------------------------------------------------------------------------------
# duplicate and move walls
for i in range(1, total_walls):
    # print(i)
    # initalize tmp structure
    outer_angle = (2*pi/total_walls)*i # angle to rotate along circle
    twist_angle = ((len(init_seq))%BP_per_twist)/BP_per_twist*(2*pi) # angle to rotate along axis of DNA
    tmp_struct = init_struct.copy()
    tmp_struct.name = chr(65 + i)

    # renumber residues
    j = total_resids - (i*len(init_seq)) + 1
    for res in tmp_struct.residue:
        for chain in res.chain:
            if chain == "A":
                res.resnum = j
                j += 1
    j = i*len(init_seq) + 1
    for res in tmp_struct.residue:
        for chain in res.chain:
            if chain == "B":
                res.resnum = j
                j += 1

    # rotate and move temp structure into position
    # transform.rotate_structure(tmp_struct, z_angle=twist_angle) # rotate along axis of DNA
    transform.translate_structure(tmp_struct, radius, y=0, z=0) # move out
    transform.rotate_structure(tmp_struct, y_angle=outer_angle) # rotate along the circle
    final_struct = final_struct.merge(tmp_struct, copy_props=True)

# final_struct.write("test_before_binding.pdb")

# ------------------------------------------------------------------------------
# bind walls together
res_num = len(init_seq)
for i in range(1, total_walls):
    final_struct = bond_walls(res_num, final_struct)
    res_num += len(init_seq)
    # print("res")
    # print(res_num)

# bind last wall to first wall
total_resids = total_walls*(len(init_seq))
final_struct = final_bond(total_resids, final_struct)
final_struct.write("models/unmin_" + filename)

# delete_hydrogens(final_struct)
# add_hydrogens(final_struct)
minimize.minimize_structure(final_struct)

final_struct.write("models/min_" + filename)
