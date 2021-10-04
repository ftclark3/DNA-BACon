from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
import build_dev as dev
import wallDev as wd
import numpy as np
import build_dna2b
import blockDev
import toml
import sys
import os
import copy

''' Main script used to build DNA-Cages. Reads input file, builds cage, modifies cage as necessary then exports '''

##------------------------------------------------------------------------------------------------------------
# function used to minimize and export cage
def min_st(st):
    min_st = minimize.Minimizer(struct=st)
    min_st.updateCoordinates(st)
    min_st.minimize()
    return min_st.getStructure()

##------------------------------------------------------------------------------------------------------------
# main junction-building function, assembles and bonds all parts functions
def build_junction(junction_settings):
    ## USER-INPUT ##
    # read junction input file
    struct = junction_settings
    block_list = []

    DNA_seq = struct["sequence"]


    # measure length of blocks
    nblock = len(DNA_seq)
    segment_length = len(DNA_seq[0][0])//4 ## length of each segment
    width = segment_length * 3

    # initialize chain naming lists
    # first is bonding chain
    cname = [[], []]
    # create chain name lists
    for iblock in range(nblock-1, -1, -1):
        cname[0].append(chr(65 + iblock + nblock))
        cname[1].append(chr(97 + iblock + nblock))

    wombat = 0 #*Jonathon loved that*#

    ## iterate over blocks to build junction
    for iblock in range(nblock-1, -1, -1):
        print('iblock: ', iblock)
        ires = iblock * 100

        tmp_block = blockDev.block(iblock, nblock, cname, width, ires, DNA_seq[iblock], segment_length, wombat)
        tmp_block.bond()
        block_list.append(tmp_block.getStructure())


    # build junction from bonds
    st = dev.junction_builder(block_list, segment_length, DNA_seq)

    # rename chains
    for i, molecule in enumerate(st.molecule):
        for residue in molecule.residue:
                residue.chain = chr(65 + i)

    for i, res in enumerate(st.residue):
        res.resnum = i


    # minimize and save junction
    if struct["min"]:
        st_min = min_st(st)
        st = st_min

    print("saving junction")
    st.write("models/junctions/" + str(nblock) + '_' + struct["name"] + ".pdb")

##------------------------------------------------------------------------------------------------------------
# main cage-building function, assembles all parts and calls additional functions
def build_cage(cage_settings,repeat_counter=1,cages=[]):
    ## USER-INPUT ##
    # read cage and linker input files
    struct = cage_settings
    if repeat_counter > 0: # condition for breaking the function repetition
        # create linker object to measure size
        linkerInfo = toml.load("linkers/linkerInformation.toml")
        linker = linkerInfo[struct["linker"]][0]

        init_link = next(structure.StructureReader(linker["file"]))
        linker_size = init_link.measure(linker["linkerS"], linker["linkerE"])

        # measurements used to build walls
        wallList = []
        ipillar = 1
        width = 3.7 * len(struct["sequence"][0][0]) + linker_size
        width2 = 3.7 * len(struct["sequence"][0][1]) + linker_size
        npillar = len(struct["sequence"]) # number of walls to make

        # do this the first time the function is called
        if not cages:
            # catch formatting errors for keep_sequence, add_name, add_num
            try:
                keep = struct["keep_sequence"][0] # complementary DNA sequence to save from deletion
            except IndexError:
                print(f"WARNING: keep_sequence in {sys.argv[1]} is an empty list\n\tall complementary DNA will be deleted")
                keep = None
            except (NameError, KeyError):
                keep = None # allows the user to revert to "default" cage by deleting the variable in the toml
            else:
                repeat_counter = len(struct["keep_sequence"]) # only triggers first time through the loop

            try:
                func = struct["add_name"][0] # get the file name for the first functional group
            except IndexError:
                print(f"Warning: add_name in {sys.argv[1]} is an empty list\n\tno functional groups will be added")
                func = ""
            except (NameError, KeyError):
                func = "" # allows the user to revert to "default cage" by deleting the variable in the toml
            else:
                try:
                    if func[-4:] != ".pdb": # func[-4:] evaluates to empty string in string is empty
                        if func:
                            print("File extension for functional group is required, and .pdb is preferred\n\tothers such as .mae may work but have not been extensively tested")
                            # StructureReader will raise its own exceptions later if needed
                except IndexError:
                    print(f"Error: invalid filename in add_name in {sys.argv[1]}")

            try:
                index = struct["add_num"][0]
            except IndexError:
                if func: # otherwise, index doesn't need to be defined
                    print(f"Error: add_num in {sys.argv[1]} is an empty list\n\tplease add integers until its length is equal to the length of add_name")
                    raise SystemExit()
            except (NameError, KeyError):
                if func:
                    print(f"NameError/KeyError: add_num is not defined in {sys.argv[1]}\n\tplease initialize it as a list of integers or see documentation for more details")
                    raise SystemExit()
                else:
                    pass # allows the user to revert to "default cage" by deleting the variable in the toml
            else:
                if type(index) != type(1):
                    if func:
                        print("Error: add_num must be a list of integers equal in length to add_name\n\tsee documentation for more details")
                        raise SystemExit()

        # do when the function calls itself
        else:
            keep = struct["keep_sequence"][-repeat_counter] # don't need to try-except because repeat_counter is len(keep_sequence)
            try:
                index = struct["add_num"][-repeat_counter]
                if len(struct["add_num"]) == 1:
                    raise IndexError
            except IndexError:
                print(f"IndexError: add_num in {sys.argv[1]} is too short\n\tplease make it equal in length to keep_sequence and add_name or see documentation for more details")
                raise SystemExit()
            if type(index) != type(1):
                print("add_num must be a list of integers equal in length to add_name and keep_sequence\n\tsee documentation for more details")
                raise SystemExit()
            try:
                func = struct["add_name"][-repeat_counter]
                if len(struct["add_name"]) == 1:
                    raise IndexError
            except IndexError:
                print(f"Warning: add_name in {sys.argv[1]} is to short\n\tno more functional groups will be added")

            # pick a new chain name for the additional complementary strands
            alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            chains = []
            global new_chain
            new_chain = []
            for chain in cages[-1].chain:
                chains.append(chain.name)
            for letter in alphabet:
                if letter not in chains:
                    new_chain.append(letter)
                    break

        ## PILLAR CONSTRUCTION ##
        ##------------------------------------------------------------------------------------------------------------
        print("Building Pillars")
        # iterate through pillars and build each wall
        for ipillar in range(npillar-1, -1, -1):
            print("pillar: " + str(ipillar))
            # make wall and append to wall list
            wall = wd.buildWall(struct, ipillar, npillar, width, width2, keep)
            wallList.append(wall)

        # build cage
        q_chain_length1 = len(struct["sequence"][0][0])
        q_chain_length2 = len(struct["sequence"][0][1])
        q_chain_length3 = len(struct["sequence"][0][2])
        q_chain_length = [q_chain_length1, q_chain_length2, q_chain_length3]

        ## CAGE CONSTRUCTION ##
        ##------------------------------------------------------------------------------------------------------------
        print("Assembling pillars into cage")
        st = dev.cage_builder(wallList, width, q_chain_length, linker["linkerSName"], linker["linkerEName"], linker["linkerSResnum"], linker["linkerEResnum"])

        # rename chain residues
        for chain in st.chain:
            i = 0
            for residue in chain.residue:
                residue.resnum = i
                i += 1

        ## CAGE FUNCTIONALIZATION ##
        ##------------------------------------------------------------------------------------------------------------
        if func:
            print("Adding functional group")
            # find the atoms to bond to the functional group
            cage_atoms = []
            for residue in st.residue:
                if residue.pdbres.strip() in ["Xa","Xc","Xg","Xt"]:
                    for atom in residue.atom:
                        if atom.element == "P":
                            for item in atom.bonded_atoms:
                                if item.element == "H":
                                    cage_atoms.append(atom)

            # collect functional group information
            func = structure.StructureReader.read(func)
            for chain in func.chain:
                chain.name = str(repeat_counter)
            cage_center = transform.get_centroid(st)
            x1 = func.atom[index].x
            y1 = func.atom[index].y
            z1 = func.atom[index].z

            # transform functional groups to the correct positions
            for c, item in enumerate(cage_atoms):
                x2 = item.x
                y2 = item.y
                z2 = item.z
                cage_x = x2 - cage_center[0]
                cage_y = y2 - cage_center[1]
                cage_z = z2 - cage_center[2]
                cage_vector = [cage_x,cage_y,cage_z]
                func_center = transform.get_centroid(func)
                func_x = func_center[0] - x1
                func_y = func_center[1] - y1
                func_z = func_center[2] - z1

                ############################################################
                # THIS IS A PLACE WHERE WE COULD DO SOME MATH TO DO A CUSTOM ORIENTATION
                # OF THE FUNCTIONAL GROUP
                # USE [-func_x, -func_y, -func_z] TO PUT THE FUNCTIONAL GROUP INSIDE
                func_vector = [func_x,func_y,func_z]
                ###########################################################

                matrix = transform.get_alignment_matrix(func_vector,cage_vector) # align func_vector with cage_vector
                transform.transform_structure(func,matrix)
                x1 = func.atom[index].x
                y1 = func.atom[index].y
                z1 = func.atom[index].z
                transform.translate_structure(func,1.5+x2-x1,1.5+y2-y1,1.5+z2-z1) # arbitrarily close

                start = len(func.residue) + 1
                for i, residue in enumerate(func.residue):
                    residue.resnum = (c+1) * start + i

                old_element = copy.deepcopy(func.atom[index].element)
                func.atom[index].element = "Rn" # unique identifier

                st = st.merge(func)
                func = structure.StructureReader.read(struct["add_name"][-repeat_counter])

                for chain in func.chain:
                    chain.name = str(repeat_counter)

                x1 = func.atom[index].x
                y1 = func.atom[index].y
                z1 = func.atom[index].z

            # find functional group atoms to bind
            func_atoms = []
            for atom in st.atom:
                if atom.element == "Rn":
                    atom.element = old_element
                    func_atoms.append(atom)

            # 1) find the cage atoms to bond to the functional group
            # (repeating the process done earlier since merging structures
            #  unbinds object references in cage_atoms)
            # 2) delete hydrogens attached to the phosphorus atoms used for bonding
            cage_atoms = []
            deleteme = []
            for residue in st.residue:
                if residue.pdbres.strip() in ["Xa","Xc","Xg","Xt"]:
                    for atom in residue.atom:
                        if atom.element == "P":
                            for item in atom.bonded_atoms:
                                if item.element == "H":
                                    deleteme.append(item.index)
                                    cage_atoms.append(atom)
            st.deleteAtoms(deleteme)

            # bond cage's complementary strand phosphorus atoms with atoms on functional group
            for i, atom in enumerate(cage_atoms):
                st.addBond(atom,func_atoms[i],1)

        # keep track of all cages created
        cages.append(st)
        repeat_counter -= 1
        build_cage(cage_settings,repeat_counter=repeat_counter,cages=cages)

    return cages

## CAGE EXPORT ##
##------------------------------------------------------------------------------------------
# function to merge all cages that have been built and delete duplicate atoms
def assemble_cages(cage_settings):
    struct = cage_settings
    temp = copy.deepcopy(struct["sequence"])
    cages = build_cage(struct)
    # merge all cages and delete unwanted atoms
    st = cages[0]
    for n, cage in enumerate(cages[1:]):
        print("Merging cages")
        deleteme = []
        for residue in cage.residue:
            if residue.pdbres.strip() not in ["Xa","Xc","Xg","Xt"] and residue.chain in "AaBCDEFGHIJKLMNOPQRSTUVWXYZ":
                # residue is not a complementary nucleotide and it is not part of a functional group
                for atom in residue.atom:
                    deleteme.append(atom.index)
        cage.deleteAtoms(deleteme)
        for chain in cage.chain:
            if type(chain.name) == type("string"):
                chain.name = new_chain[n]

    for cage in cages[1:]:
        st = st.merge(cage)

    # make sure that each atom has the correct number of hydrogens
    delete_hydrogens(st)
    add_hydrogens(st)

    # minimize cage
    if struct["min"]:
        print("Minimizing final cage")
        st_min = min_st(st)
        st = st_min

    # save cage
    to_save = 'models/cages/' + struct["name"] + ".pdb"
    st.write(to_save)
    print("Final cage saved: " + to_save)

    # return st
    return [to_save, struct["name"]]

##------------------------------------------------------------------------------------------------------------
# main function, calls assemble_cages and build_junction functions, capable of building many cages
def main():
    to_run = []
    try:
        config = toml.load(sys.argv[1])
    except IndexError:
        config = toml.load(input("You did not pass the input file name as a command-line argument,\nbut you can enter the relative path here:"))
    if "cage" in config:
        for cage_config in config["cage"]:
            to_run.append(assemble_cages(cage_config))
    if "junction" in config:
        for junction_config in config["junction"]:
            to_run.append(build_junction(junction_config))

    return to_run

##------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
