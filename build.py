import os
import sys
import numpy as np
import wallDev as wd
import cage_dev as dev
from schrodinger import structure
from inputDev import fileReader, dnaStruct
from schrodinger.structutils import transform, analyze, minimize

# Main script used to build DNA-Cages. Reads input file, builds cage, modifies cage as necessary then exports ##
##------------------------------------------------------------------------------------------------------------
# function used to minimize and export cage
def min_cage(st):
    min = minimize.Minimizer(struct=st)
    min.updateCoordinates(st)
    min.minimize()
    st_min = min.getStructure()
    return st_min

##------------------------------------------------------------------------------------------------------------
# function to place counter ions
def place_mg(struct, dist, num_walls, mg_conc, resnum=1500, buff_dist=10):
#    k=st_min.formal_charge
    #random.seed(17)
    # initiliaze constants / Mg
    mg = next(structure.StructureReader("Mg.pdb"))
    avocado = 6.0221*10**(23)
    angle = 2*np.pi/num_walls
    radius = dist/(2*np.tan(angle/2))
    box_size = np.power(2*(dist/(2*np.sin(angle/2))+buff_dist), 3)
    box_size_L = box_size*10**(-27)
    atoms = int(np.around(mg_conc/1000*avocado*box_size_L))

    # set coordinates
    x = np.random.permutation(np.linspace(-radius/2,radius/2,num=atoms))
    y = np.random.permutation(np.copy(x))
    z = np.random.permutation(np.copy(x))

    # place Magnesiums
    for i in range(atoms):
        mg.atom[1].x = x[i]
        mg.atom[1].y = y[i]
        mg.atom[1].z = z[i]
        mg.atom[1].resnum = i + resnum
        if i == 0:
            maggedst=struct.merge(mg)
        else:
            maggedst=maggedst.merge(mg)
        #maggedst.write("mag_me.mae")
        return maggedst

## INPUT ##
def main(input_file=None):
    ##------------------------------------------------------------------------------------------------------------
    # read input file and construct list of dnaStructs
    if input_file == None:
        input_file = sys.argv[1]
    struct = fileReader(input_file)

    ## CAGE CREATION ##
    ##------------------------------------------------------------------------------------------------------------
    # create linker object to measure size
    init_link = next(structure.StructureReader(struct.linker))
    linker_size = init_link.measure(struct.linkerS, struct.linkerE)
    # measurements used to build walls
    wallList = []
    ipillar = 1
    width = 3.7 * len(struct.dnaSeqs[0][1]) + linker_size
    width2 = 3.7 * len(struct.dnaSeqs[0][2]) + linker_size
    npillar = len(struct.dnaSeqs) # number of walls to make

    print("Building Pillars")
    # iterate through pillars and build each wall
    for ipillar in range(npillar-1, -1, -1):
        print("pillar: " + str(ipillar))
        # make wall and append to wall list
        wall = wd.buildWall(struct, ipillar, npillar, width, width2)
        #wall.write("wall.pdb")
        wallList.append(wall)

    # build cage
    q_chain_length1 = len(struct.dnaSeqs[0][1])
    q_chain_length2 = len(struct.dnaSeqs[0][2])
    q_chain_length3 = len(struct.dnaSeqs[0][3])
    q_chain_length = [q_chain_length1, q_chain_length2, q_chain_length3]
    print("Assembling pillars into cage")
    st = dev.cage_builder(wallList, width, q_chain_length, struct.linkerSName, struct.linkerEName, struct.linkerSResnum, struct.linkerEResnum)

    # rename chain residues
    for chain in st.chain:
        i = 0
        for residue in chain.residue:
            residue.resnum = i
            i += 1

    ## CAGE MODIFICATION AND EXPORT ##
    ##------------------------------------------------------------------------------------------------------------
    # create save path
    savepath = 'models/cages/'  ## path to output
    outfolder = struct.name + '/'
    #print(savepath + outfolder + struct.name + ".pdb")
    if not os.path.isdir(savepath + outfolder):
        os.makedirs(savepath + outfolder)

    # add magnesium and save cage
    if struct.mg:
        st_mg = place_mg(st, width, npillar, 15) # can manually change Mg concentration, currently at 15 mM (I think)
        st_mg.write(savepath + outfolder +  struct.name + "_Mg_full.pdb")
        st = st_mg

    # minimize and save cage
    if struct.min:
        st_min = min_cage(st)
        st_min.write(savepath + outfolder + struct.name + "_min.pdb")
        st = st_min

    st.write(savepath + outfolder + struct.name + "_full.pdb")

    # add functionals to cage
    #if struct.add:
    #    st = dev.add_func(st, struct, npillar)
    #else:
    #    # mark all comp atoms for deletion
    #    deleteme = []
    #    for atom in st.atom:
    #        if atom.chain >= chr(65 + npillar) and atom.chain < chr(97): # see if this includes linkers
    #            deleteme.append(atom.index)
    #    st.deleteAtoms(deleteme)

    ## save final version
    #to_save = savepath + outfolder + struct.name + "_final.pdb"
    #print(to_save)
    #st.write(to_save)
    #print("Cage saved")

    # return st
    #return [to_save, struct.name]

main()
