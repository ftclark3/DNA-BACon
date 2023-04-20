import schrodinger
import random
import sys, numpy
import build_dna2b
import dev
import os
from seq_writer import gen_DNA
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from sklearn.externals import joblib

##------------------------------------------------------------------------------------------------------------
## read the input file
inputfile = sys.argv[1]
savepath = 'models/junctions/'  ## path to output
if not os.path.isdir(savepath):
    os.makedirs(savepath)

## initialize DNA structure
DNA_seq = []
block_list = []

##------------------------------------------------------------------------------------------------------------
## create DNA strand as structure
try:
    dna_structures = dev.file_reader(inputfile)
except:
    dna_structures = dev.file_reader(input("It appears you did not input a valid file, please type it now: "))

for dna_struct in dna_structures:
    DNA_seq = dna_struct.DNA_seq
    print(DNA_seq[0])
    print('the length is' + str(len(DNA_seq[0][0])//2))

    ## measure length of blocks
    nblock = len(DNA_seq)
    segment_length = len(DNA_seq[0][0])//4 ## length of each segment
    print(segment_length)
    print(len(DNA_seq[0][0])/4)
    width = segment_length * 3


    ## initialize chain naming lists
    ## first is bonding chain
    cname = [[], []]
    ## create chain name lists
    for iblock in range(nblock-1, -1, -1):
        cname[0].append(chr(65 + iblock + nblock))
        cname[1].append(chr(97 + iblock + nblock))

    #print(cname)

    wombat = 0

    ## iterate over blocks to build junction
    for iblock in range(nblock-1, -1, -1):
        print('iblock: ', iblock)
        ires = iblock * 100

        ##------------------------------------------------------------------------------------------------------------
        ## create and move top pillar
        top = build_dna2b.process_sequence(DNA_seq[iblock][0][:(segment_length*2 + 1)])

        ## rotate pillar
        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([1, 0, 0]))
        atom_index_list = range(1, top.atom_total+1)
        transform.transform_structure(top, transform_matrix, atom_index_list)

        ## get centroid of top pillar
        centroid = transform.get_centroid(top, atom_list=[atom.index for atom in top.atom if atom.element=="P"])
        newpos = numpy.array([0.0, -width, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(top, x, y, z, atom_index_list)

        ## renumber residues
        i = 0
        for res in top.residue:
            res.resnum = 1000 + ires + i
            i += 1

        ## rename chain
        ## capital for outer chain / lower for inner
        for chain in top.chain:
            if chain.name == 'A':
                if iblock == (nblock-1):
                    chain.name = cname[0][nblock-1]
                else:
                    chain.name = cname[0][iblock]

            else:
                if iblock == (nblock-1):
                    chain.name = cname[1][nblock-1]
                else:
                    chain.name = cname[1][iblock]

        ##------------------------------------------------------------------------------------------------------------
        ## create and move middle pillar
        mid = build_dna2b.process_sequence(gen_DNA(segment_length))

        ## delete second strand
        deleteme = []
        for chain in mid.chain:
            if chain.name == 'A':
                for atom in chain.atom:
                    deleteme.append(atom.index)
        mid.deleteAtoms(deleteme)

        ## rotate
        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, -1, 0]))
        atom_index_list = range(1, mid.atom_total+1)
        transform.transform_structure(mid, transform_matrix, atom_index_list)

        ## get centroid
        centroid = transform.get_centroid(mid, atom_list=[atom.index for atom in mid.atom if atom.element=="P"])
        atom_index_list = range(1, mid.atom_total+1)
        newpos = numpy.array([0.0, 0.0, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(mid, x, y, z, atom_index_list)

        ## rename residues
        lastres = 0
        for res in mid.residue:
            lastres = res.resnum

        ## rename chain
        for chain in mid.chain:
            chain.name = chr(109+wombat)
        wombat+=1

        ##------------------------------------------------------------------------------------------------------------
        ## create and move bottom pillar
        #if iblock == (nblock-1):
        #    bottom = build_dna2b.process_sequence(DNA_seq[0][0][(segment_length*2 + 1):])
        #else:
        #    bottom = build_dna2b.process_sequence(DNA_seq[iblock + 1][0][(segment_length*2 + 1):])
        bottom = build_dna2b.process_sequence(DNA_seq[iblock][0][(segment_length*2 + 1):])
        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([-1, 0, 0]))
        atom_index_list = range(1, bottom.atom_total+1)
        transform.transform_structure(bottom, transform_matrix, atom_index_list)

        ## get centroid
        centroid = transform.get_centroid(bottom, atom_list=[atom.index for atom in bottom.atom if atom.element=="P"])
        newpos = numpy.array([0.0, width, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(bottom, x, y, z, atom_index_list)

        ## rename residues
        i = 0
        for res in bottom.residue:
            res.resnum = 2000 + ires + i
            i += 1

        ## rename chains
        for chain in bottom.chain:
            if chain.name == 'A':
                if iblock == 0:
                    chain.name = cname[0][nblock-1]
                else:
                    chain.name = cname[0][iblock-1]

            else:
                if iblock == 0:
                    chain.name = cname[1][nblock-1]
                else:
                    chain.name = cname[1][iblock-1]


        ##------------------------------------------------------------------------------------------------------------
        ## merge into one structure
        st = top.merge(mid, copy_props=True)
        st = st.merge(bottom, copy_props=True)
        st.write(str(nblock) + "_block_" + dna_struct.name + ".pdb")
        delete_hydrogens(st)

        ##------------------------------------------------------------------------------------------------------------
        ## bond pillars together
        deleteme = []
        deletebond_top = []
        deletebond_bottom = []
        convertdouble = []
        connectme_mid = []
        connectme_out = []
        check1 = []
        check2 = []

        ## make resnums variables(?)
        for chain in st.chain:
            ## middle
            if chain.name == chr(109+wombat-1):
                #print(lastres)
                for atom in chain.atom:
                    ## delete
                    if atom.resnum == 1:
                        if atom.pdbname.strip() == 'O3T':
                            deleteme.append(atom.index)
                    ## goes to top
                    if atom.resnum == 2:
                        if atom.pdbname.strip() == 'O3T' or atom.pdbname.strip() == 'OP1':
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'OP2':
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == "O5'":
                            connectme_mid.append(atom.index)
                    ## goes to bottom
                    if atom.resnum == lastres:
                        if atom.pdbname.strip() == 'O5T' or atom.pdbname.strip() == 'O1P':
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'O2P':
                            deleteme.append(atom.index)
                    if atom.resnum == lastres - 1:
                        if atom.pdbname.strip() == "O3'":
                            connectme_mid.append(atom.index)

            ## top
            if chain.name == cname[0][iblock] or chain.name == cname[0][iblock-1]:
                for atom in chain.atom:
                    if atom.resnum == 1000 + ires + segment_length:
                        if atom.pdbname.strip() == 'P':
                            connectme_out.append(atom.index)
                            deletebond_top.append(atom.index)
                            check1.append(atom.index)
                    if atom.resnum == 1000 + ires + segment_length - 1:
                        if atom.pdbname.strip() == "O3'":
                            deletebond_top.append(atom.index)
                            check2.append(atom.index)
                    if atom.resnum == 1000 + ires + 1:
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'OP2' or atom.pdbname.strip() == 'OP1' or atom.pdbname.strip() == 'O3T':
                            deleteme.append(atom.index)
                    if atom.resnum == 1000 + ires:
                        if atom.pdbname.strip() == 'O3T':
                            deleteme.append(atom.index)


            ## bottom
            if chain.name == cname[0][iblock] or chain.name == cname[0][iblock-1]:
                for atom in chain.atom:
                    if atom.resnum == 2000 + ires + segment_length:
                        if atom.pdbname.strip() == 'P':
                            connectme_out.append(atom.index)
                            deletebond_bottom.append(atom.index)
                    if atom.resnum == 2000 + ires + segment_length - 1:
                        if atom.pdbname.strip() == "O3'":
                            deletebond_bottom.append(atom.index)
                    if atom.resnum == 2000 + ires + len(DNA_seq[0][0])//2+1:
                        if atom.pdbname.strip() == 'O5T':
                            print(atom.index)
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'O2P':
                            # print('hit' + atom.pdbname.strip())
                            convertdouble.append(atom.index)
                    # if atom.resnum == 2000 + ires + len(DNA_seq[0][0]//2+1:
                    #     if atom.pdbname.strip == "O5'" or atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'OP2' or atom.pdbname.strip() == 'OP1' or atom.pdbname.strip() == 'O3T':


        #    st.write('junct_testing.pdb')
        ## delete double bond and create single bond between P and O2P
        print('checking bonds')
        print(check1, check2)
        print(st.atom[check1[0]].pdbname)
        print(st.atom[check2[0 ]].pdbname)
        # for i in range(len(check1)):
        #     print(st.areBound(check1[i], check2[i]))
        # print(connectme_out, connectme_mid)
        # print('to delete')
        # print(deletebond_top, deletebond_bottom, deleteme)
        try:
            st.deleteBond(convertdouble[0], convertdouble[1])
            st.addBonds([(convertdouble[0], convertdouble[1], 1)])
        except:
            print("could not delete double")
        #st.addBond(convertdouble[0], convertdouble[1], 3)
        #st.addBonds()
        ## bond pillars to form block
        if iblock != 0:
            st.deleteBond(deletebond_top[0], deletebond_top[1])
            st.deleteBond(deletebond_bottom[0], deletebond_bottom[1])
            st.addBonds([(connectme_out[0], connectme_mid[1], 1),
                         (connectme_out[1], connectme_mid[0], 1)])
        else:
            st.deleteBond(deletebond_top[0], deletebond_top[1])
            st.deleteBond(deletebond_bottom[0], deletebond_bottom[1])
            st.addBonds([(connectme_out[0], connectme_mid[0], 1),
                         (connectme_out[1], connectme_mid[1], 1)])

        st.deleteAtoms(deleteme)

        print('checking bonds')
        print(check1, check2)
        for i in range(len(check1)):
            print(st.areBound(check1[i], check2[i]))
        #st.write('nucTest.pdb')
        # print('hellosssss')
        # print(st.atom_total)

        block_list.append(st)
        # print(block_list)

        del st

    ## create junction as object
    outfolder = str(nblock) + '_' + dna_struct.name + '/'
    print(savepath + outfolder + dna_struct.name + ".mae")
    if not os.path.isdir(savepath + outfolder):
        os.makedirs(savepath + outfolder)
    st = dev.junction_builder(block_list, segment_length, DNA_seq)
    print(dna_struct.Min)
    ## minimize structure and export
    if dna_struct.Min:
        min = minimize.Minimizer(struct=st)

        min.updateCoordinates(st)
        min.minimize()
        st_min = min.getStructure()
    else:
        st_min=st.copy()

    st_min.write(savepath + outfolder + str(nblock) + '_' + dna_struct.name + ".mae")

    ## add Mg to balance charge
    if dna_struct.Mg:
        k=st_min.formal_charge
        random.seed(17)
        mg = next(structure.StructureReader("Mg.pdb"))
        #print(i)
        j = 0
        angle = 2*numpy.pi/nblock
        radius = width*10/(2*numpy.tan(angle/2))

        transform.translate_to_origin(st_min)
        for j in range(-st_min.formal_charge//2 + 1):
            i = i + 1
            mg.atom[1].x = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].y = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].z = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].resnum = i
            if j == 0:
                maggedst=st_min.merge(mg)
                #print("hi")
            else:
                maggedst=maggedst.merge(mg)
        maggedst.write(savepath + outfolder + str(nblock) + "_mg_" + dna_struct.name + ".mae")
