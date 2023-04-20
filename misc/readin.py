import schrodinger
import random
import sys, numpy
import build_dna2b
import os
import dev
from math import pi
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from sklearn.externals import joblib

## main function - to call from main
def readin(inputfile, linkerfile):
    ##------------------------------------------------------------------------------------------------------------
    savepath = 'models/cages/'  ## path to output
    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    ## set linker variables
    if linkerfile == 'HEG.pdb':
        HEG = True
        linkerS = 11
        linkerS_name = "C8"
        linkerS_resnum = 1
        linkerE = 17
        linkerE_name = "C12"
        linkerE_resnum = 1
        linker_center = "O1"

    elif linkerfile == '4T.pdb':
        HEG = False
        linkerS = 5
        linkerS_name = "C5*"
        linkerS_resnum = 3
        linkerE = 104
        linkerE_name = "C3*"
        linkerE_resnum = 6
        linker_center = "P*"

    elif linkerfile == 'benzene.pdb':
        HEG = False
        linkerS = 7
        linkerS_name = "C7"
        linkerS_resnum = 1
        linkerE = 13
        linkerE_name = "C13"
        linkerE_resnum = 1
        linker_center = "C5"

    else:
        ## for benzene:
            # 13
            # C13
            # 1
            # 7
            # C7
            # 1
        linkerS = input('linker start atom number: ')
        linkerS_name = str(input('linker start name: '))
        linkerS_resnum = input('linker start resnum: ')
        linkerE = input('linker end atom number: ')
        linkerE_name = str(input('linker end name: '))
        linkerE_resnum = input('linker end resnum :')
        linker_center = str(input('linker center atom: '))

    ## initialize DNA structure
    DNA_seq = []
    wall_list = []

    ## initialize linker
    init_link = next(structure.StructureReader(linkerfile))
    linkerSize = init_link.measure(linkerS,linkerE)
    iterator = 5000

    ##------------------------------------------------------------------------------------------------------------
    ## create DNA strand as structure
    dna_structures = dev.file_reader(inputfile)
    #except:
        #dna_structures = dev.file_reader(input("It appears you did not input a valid file, please type it now: "))

    for dna_struct in dna_structures:
        ## measure number of chains and lengths
        DNA_seq = dna_struct.DNA_seq
        npillar = len(DNA_seq) # number of lines in DNA seq file determines vertices
        chain_length = len(DNA_seq[0][1])
        chain_length2 = len(DNA_seq[0][2])
        chain_length3 = len(DNA_seq[0][3])
        q_chain_length = []
        q_chain_length.append(len(DNA_seq[0][1]))
        q_chain_length.append(len(DNA_seq[0][2]))
        q_chain_length.append(len(DNA_seq[0][3]))
        #print(chain_length2)
        width = 3.7 * chain_length + linkerSize
        width2 = 3.7 * chain_length2 + linkerSize

        ## if specified, initialize addition
        Add = False
        if dna_struct.Add:
            Add = True
            AddSeqs = dna_struct.AddSeqs
            AddStructs = dna_struct.AddStructs
            # to_add_comp = dev.comp_seq(to_add)
            # to_del_comp = [0,0]


        ## build cage by iterating over pillars
        for ipillar in range(npillar-1, -1, -1):
            #print('pillar: ', ipillar)
            chain_name = chr(65+ipillar)
            #print('chain name: ', chain_name)


            ##------------------------------------------------------------------------------------------------------------
            ## move first column
            ct_1 = build_dna2b.process_sequence(DNA_seq[ipillar][1])

            # if addon necc
            #ct_1.write("ct_1_0.pdb")

            atom_index_list = range(1, ct_1.atom_total+1)
            #transform_matrix = transform.get_rotation_matrix_from_eulers(pi, 0, 0)
            #transform.transform_structure(ct_1, transform_matrix, atom_index_list)
            #schrodinger.structutils.transform.rotate_structure(ct_1, x_angle=0, y_angle=0.5*pi, z_angle=0, rot_center=None)


            ## rotate into posn
            transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([1, 0, 0]))
            transform.transform_structure(ct_1, transform_matrix, atom_index_list)


            ## get centroid move
            centroid = transform.get_centroid(ct_1, atom_list=[atom.index for atom in ct_1.atom if atom.element=="P"])
            newpos = numpy.array([0.0, 0.5*width2, 0.0, 0.0])
            (x, y, z, q) = newpos - centroid
            transform.translate_structure(ct_1, x, y, z, atom_index_list)

            for chain in ct_1.chain:
                if chain.name == "A":
                    chain.name = chain_name
                    #print('hit ' + chain.name)
                else:
                    # check to see if comp chain to addn
                    chain.name = chr(65+(npillar*2)+ipillar)
                    #print('else hit ' + chain.name)
            ## check if add
            # if Add:
            #     if to_add_comp in DNA_seq[ipillar][1]:
            #         for chain in ct_1,chain:
            #             chain.name = "+"


            for res in ct_1.residue:
                print("uh")
                print(DNA_seq[ipillar][0])
                res.resnum = len(DNA_seq[ipillar][0])+res.resnum-1

                #ct_1.write("ct_1.pdb")

            ##------------------------------------------------------------------------------------------------------------
            ## move third column
            ct_3 = build_dna2b.process_sequence(DNA_seq[ipillar][3])

            transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([-1, 0, 0]))
            atom_index_list = range(1, ct_3.atom_total+1)
            transform.transform_structure(ct_3, transform_matrix, atom_index_list)

            centroid = transform.get_centroid(ct_3, atom_list=[atom.index for atom in ct_3.atom if atom.element=="P"])
            #print(centroid)

            newpos = numpy.array([0.0, -0.5*width2, 0.0, 0.0])
            #print(newpos)

            (x, y, z, q) = newpos - centroid
            #print(x)
            transform.translate_structure(ct_3, x, y, z, atom_index_list)

            for chain in ct_3.chain:
                if chain.name == "A":
                    chain.name = chain_name
                    #print('hit ', chain.name)
                else:
                    #chain.name = chr(65+ipillar+2*npillar)
                    chain.name = chr(65+npillar+ipillar)




            for res in ct_3.residue:
                res.resnum = len(DNA_seq[ipillar][0]+DNA_seq[ipillar][1]+DNA_seq[ipillar][2])+res.resnum-1

                #ct_3.write("ct_3.pdb")

            ##------------------------------------------------------------------------------------------------------------
            ## move second/middle column
            ct_2 = build_dna2b.process_sequence(DNA_seq[ipillar][2])


            transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, -1, 0]))
            atom_index_list = range(1, ct_2.atom_total+1)
            transform.transform_structure(ct_2, transform_matrix, atom_index_list)
            schrodinger.structutils.transform.rotate_structure(ct_2, x_angle=0, y_angle=0.5*pi, z_angle=0, rot_center=None)


            centroid = transform.get_centroid(ct_2, atom_list=[atom.index for atom in ct_2.atom if atom.element=="P"])
            #print(centroid)

            newpos = numpy.array([ 0.5*width, 0.0, 0.0, 0.0])
            #print(newpos)

            (x, y, z, q) = newpos - centroid
            #print(x)
            transform.translate_structure(ct_2, x, y, z, atom_index_list)

            for chain in ct_2.chain:
                if chain.name == "A":
                    chain.name = chain_name
                    for res in chain.residue:
                        res.resnum = len(DNA_seq[ipillar][0]+DNA_seq[ipillar][1])+res.resnum-1
                else:
                    if ipillar > 0:
                        chain.name = chr(65+ipillar-1)
                    else:
                        chain.name = chr(65+npillar-1)

                    n1 = len(DNA_seq[ipillar][0])
                    n2 = len(DNA_seq[ipillar][0]+DNA_seq[ipillar][1]+DNA_seq[ipillar][2]+DNA_seq[ipillar][3])

                    for res in chain.residue:
                        if res.resnum < n1+2:
                            res.resnum = n2+res.resnum-1
                        else:
                            res.resnum = res.resnum-n1-1
                        if res.resnum == 1:
                            atom1 = [atom for atom in res.atom if atom.pdbname.strip() == "P"]
                            atom2 = [atom for atom in atom1[0].bonded_atoms if atom.resnum >1]
                            #print(atom1[0].resnum, atom1[0].pdbname, atom2[0].resnum, atom2[0].pdbname)
                            ct_2.deleteBond(atom1[0], atom2[0])
                            add_hydrogens(ct_2)
                            atom3 = [atom for atom in atom1[0].bonded_atoms if atom.element == "H"]
                            atom3[0].element = "O"
                            atom3[0].pdbname = " OP3"
                            atom3[0].formal_charge=-1


        #    ct_2.write("ct_2.pdb")


            st = ct_1.merge(ct_2,copy_props=True)
            st = st.merge(ct_3,copy_props=True)
        #    st.write("wall_no_linker.pdb")

            ##------------------------------------------------------------------------------------------------------------
            ## linker constants
            end = linkerE_resnum * 10 + ipillar
            start = linkerS_resnum * 10 + ipillar
            HEGaddn = linkerE_resnum * 10 + ipillar

            ## create initial linker object and atom index list
            link = init_link.copy()
            atom_index_list = range(1, link.atom_total+1)

            ## name initial residues
            for residue in link.residue:
                ## if using HEG linker
                if linkerE_resnum == linkerS_resnum:
                    residue.resnum = 2000 + HEGaddn
                    chain.name = '5'     # why am I using 5...?

                ## if on last linker resnum
                elif residue.resnum == linkerE_resnum:
                    residue.resnum = 2000 + end
                    chain.name = '5'

                ## if on first linker resnum
                elif residue.resnum == linkerS_resnum:
                    residue.resnum = 2000 + start
                    chain.name = '5'

                ## if on middle linker resnum
                else:
                    residue.resnum = (residue.resnum) + 2000 + (ipillar *100)

            ##------------------------------------------------------------------------------------------------------------
            ## move four linkers into place
            #print(schrodinger.structure.Structure.__file__)

            ## link up four columns to create a complete side (pillar)
            for atom in link.atom:
                #print(atom.pdbname.strip())
                if atom.pdbname.strip() == linker_center:
                        ## initial linker (top right)
                        schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=-0.5*pi, z_angle=0, rot_center=None)
                        (a,b,c) = numpy.array([0.45*width, 0.52*width2, 0.2*width]) - atom.xyz  # new linker position
                        transform.translate_structure(link, a, b, c, atom_index_list)   # move linker

                        st = st.merge(link)

                        ## rename second linker residues (bottom right)
                        for residue in link.residue:
                            #    print(residue.resnum)
                                if residue.resnum == 2000 + HEGaddn:
                                    residue.resnum = 3000 + HEGaddn
                                elif residue.resnum == 2000 + end:
                                    residue.resnum = 3000 + end
                                elif residue.resnum == 2000 + start:
                                    residue.resnum = 3000 + start
                                else:
                                    residue.resnum = residue.resnum + 1000 + ipillar

                        #transform_matrix = transform.get_alignment_matrix(numpy.array([0, 1, 0]), numpy.array([1, 0, 0]))
                        #transform.transform_structure(link, transform_matrix, atom_index_list)
                        schrodinger.structutils.transform.rotate_structure(link, x_angle=0.25*pi, y_angle=0, z_angle=0, rot_center=None)
                        (a,b,c) = numpy.array([0.45*width, -0.5*width2, 0.15*width]) - atom.xyz
                        transform.translate_structure(link, a, b, c, atom_index_list)   # move linker

                        st = st.merge(link)

                        ## rename third linker residues (top left)
                        for residue in link.residue:
                                if residue.resnum == 3000 + HEGaddn:
                                    residue.resnum = 1000 + HEGaddn
                                elif residue.resnum == 3000 + end:
                                    residue.resnum = 1000 + end
                                elif residue.resnum == 3000 + start:
                                    residue.resnum = 1000 + start
                                else:
                                    residue.resnum = (residue.resnum) + 2000 + ipillar

                        #transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, 0, -1]))
                        #transform.transform_structure(link, transform_matrix, atom_index_list)
                        if linkerfile == "4T.pdb":
                                #print('4 what')
                                schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=0, z_angle=pi/2, rot_center=None)
                                (a,b,c) = numpy.array([-0.5*width, 0.5*width2, 0.25*width]) - atom.xyz
                                transform.translate_structure(link, a, b, c, atom_index_list)
                        else:
                                schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=0, z_angle=pi, rot_center=None)
                                (a,b,c) = numpy.array([-0.45*width, 0.5*width2, 0.15*width]) - atom.xyz
                                transform.translate_structure(link, a, b, c, atom_index_list)
                        st = st.merge(link)

                        ## rename fourth linker residues (bottom left)
                        for residue in link.residue:
                                if residue.resnum == 1000 + HEGaddn:
                                    residue.resnum = 4000 + HEGaddn
                                elif residue.resnum == 1000 + end:
                                    residue.resnum = 4000 + end
                                elif residue.resnum == 1000 + start:
                                    residue.resnum = 4000 + start
                                else:
                                    residue.resnum = (residue.resnum) + 3000 + ipillar

                        if linkerfile == "4T.pdb":
                                #print('4 what')
                                schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=0, z_angle=pi/2, rot_center=None)
                                (a,b,c) = numpy.array([-0.4*width, -0.5*width2, 0.25*width]) - atom.xyz
                                transform.translate_structure(link, a, b, c, atom_index_list)
                        else:
                                (a,b,c) = numpy.array([-0.4*width, -0.5*width2, 0.15*width]) - atom.xyz
                                transform.translate_structure(link, a, b, c, atom_index_list)

                        st = st.merge(link)
                        #st.write("wall.pdb")


            ##------------------------------------------------------------------------------------------------------------
            ## initialize bonding lists
            deleteme = []
            connectme = []
            linkmeS = []
            linkmeE = []
            addme = []
            connectme_add = []
            i = 1

            ## delete atoms (ask)
            for atom in st.atom:
                if atom.pdbres.strip() == "POT" or atom.pdbres.strip() == "HXL":
                        if atom.chain < chr(65+npillar+ipillar):
                                deleteme.append(atom.index)
            st.deleteAtoms(deleteme)
            deleteme = []

            ## bond or delete residues
            for residue in st.residue:
                ## if not a linker residue, delete extraneous atoms
                if residue.chain == chain_name and residue.chain != '5':
                    if i == 1:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                            deleteme.append(atom.index)
                                    if atom.pdbname.strip()=="O5'":
                                            connectme.append(atom.index)

                    if i == chain_length:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="O3'":
                                            connectme.append(atom.index)

                    if i == 1 + chain_length:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                            deleteme.append(atom.index)
                                    if atom.pdbname.strip()=="O5'":
                                            connectme.append(atom.index)

                    if i == chain_length + chain_length2:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="O3'":
                                            connectme.append(atom.index)

                    if i == 1 + chain_length + chain_length2:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                            deleteme.append(atom.index)
                                    if atom.pdbname.strip()=="O5'":
                                            connectme.append(atom.index)

                    if i == chain_length + chain_length2 + chain_length3:
                            for atom in residue.atom:
                                    if atom.pdbname.strip()=="O3'":
                                            connectme.append(atom.index)

                    i = i + 1

                ## link atoms from first linker to columns
                if residue.resnum == 1000 + ipillar or residue.resnum == 1000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkmeE.append(atom.index)


                if residue.resnum == 1000 + ipillar or residue.resnum == 1000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkmeS.append(atom.index)

                ## link atoms from second linker to columns
                if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkmeE.append(atom.index)

                if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkmeS.append(atom.index)

                ## link atoms from third linker to columns
                if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkmeE.append(atom.index)

                if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkmeS.append(atom.index)


                ## link atoms from fourth linker to columns
                if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkmeE.append(atom.index)


                if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkmeS.append(atom.index)



            ## add bonds dont use linkme 1/6
            ## should be c8 to o5
            st.addBonds([(linkmeS[2], connectme[4], 1), # c12 o5
                         (linkmeS[1], connectme[2], 1), # c8 o3
                         (linkmeE[3], connectme[5], 1), # c12 o5
                         (linkmeS[0], connectme[0], 1), # c8 o3
                         (linkmeE[1], connectme[1], 1), # c12 o5
                         (linkmeE[2], connectme[3], 1)
                         ]) # c8 o3

        #    st.write('wall_linkers.pdb')

            ## delete extraneous atoms
            st.deleteAtoms(deleteme)
            connectme = []
            deleteme = []
            linkme = []

            ## add wall to entire structure
            wall_list.append(st)
            #st.write("pillar"+chain_name+".pdb")
            del st
            #sys.exit(1)

        ##------------------------------------------------------------------------------------------------------------
        print('start')
        #for residue in wall_list[0].residue:
         #   print(residue.resnum)
        dist = width #change to increase corner distances

        ## create cage as object
        st = dev.cage_builder(wall_list, dist, q_chain_length, linkerS_name, linkerE_name, linkerS_resnum, linkerE_resnum)

        outfolder = str(npillar) + '_' + dna_struct.name + '/'
        print(savepath + outfolder + dna_struct.name + ".pdb")
        if not os.path.isdir(savepath + outfolder):
            os.makedirs(savepath + outfolder)

        # rename chains
        for chain in st.chain:
            i = 0
            for residue in chain.residue:
                residue.resnum = i
                i += 1

        deleteme = []
        # mark all comp atoms for deletion
        for atom in st.atom:
            if atom.chain >= chr(65 + npillar): # see if this includes linkers
                deleteme.append(atom.index)


        # check if user wants addition(s)
        if Add:
            keepme = []
            keepme_resnum = []
            coords_final = []
            linkchain = []
            i = 0
            # iterate through chains in order to find complimentary residue to keep
            for chain in st.chain:
                # check if chain is a complimentary chain
                if chain.name > chr(65 + npillar -1):
                #    print(chain.name)
                    chain_seq = ""
                #    print("hit")
                    for residue in chain.residue:
                        keepme_resnum.append(residue.resnum)
                        if "D" in residue.pdbres:
                            chain_seq += residue.pdbres.lstrip(" D").strip(" ")
                    # check if chain_seq matches any additions
                    for addn in AddSeqs:
                        if addn in chain_seq:
                            #print("its in")
                            # unmark residue for deletion
                            start_resnum = keepme_resnum[chain_seq.find(addn)]
                            end_resnum = keepme_resnum[len(addn) + start_resnum]
                            for residue in chain.residue:
                                if (residue.resnum >= start_resnum) & (residue.resnum <= end_resnum):
                                #    print("resnum" + str(residue.resnum))
                                    for atom in residue.atom:
                                        keepme.append(atom.index)\
                                # remember chain name for later
                                if residue.resnum == start_resnum:
                                    for atom in residue.atom:
                                        #print(atom.pdbname.strip())
                                        if atom.pdbname.strip() == "O3T":
                                            coords_final.append(atom.xyz)
                                            deleteme.append(atom.index)

                                if residue.resnum == start_resnum + 1:
                                    linkchain.append(chain.name)
            # iterate thru all additions
            i = 0
            for decoration in AddStructs:
                # create 5T linker and decoration as structures
                to_add = next(structure.StructureReader(decoration))

                # bind addition to complimentary strand
                bindme = []
                for atom in to_add.atom:
                    if atom.index == 1:
                        bindme.append(atom.index)
                        coords_initial = atom.xyz

                # move addition into place
                for atom in to_add.atom:
                    if atom.index == 0:
                        coords_initial = atom.xyz
                # get updated coordinates
                (a, b, c) = dev.add_coords(coords_final[i][0], coords_initial[0]), dev.add_coords(coords_final[i][1], coords_initial[1]), dev.add_coords(coords_final[i][2], coords_initial[2])
                transform.translate_structure(to_add, a, b, c)

                # rename chain
                for chain in to_add.chain:
                    chain.name = "+"

                st = st.merge(to_add)
                i += 1


        deleteme = [x for x in deleteme if x not in keepme]
        st.deleteAtoms(deleteme)

        # bind decoration to cube
        linkme = [] # cage
        bindme = [] # decoration
        decoration_start = 1

        for chain in st.chain:
            if chain.name in linkchain:
                for residue in chain.residue:
                    if residue.resnum == start_resnum + 1:
                        for atom in residue.atom:
                            if atom.pdbname.strip() == "P":
                                linkme.append(atom.index)
            if chain.name == "+":
                for residue in chain.residue:
                    if residue.resnum == decoration_start:
                        for atom in residue.atom:
                            if atom.pdbname.strip() == "C1":
                                bindme.append(atom.index)
        #print(linkme, bindme)

        for i in range(len(bindme)):
            st.addBond(bindme[i], linkme[i], 1)





        ## minimize structure and export
        if dna_struct.Min:
            min = minimize.Minimizer(struct=st)

            min.updateCoordinates(st)
            min.minimize()
            st_min = min.getStructure()
        else:
            st_min=st.copy()

        st_min.write(savepath + outfolder + str(npillar) + "_full_" + dna_struct.name + ".pdb")


        ## add Mg to balance charge
        if dna_struct.Mg:
        #    print("mg running")
            k=st_min.formal_charge
            random.seed(17)
            mg = next(structure.StructureReader("Mg.pdb"))
            #print(i)
            j = 0
            angle = 2*numpy.pi/len(wall_list)
            radius = dist/(2*numpy.tan(angle/2))

            transform.translate_to_origin(st_min)
            for j in range(-st_min.formal_charge//2 + -st_min.formal_charge%2):
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
            maggedst_path = savepath + outfolder + "mg_" + dna_struct.name + ".mae"
            print(maggedst_path)
            maggedst.write(maggedst_path)
            #maggedst.write("~/HDD/Maestro/jobs/" + dna_struct.name + "-setup/" + "mg_" + dna_struct.name + ".mae")


        # return structure
        return [maggedst_path, dna_struct.name]
