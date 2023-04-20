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
    addendum = init_link.measure(linkerS,linkerE)
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
        width = 3.7 * chain_length + addendum
        width2 = 3.7 * chain_length2 + addendum

        ## if specified, initialize addition
        Add = False
        if dna_struct.Add:
            Add = True
            to_add = next(structure.StructureReader('O3_5T.mae'))

        ## build cage by iterating over pillars
        for ipillar in range(npillar-1, -1, -1):
            print('pillar: ', ipillar)
            chain_name = chr(65+ipillar)
            print('chain name: ', chain_name)


            ##------------------------------------------------------------------------------------------------------------
            ## move first column
            ct_1 = build_dna2b.process_sequence(DNA_seq[ipillar][1])

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
                    #chain.name = chr(65+ipillar+npillar)
                    chain.name = chr(65+(npillar*2)+ipillar)
                    #print('else hit ' + chain.name)

            for res in ct_1.residue:
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
                    print('here: '+chain.name)
                    #chain.name = chr(65+ipillar+2*npillar)
                    chain.name = chr(65+npillar+ipillar)

                    #print('else hit ', chain.name)

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


        #                    ct_2.write("ct_2.pdb")


            st = ct_1.merge(ct_2,copy_props=True)
            st = st.merge(ct_3,copy_props=True)

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
            ## link up four columns to create a side (pillar)
            for atom in link.atom:
                print(atom.pdbname.strip())
                if atom.pdbname.strip() == linker_center:
                        ## initial linker (top right)
                        schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=-0.5*pi, z_angle=0, rot_center=None)
                        (a,b,c) = numpy.array([0.45*width, 0.5*width, 0.15*width]) - atom.xyz  # new linker position
                        transform.translate_structure(link, a, b, c, atom_index_list)   # move linker
                        st = st.merge(link)

    #[0.75*width, 0.5*width2, 10]
                        ## rename second linker residues
                        for residue in link.residue:
                                print(residue.resnum)
                                if residue.resnum == 2000 + HEGaddn:
                                    residue.resnum = 3000 + HEGaddn
                                elif residue.resnum == 2000 + end:
                                    print('hit 2000')
                                    residue.resnum = 3000 + end
                                elif residue.resnum == 2000 + start:
                                    residue.resnum = 3000 + start
                                else:
                                    residue.resnum = residue.resnum + 1000 + ipillar

                        ## move first linker (bottom right)
                        #transform_matrix = transform.get_alignment_matrix(numpy.array([0, 1, 0]), numpy.array([1, 0, 0]))
                        #transform.transform_structure(link, transform_matrix, atom_index_list)
                        schrodinger.structutils.transform.rotate_structure(link, x_angle=0.25*pi, y_angle=0, z_angle=0, rot_center=None)
                        (a,b,c) = numpy.array([0.45*width, -0.5*width2, 0.15*width]) - atom.xyz
                        transform.translate_structure(link, a, b, c, atom_index_list)
                        st = st.merge(link)

                        ## rename third linker residues
                        for residue in link.residue:
                                if residue.resnum == 3000 + HEGaddn:
                                    residue.resnum = 1000 + HEGaddn
                                elif residue.resnum == 3000 + end:
                                    residue.resnum = 1000 + end
                                elif residue.resnum == 3000 + start:
                                    residue.resnum = 1000 + start
                                else:
                                    residue.resnum = (residue.resnum) + 2000 + ipillar

                        ## move second linker (top left)
                        #transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, 0, -1]))
                        #transform.transform_structure(link, transform_matrix, atom_index_list)
                        schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=0, z_angle=pi, rot_center=None)
                        (a,b,c) = numpy.array([-0.45*width, 0.5*width2, 0.15*width]) - atom.xyz
                        transform.translate_structure(link, a, b, c, atom_index_list)
                        st = st.merge(link)

                        ## rename fourth linker residues
                        for residue in link.residue:
                                if residue.resnum == 1000 + HEGaddn:
                                    residue.resnum = 4000 + HEGaddn
                                elif residue.resnum == 1000 + end:
                                    residue.resnum = 4000 + end
                                elif residue.resnum == 1000 + start:
                                    residue.resnum = 4000 + start
                                else:
                                    residue.resnum = (residue.resnum) + 3000 + ipillar

                        ## move fourth linker residues
                        #schrodinger.structutils.transform.rotate_structure(link, x_angle=0, y_angle=0, z_angle=0, rot_center=None)
                        #transform_matrix = transform.get_alignment_matrix(numpy.array([0, -1, 0]), numpy.array([1, 0, 0]))
                        #transform.transform_structure(link, transform_matrix, atom_index_list)
                        (a,b,c) = numpy.array([-0.52*width, -0.5*width2, 0.15*width]) - atom.xyz
                        transform.translate_structure(link, a, b, c, atom_index_list)
                        st = st.merge(link)

            ##------------------------------------------------------------------------------------------------------------
            ## if specified, create and move additions
            if Add:
                ## create and flip first addition
                add1 = to_add.copy()

                ## rotate by 90 degrees
                atom_index_list = range(1, add1.atom_total+1)
                transform_matrix = transform.get_rotation_matrix_from_eulers(-0.5*pi, -0.5*pi, 0)
                transform.transform_structure(add1, transform_matrix, atom_index_list)

                ## move into position
                centroid = transform.get_centroid(add1, atom_list=[atom.index for atom in add1.atom if atom.element=="P"])
                newpos = numpy.array([-0.25*width2, 2.5*chain_length3, 0.4*width2, 0.0])
                (x, y, z, q) = newpos - centroid
                transform.translate_structure(add1, x, y, z, atom_index_list)

                ## rename chain
                for chain in add1.chain:
                    chain.name = chr(119)

                ## create and rotate second addition
                add2 = to_add.copy()
                atom_index_list = range(1, add2.atom_total+1)
                transform_matrix = transform.get_rotation_matrix_from_eulers(-0.5*pi, -0.5*pi, 0)
                transform.transform_structure(add2, transform_matrix, atom_index_list)

                ## move
                centroid = transform.get_centroid(add2, atom_list=[atom.index for atom in add1.atom if atom.element=="P"])
                newpos = numpy.array([0.25*width2, 2.5*chain_length3, 0.4*width2, 0.0])
                (x, y, z, q) = newpos - centroid
                transform.translate_structure(add2, x, y, z, atom_index_list)

                ## rename chain
                for chain in add2.chain:
                    chain.name = chr(120)

                add1 = add1.merge(add2, copy_props=True)

                ## repeat with second pair
                add_ct3 = add1.copy()

                #add_ct3.write('ct3_test.pdb')

                ## rotate by 90 degrees
                atom_index_list = range(1, add_ct3.atom_total+1)
                transform_matrix = transform.get_rotation_matrix_from_eulers(0, 0, pi)
                transform.transform_structure(add1, transform_matrix, atom_index_list)

                additions = add1.merge(add_ct3, copy_props=True)
                st = st.merge(additions, copy_props=True)



            ##------------------------------------------------------------------------------------------------------------
            ## initialize bonding lists
            deleteme = []
            connectme = []
            linkme = []
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

            ## if addition specified, get residues to delete
            if Add:
                res_top = []
                res_bottom = []
                for chain in st.chain:
                    if chain.name == chr(65+(npillar*2)+ipillar):
                        for res in chain.residue:
                            res_top.append(res)
                    if chain.name == chr(65+npillar+ipillar):
                        for res in chain.residue:
                            res_bottom.append(res)

            ## bond or delete residues
            for residue in st.residue:
                ## if addition specified, delete first/last 3 bp and bond
                if Add:
                    if residue.chain == chr(120) and residue.resnum == 1:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2" or atom.pdbname.strip()=="OP3":
                                deleteme.append(atom.index)
                            if atom.pdbname.strip()=="O5'":
                                addme.append(atom.index)
                    if residue.chain == chr(119) and residue.resnum == 1:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="OP2":
                                deleteme.append(atom.index)
                            if atom.pdbname.strip()=="P":
                                addme.append(atom.index)
                    ## 3' end
                    if residue in res_bottom[-1:-5:-1]:
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    ## 5' end
                    if residue in res_bottom[0:4]:
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    if residue in res_top[-1:-5:-1]:
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    if residue in res_top[0:4]:
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    if residue == res_bottom[-5]:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="O3'":
                                connectme_add.append(atom.index)
                    if residue == res_bottom[4]:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="P":
                                connectme_add.append(atom.index)
                    if residue == res_top[-5]:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="O3'":
                                connectme_add.append(atom.index)
                    if residue == res_top[4]:
                        for atom in residue.atom:
                            if atom.pdbname.strip()=="P":
                                connectme_add.append(atom.index)

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
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="O3*":
                                #         deleteme.append(atom.index)


                if residue.resnum == 1000 + ipillar or residue.resnum == 1000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="O1P" or atom.pdbname.strip()=="O2P" or atom.pdbname.strip()=="O5*":
                                #         deleteme.append(atom.index)


                ## link atoms from second linker to columns
                if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="O3*":
                                #         deleteme.append(atom.index)
                if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="O1P" or atom.pdbname.strip()=="O2P" or atom.pdbname.strip()=="O5*":
                                #         deleteme.append(atom.index)

                ## link atoms from third linker to columns
                if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="O3*":
                                #         deleteme.append(atom.index)
                if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="O1P" or atom.pdbname.strip()=="O2P" or atom.pdbname.strip()=="O5*":
                                #         deleteme.append(atom.index)

                ## link atoms from fourth linker to columns
                if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + end:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerE_name:
                                    linkme.append(atom.index)
                                # if not HEG:
                                #     if atom.pdbname.strip()=="O3*":
                                #         deleteme.append(atom.index)

                if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + start:
                        for atom in residue.atom:
                                if atom.pdbname.strip() == linkerS_name:
                                    linkme.append(atom.index)
                                # # if not HEG:
                                #     if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="O1P" or atom.pdbname.strip()=="O2P" or atom.pdbname.strip()=="O5*":
                                #         deleteme.append(atom.index)


            if Add:
                st.addBonds([(addme[0], connectme_add[1], 1),
                             (addme[1], connectme_add[3], 1),
                             (addme[2], connectme_add[0], 1),
                             (addme[3], connectme_add[2], 1)])



            print('linkme', linkme)
            for atom in st.atom:
                   if atom.index in linkme:
                       print(atom.pdbname + ', ' + str(atom.index))


            print('connectme', connectme)
            for atom in st.atom:
                   if atom.index in connectme:
                       print(atom.pdbname + ', ' + str(atom.index))


            ## add bonds dont use linkme 1/6
            ## should be c8 to o5
            if linkerfile == '4T.pdb':
                st.addBonds([(linkme[0], connectme[0], 1), #
                             (linkme[3], connectme[1], 1), #
                             (linkme[2], connectme[2], 1), # fix
                             (linkme[5], connectme[3], 1), #
                             (linkme[4], connectme[4], 1), # fix
                             (linkme[7], connectme[5], 1)
                             ]) # c8 o3
            else:
                st.addBonds([(linkme[1], connectme[0], 1), # c12 o5
                             (linkme[3], connectme[2], 1), # c8 o3
                             (linkme[2], connectme[1], 1), # c12 o5
                             (linkme[5], connectme[4], 1), # c8 o3
                             (linkme[4], connectme[3], 1), # c12 o5
                             (linkme[6], connectme[5], 1)
                             ]) # c8 o3




        #    for j in range(len(linkme)):
        #	struct.addBond(connectme[j],linkme[j],1)
            st.write('temp.pdb')
            ## delete extraneous atoms
            st.deleteAtoms(deleteme)
            connectme = []
            deleteme = []
            linkme = []


            #for residue in st.residue:
             #   if residue.chain == chain_name:
              #          residue.resnum = iterator
               #         iterator = iterator + 1

        #    ct_2.write("test.pdb")
        #    break

            ## add wall to entire structure
            wall_list.append(st)
        #   st.write("st"+chain_name+".pdb")
            del st
        #    sys.exit(1)

        ##------------------------------------------------------------------------------------------------------------
        print('start')
        #for residue in wall_list[0].residue:
         #   print(residue.resnum)
        dist = width #change to increase corner distances

        ## create cage as object
        st = dev.cage_builder(wall_list, dist, q_chain_length, linkerS_name, linkerE_name, linkerS_resnum, linkerE_resnum)

        ## rename chains
        #i = 65
        #for chain in st.chain:
        #    chain.name = chr(i)
        #    i += 1

        outfolder = str(npillar) + '_' + dna_struct.name + '/'
        print(savepath + outfolder + dna_struct.name + ".pdb")
        if not os.path.isdir(savepath + outfolder):
            os.makedirs(savepath + outfolder)

        for chain in st.chain:
            i = 0
            for residue in chain.residue:
                residue.resnum = i
                i += 1

        ## minimize structure and export
        if dna_struct.Min:
            min = minimize.Minimizer(struct=st)

            min.updateCoordinates(st)
            min.minimize()
            st_min = min.getStructure()
        else:
            st_min=st.copy()

        st_min.write(savepath + outfolder + str(npillar) + "_full_" + dna_struct.name + ".pdb")

        ## delete extra atoms and export
        #deleteme=[]
        #k = 0
        #for atom in st_min.atom:
        #    if atom.chain > chr(65+len(wall_list)-1):
        #        deleteme.append(atom.index)

        st_min.deleteAtoms(deleteme)

        st_min.write(savepath + outfolder + str(npillar) + '_' + dna_struct.name + ".pdb")

        ## add Mg to balance charge
        if dna_struct.Mg:
            k=st_min.formal_charge
            random.seed(17)
            mg = next(structure.StructureReader("Mg.pdb"))
            #print(i)
            j = 0
            angle = 2*numpy.pi/len(wall_list)
            radius = dist/(2*numpy.tan(angle/2))

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
            maggedst.write(savepath + outfolder + "mg_" + dna_struct.name + ".pdb")
