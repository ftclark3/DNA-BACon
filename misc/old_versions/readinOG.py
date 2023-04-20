import random
import sys, numpy
import build_dna2b
import dev
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from sklearn.externals import joblib


    
##-------------------

## read the input file
inputfile = sys.argv[1]

DNA_seq = []
wall_list = []


#linkers = input("Linker start atom number: ")
#linkers_name = raw_input("Linker start atom name: ")
#linkere = input("Linker end atom number: ")
#linkere_name = raw_input("Linker end atom name: ")
#linker_center = raw_input("Linker center atom name: ")

linkerfile = sys.argv[2]

#4T
#linkers = 5
#linkers_name = "O5*"
#linkere = 104
#linkere_name = "C3*"
#linker_center = "P"

#HEG
linkers = 11
linkers_name = "C8"
linkere = 17
linkere_name = "C12"
linker_center = "O1"

init_link = structure.StructureReader(linkerfile).next()
addendum = init_link.measure(linkers, linkere)
iterator = 5000

try:
    dna_structures = dev.file_reader(sys.argv[1])
except:
    dna_structures = dev.file_reader(raw_input("It appears you did not input a valid file, please type it now: "))

for dna_struct in dna_structures:
    DNA_seq = dna_struct.DNA_seq

    npillar = len(DNA_seq)
    chain_length = len(DNA_seq[0][1])
    chain_length2 = len(DNA_seq[0][2])
    chain_length3 = len(DNA_seq[0][3])
    q_chain_length = []
    q_chain_length.append(len(DNA_seq[0][1]))
    q_chain_length.append(len(DNA_seq[0][2]))
    q_chain_length.append(len(DNA_seq[0][3]))
    print chain_length2
    width = 3.7 * chain_length + addendum
    width2 = 3.7 * chain_length2 + addendum

    ## create the cage
    for ipillar in range(npillar-1, -1, -1):
        print ipillar
        chain_name = chr(65+ipillar)
        print chain_name, ipillar

        ##-----
        ct_1 = build_dna2b.process_sequence(DNA_seq[ipillar][1])

        #ct_1.write("ct_1_0.pdb")
      
        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([1, 0, 0]))    
        atom_index_list = range(1, ct_1.atom_total+1)
        transform.transform_structure(ct_1, transform_matrix, atom_index_list)    

        centroid = transform.get_centroid(ct_1, atom_list=[atom.index for atom in ct_1.atom if atom.element=="P"])
        #print centroid

        newpos = numpy.array([0.0, 0.5*width2, 0.0, 0.0])
        #print newpos

        (x, y, z, q) = newpos - centroid
        #print x
        transform.translate_structure(ct_1, x, y, z, atom_index_list)

        for chain in ct_1.chain:
            if chain.name == "A":
                chain.name = chain_name
            else:
                chain.name = chr(65+ipillar+npillar)
        
        for res in ct_1.residue:
            res.resnum = len(DNA_seq[ipillar][0])+res.resnum-1
        
            #ct_1.write("ct_1.pdb")    
               
        ##------
        ct_3 = build_dna2b.process_sequence(DNA_seq[ipillar][3])

        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([-1, 0, 0]))    
        atom_index_list = range(1, ct_3.atom_total+1)
        transform.transform_structure(ct_3, transform_matrix, atom_index_list)    

        centroid = transform.get_centroid(ct_3, atom_list=[atom.index for atom in ct_3.atom if atom.element=="P"])
        #print centroid

        newpos = numpy.array([0.0, -0.5*width2, 0.0, 0.0])
        #print newpos

        (x, y, z, q) = newpos - centroid
        #print x
        transform.translate_structure(ct_3, x, y, z, atom_index_list)
      
        for chain in ct_3.chain:
            if chain.name == "A":
                chain.name = chain_name
            else:
                chain.name = chr(65+ipillar+2*npillar)
        
        for res in ct_3.residue:
            res.resnum = len(DNA_seq[ipillar][0]+DNA_seq[ipillar][1]+DNA_seq[ipillar][2])+res.resnum-1

            #ct_3.write("ct_3.pdb")

        ##------

        ct_2 = build_dna2b.process_sequence(DNA_seq[ipillar][2])


        transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, -1, 0]))    
        atom_index_list = range(1, ct_2.atom_total+1)
        transform.transform_structure(ct_2, transform_matrix, atom_index_list)    

        centroid = transform.get_centroid(ct_2, atom_list=[atom.index for atom in ct_2.atom if atom.element=="P"])
        #print centroid

        newpos = numpy.array([ 0.5*width, 0.0, 0.0, 0.0])
        #print newpos

        (x, y, z, q) = newpos - centroid
        #print x
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
                        print atom1[0].resnum, atom1[0].pdbname, atom2[0].resnum, atom2[0].pdbname
                        ct_2.deleteBond(atom1[0], atom2[0])
                        add_hydrogens(ct_2)
                        atom3 = [atom for atom in atom1[0].bonded_atoms if atom.element == "H"]
                        atom3[0].element = "O"
                        atom3[0].pdbname = " OP3"
                        atom3[0].formal_charge=-1
                                           
    #                    ct_2.write("ct_2.pdb")

        
        st = ct_1.merge(ct_2,copy_props=True)
        st = st.merge(ct_3,copy_props=True)
        ##-------------------------------------------------------------------


        link =init_link.copy()
        print('link created')
        atom_index_list = range(1, link.atom_total+1)

        for residue in link.residue:
            print('residue testing')
            if residue.resnum == linkers:
                residue.resnum = 2000 + ipillar + linkers
                print("initial firing")
            elif residue.resnum == linkere:
                residue.resnum = 2000 + ipillar + linkere
                print('initial firing II')
            else:
                residue.resnum=2000 + ipillar

        for atom in link.atom:
            if atom.pdbname.strip()==linker_center:
                    (a,b,c) = numpy.array([0.75*width, 0.5*width2, -10]) - atom.xyz
                    transform.translate_structure(link, a, b, c, atom_index_list)
                    st = st.merge(link)
                    
                    for residue in link.residue:
                        if residue.resnum == 2000 + ipillar + linkers:
                            residue.resnum = 3000 + ipillar + linkers
                        elif residue == 2000 + ipillar + linkere:
                            residue.resnum = 3000 + ipillar + linkere
                        else:
                            residue.resnum = 3000 + ipillar
                    
                    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 1, 0]), numpy.array([1, 0, 0]))    
                    transform.transform_structure(link, transform_matrix, atom_index_list)
                    (a,b,c) = numpy.array([0.75*width, -0.5*width2, -5.0]) - atom.xyz
                    transform.translate_structure(link, a, b, c, atom_index_list)
                    st = st.merge(link)

                    for residue in link.residue:
                        for residue in link.residue:
                            if residue.resnum == 2000 + ipillar + linkers:
                                residue.resnum = 1000 + ipillar + linkers
                            elif residue == 2000 + ipillar + linkere:
                                residue.resnum = 1000 + ipillar + linkere
                            else:
                                residue.resnum = 1000 + ipillar

                    
                    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, 0, -1]))    
                    transform.transform_structure(link, transform_matrix, atom_index_list)
                    (a,b,c) = numpy.array([-0.75*width, 0.5*width2, -5.0]) - atom.xyz
                    transform.translate_structure(link, a, b, c, atom_index_list)
                    st = st.merge(link)


                    for residue in link.residue:
                        for residue in link.residue:
                            if residue.resnum == 2000 + ipillar + linkers:
                                residue.resnum = 4000 + ipillar + linkers
                            elif residue == 2000 + ipillar + linkere:
                                residue.resnum = 4000 + ipillar + linkere
                            else:
                                residue.resnum = 4000 + ipillar
                    
                    transform_matrix = transform.get_alignment_matrix(numpy.array([0, -1, 0]), numpy.array([1, 0, 0]))    
                    transform.transform_structure(link, transform_matrix, atom_index_list)    
                    (a,b,c) = numpy.array([-0.75*width, -0.5*width2, -5.0]) - atom.xyz
                    transform.translate_structure(link, a, b, c, atom_index_list)
                    st = st.merge(link)



        deleteme = []
        connectme = []
        linkme = []

    #addendum = link.measure(11,17)
        i = 1


        for atom in st.atom:
            if atom.pdbres.strip() == "POT" or atom.pdbres.strip() == "HXL":
                    if atom.chain < chr(65+npillar):
                            deleteme.append(atom.index)
        st.deleteAtoms(deleteme)
        deleteme = []
        for residue in st.residue:
            if residue.chain == chain_name:
    #		if i == 0:
    #			for atom in residue.atom:
    #				deleteme.append(atom.index)
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
            
            if residue.resnum == 1000 + ipillar + linkere:
                print('A1')
                linkme.append(atom.index)
            if residue.resnum == 1000 + ipillar + linkers:
                print('A2')
                linkme.append(atom.index)


            if residue.resnum == 2000 + ipillar + linkere:
                print("A")
                linkme.append(atom.index)
            if residue.resnum == 2000 + ipillar + linkers:
                print("B")
                linkme.append(atom.index)


            if residue.resnum == 3000 + ipillar + linkere:
                print("C")
                linkme.append(atom.index)
            if residue.resnum == 3000 + ipillar + linkers:
                print("D")
                linkme.append(atom.index)


            if residue.resnum == 4000 + ipillar + linkere:
                linkme.append(atom.index)
                print("E")
            if residue.resnum == 4000 + ipillar + linkers:
                linkme.append(atom.index)
                print("F")



        #print(type(linkme[0]))
        #print(type(connectme[0]))
        #connectme[1],connectme[2], connectme[3]. connectme[4], connectme[5]))
            
        st.addBonds([
                    (linkme[0], connectme[0], 1),
                    (linkme[3], connectme[1], 1),
                    (linkme[2], connectme[2], 1),
                    (linkme[5], connectme[3], 1),
                    (linkme[4], connectme[4], 1),
                    (linkme[7], connectme[5], 1)])
        
    

    #    for j in range(len(linkme)):
    #	struct.addBond(connectme[j],linkme[j],1)
        st.deleteAtoms(deleteme)
        connectme = []
        deleteme = []
        linkme = []
                    
        for residue in st.residue:
            if residue.chain == chain_name:
                    residue.resnum = iterator	
                    iterator = iterator + 1

    #    ct_2.write("test.pdb")
    #    break

        ##-------------------------------------------------------------------

        wall_list.append(st)
    #    st.write("test.pdb")
    #    st.write("st"+chain_name+".pdb")
        del st
        #sys.exit(1)


    #------------------------------------------------------------------------------------------------------------
    dist = width #change to increase corner distances

    print(type(wall_list))
    print(type(dist))
    print(type(q_chain_length))

    st = dev.cage_builder(wall_list,dist,q_chain_length, linkers_name, linkere_name)

    if dna_struct.Min:
        min = minimize.Minimizer(struct=st)
           
        min.updateCoordinates(st)
        min.minimize()
        st_min = min.getStructure()
    else:
        st_min=st.copy()

    st_min.write("full" + dna_struct.name + ".pdb")

    deleteme=[]
    k = 0
    for atom in st_min.atom:
        if atom.chain > chr(65+len(wall_list)-1):
            deleteme.append(atom.index)

    st_min.deleteAtoms(deleteme)

    st_min.write(dna_struct.name + ".pdb")

    if dna_struct.Mg:
        k=st_min.formal_charge
        random.seed(17)
        mg = structure.StructureReader("Mg.pdb").next()
        print i
        j = 0
        angle = 2*numpy.pi/len(wall_list)
        radius = dist/(2*numpy.tan(angle/2))

        transform.translate_to_origin(st_min)
        for j in range(-st_min.formal_charge/2 + 1):
            i = i + 1
            mg.atom[1].x = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].y = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].z = random.uniform(-(radius)/2.0,(radius)/2.0)
            mg.atom[1].resnum = i
            if j == 0:
                maggedst=st_min.merge(mg)
                print "hi"
            else:
                maggedst=maggedst.merge(mg)
        maggedst.write("mg" + dna_struct.name + ".pdb")

