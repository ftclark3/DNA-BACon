import random
import sys, numpy
import build_dna2b
import schrodinger
import os
from math import tan, pi
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize
from sklearn.externals import joblib

# tools to build wall using input object and pillar number
##------------------------------------------------------------------------------------------------------------
# class to hold column information
class column:
    def __init__(self, dnaSeqs, ipillar, npillar, columnNum):
        self._dnaSeqs = dnaSeqs
        self._ipillar = ipillar
        self._npillar = npillar
        self._chainName = chr(65+ipillar)
        self._columnNum = columnNum

        # print("DNA: " + str(self._dnaSeqs))
        # print("Ipillar: " + str(self._ipillar))
        # print("Column: " + str(self._columnNum))

    # method to make column and move it
    def makeColumn(self, transform_matrix, newpos, rotation=None):
        ## build column
        # print("columnNum: " + str(self._columnNum))
        # print("seq: " + str(self._dnaSeqs[self._columnNum]))
        ct = build_dna2b.process_sequence(self._dnaSeqs[self._columnNum])
        ## calculate atom index
        atom_index_list = range(1, ct.atom_total+1)
        ## rotate into posn
        transform.transform_structure(ct, transform_matrix, atom_index_list)
        ## rotate if necessary
        if rotation is not None:
            schrodinger.structutils.transform.rotate_structure(ct, x_angle=rotation[0], y_angle=rotation[1], z_angle=rotation[2], rot_center=rotation[3])
        ## translate into position
        centroid = transform.get_centroid(ct, atom_list=[atom.index for atom in ct.atom if atom.element=="P"])
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(ct, x, y, z, atom_index_list)
        ## rename chains
        self._chainsResnums(ct)
        return ct

    # method to rename chains and renumber resnums according to columnNum
    def _chainsResnums(self, ct):
        # if first column
        if self._columnNum == 1:
            for chain in ct.chain:
                if chain.name == "A":
                    chain.name = self._chainName
                else:
                    chain.name = chr(65+(self._npillar*2)+self._ipillar)
            for res in ct.residue:
                res.resnum = len(self._dnaSeqs[0])+res.resnum-1
        # if third column
        elif self._columnNum == 3:
            for chain in ct.chain:
                if chain.name == "A":
                    chain.name = self._chainName
                else:
                    chain.name = chr(65+self._npillar+self._ipillar)
            for res in ct.residue:
                res.resnum = len(self._dnaSeqs[0]+self._dnaSeqs[1]+self._dnaSeqs[2])+res.resnum-1
        # if second column
        elif self._columnNum == 2:
            for chain in ct.chain:
                if chain.name == "A":
                    chain.name = self._chainName
                    for res in chain.residue:
                        res.resnum = len(self._dnaSeqs[0]+self._dnaSeqs[1])+res.resnum-1
                else:
                    if self._ipillar > 0:
                        chain.name = chr(65+self._ipillar-1)
                    else:
                        chain.name = chr(65+self._npillar-1)
                    n1 = len(self._dnaSeqs[0])
                    n2 = len(self._dnaSeqs[0]+self._dnaSeqs[1]+self._dnaSeqs[2]+self._dnaSeqs[3])
                    for res in chain.residue:
                        if res.resnum < n1+2:
                            res.resnum = n2+res.resnum-1
                        else:
                            res.resnum = res.resnum-n1-1
                        if res.resnum == 1:
                            atom1 = [atom for atom in res.atom if atom.pdbname.strip() == "P"]
                            atom2 = [atom for atom in atom1[0].bonded_atoms if atom.resnum >1]
                            ct.deleteBond(atom1[0], atom2[0])
                            add_hydrogens(ct)
                            atom3 = [atom for atom in atom1[0].bonded_atoms if atom.element == "H"]
                            atom3[0].element = "O"
                            atom3[0].pdbname = " OP3"
                            atom3[0].formal_charge=-1

##------------------------------------------------------------------------------------------------------------
# class to hold linker information
class linker:
    def __init__(self, end, start, center, linkerfile, linkerNum):
        self._endResnum = end
        self._startResnum = start
        self._centerAtom = center
        self._linkerNum = linkerNum
        self._link = next(structure.StructureReader(linkerfile))

    # method to make linker and move it
    def makeLinker(self, newpos, rotation, ipillar):
        # get atom index list
        atom_index_list = range(1, self._link.atom_total+1)
        # iterate through atoms to find the center
        for atom in self._link.atom:
            if atom.pdbname.strip() == self._centerAtom:
                schrodinger.structutils.transform.rotate_structure(self._link, x_angle=rotation[0], y_angle=rotation[1], z_angle=rotation[2], rot_center=None)
                (a, b, c) = newpos - atom.xyz
                transform.translate_structure(self._link, a, b, c, atom_index_list)
        # update resnum and chain.name
        self._chainsResnums(ipillar)
        # for res in self._link.residue:
        #     print(res.resnum)
        return self._link

    # method to rename chains and renumber resnums according to linkerNum
    def _chainsResnums(self, ipillar):
        # initialize resnum positional additions
        end = self._endResnum * 10 + ipillar
        start = self._startResnum * 10 + ipillar
        hegAddn = self._startResnum * 10 + ipillar
        # iterate through chains and name them
        for chain in self._link.chain:
            chain.name = 'a'
        # iterate through residues
        for residue in self._link.residue:
            # top right
            if self._linkerNum == 2:
                ## if using HEG linker
                if self._endResnum == self._startResnum:
                    residue.resnum = 2000 + hegAddn
                ## if on last linker resnum
                elif residue.resnum == self._endResnum:
                    residue.resnum = 2000 + end
                ## if on first linker resnum
                elif residue.resnum == self._startResnum:
                    residue.resnum = 2000 + start
                ## if on middle linker resnum
                else:
                    residue.resnum = (residue.resnum) + 2000 + (ipillar*100)
            # bottom right linker
            elif self._linkerNum == 3:
                if self._endResnum == self._startResnum:
                    residue.resnum = 3000 + hegAddn
                ## if on last linker resnum
                elif residue.resnum == self._endResnum:
                    residue.resnum = 3000 + end
                ## if on first linker resnum
                elif residue.resnum == self._startResnum:
                    residue.resnum = 3000 + start
                ## if on middle linker resnum
                else:
                    residue.resnum = residue.resnum + 3000 + ipillar
            # top left
            elif self._linkerNum == 1:
                if self._endResnum == self._startResnum:
                    residue.resnum = 1000 + hegAddn
                ## if on last linker resnum
                elif residue.resnum == self._endResnum:
                    residue.resnum = 1000 + end
                ## if on first linker resnum
                elif residue.resnum == self._startResnum:
                    residue.resnum = 1000 + start
                ## if on middle linker resnum
                else:
                    residue.resnum = residue.resnum + 1000 + ipillar
            # bottom left
            elif self._linkerNum == 4:
                if self._endResnum == self._startResnum:
                    residue.resnum = 4000 + hegAddn
                ## if on last linker resnum
                elif residue.resnum == self._endResnum:
                    residue.resnum = 4000 + end
                ## if on first linker resnum
                elif residue.resnum == self._startResnum:
                    residue.resnum = 4000 + start
                ## if on middle linker resnum
                else:
                    residue.resnum = residue.resnum + 4000 + ipillar

##------------------------------------------------------------------------------------------------------------
# function to build columns of the wall
def buildColumns(seqs, ipillar, npillar, width, width2):
    # first column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([1, 0, 0]))
    newpos = numpy.array([0.0, 0.5*width2, 0.0, 0.0])
    c1 = column(seqs, ipillar, npillar, 1)
    c1 = c1.makeColumn(transform_matrix, newpos)

    # third column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([-1, 0, 0]))
    newpos = numpy.array([0.0, -0.5*width2, 0.0, 0.0])
    c3 = column(seqs, ipillar, npillar, 3)
    c3 = c3.makeColumn(transform_matrix, newpos)

    # second (middle) column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, -1, 0]))
    newpos = numpy.array([0.5*width, 0.0, 0.0, 0.0])
    rotation = [0, 0.5*pi,0, None]
    c2 = column(seqs, ipillar, npillar, 2)
    c2 = c2.makeColumn(transform_matrix, newpos, rotation)

    # merge columns
    columns = c1.merge(c2, copy_props=True)
    columns = columns.merge(c3, copy_props=True)
    return columns

##------------------------------------------------------------------------------------------------------------
# function to build linkers of the wall
def buildLinkers(endResnum, startResnum, centerAtom, linkerFile, ipillar, width, width2):
    # MAY NOT BE WORKING
    # first linker (top left)
    l1 = linker(endResnum, startResnum, centerAtom, linkerFile, 1)
    if linkerFile == "linkers/4T.pdb":
        l1 = l1.makeLinker(numpy.array([-0.5*width, 0.5*width2, 0.25*width]), [0, 0, pi/2], ipillar)
    else:
        l1 = l1.makeLinker(numpy.array([-0.45*width, 0.5*width2, 0.15*width]), [0, 0, pi], ipillar)

    # second linker (top right)
    l2 = linker(endResnum, startResnum, centerAtom, linkerFile, 2)
    l2 = l2.makeLinker(numpy.array([0.45*width, 0.52*width2, 0.2*width]), [0, -0.5*pi, 0], ipillar)

    # third linker (bottom right)
    l3 = linker(endResnum, startResnum, centerAtom, linkerFile, 3)
    if linkerFile =="linkers/4T.pdb":
        l3 = l3.makeLinker(numpy.array([0.45*width, -0.5*width2, 0.15*width]), [-0.25*pi, 0, 0.25*pi], ipillar)
    else:
        l3 = l3.makeLinker(numpy.array([0.45*width, -0.5*width2, 0.15*width]), [0.25*pi, 0, 0], ipillar)

    # fourth linker (bottom left)
    l4 = linker(endResnum, startResnum, centerAtom, linkerFile, 4)
    if linkerFile == "linkers/4T.pdb":
        l4 = l4.makeLinker(numpy.array([-0.4*width, -0.5*width2, 0.25*width]), [pi, 0, pi/2], ipillar)
    else:
        l4 = l4.makeLinker(numpy.array([-0.4*width, -0.5*width2, 0.15*width]), [0, 0, 0], ipillar)

    # merge linkers
    linkers = l1.merge(l2, copy_props=True)
    linkers = linkers.merge(l3, copy_props=True)
    linkers = linkers.merge(l4, copy_props=True)
    # for res in linkers.residue:
    #     print(res.resnum)
    return linkers

##------------------------------------------------------------------------------------------------------------
# function to take build wall out of linkers and columns
def buildWall(struct, ipillar, npillar, width, width2):
    # extract info from wall input
    seqs = struct.dnaSeqs[ipillar]
    endResnum = struct.linkerEResnum
    startResnum = struct.linkerSResnum
    centerAtom = struct.linkerCenter
    linkerFile = struct.linker
    linkerEName = struct.linkerEName
    linkerSName = struct.linkerSName
    linkerEResnum = struct.linkerEResnum
    linkerSResnum = struct.linkerSResnum

    # measurement variables
    chainName = chr(65+ipillar)
    # print(seqs)
    chainLength = len(seqs[1])
    chainLength2 = len(seqs[2])
    chainLength3 = len(seqs[3])
    # print("ChainLength1: " + str(chainLength))
    # print("ChainLength2: " + str(chainLength2))
    # print("ChainLength3: " + str(chainLength3))
    end = linkerEResnum * 10 + ipillar
    start = linkerSResnum * 10 + ipillar

    # generate columns and linkers
    columns = buildColumns(seqs, ipillar, npillar, width, width2)
    linkers = buildLinkers(endResnum, startResnum, centerAtom, linkerFile, ipillar, width, width2)
    wall= columns.merge(linkers, copy_props=True)

    # bond columns to linkers
    ## initialize bonding lists
    deleteMe = []
    connectMe = []
    linkMeS = []
    linkMeE = []
    i = 1

    # delete ends(?)
    for atom in wall.atom:
        if atom.pdbres.strip() == "POT" or atom.pdbres.strip() == "HXL":
                if atom.chain < chr(65+npillar+ipillar):
                        deleteMe.append(atom.index)
    wall.deleteAtoms(deleteMe)
    deleteMe = []

    ## bond or delete residues
    for residue in wall.residue:
        ## if not a linker residue, delete extraneous atoms
        if residue.chain == chainName and residue.chain != 'a':
            if i == 1:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                    deleteMe.append(atom.index)
                            if atom.pdbname.strip()=="O5'":
                                    connectMe.append(atom.index)
            if i == chainLength:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="O3'":
                                    connectMe.append(atom.index)
            if i == 1 + chainLength:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                    deleteMe.append(atom.index)
                            if atom.pdbname.strip()=="O5'":
                                    connectMe.append(atom.index)
            if i == chainLength + chainLength2:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="O3'":
                                    connectMe.append(atom.index)
            if i == 1 + chainLength + chainLength2:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="P" or atom.pdbname.strip()=="OP1" or atom.pdbname.strip()=="OP2":
                                    deleteMe.append(atom.index)
                            if atom.pdbname.strip()=="O5'":
                                    connectMe.append(atom.index)
            if i == chainLength + chainLength2 + chainLength3:
                    for atom in residue.atom:
                            if atom.pdbname.strip()=="O3'":
                                    connectMe.append(atom.index)
            i = i + 1

        ## link atoms from first linker to columns
        if residue.resnum == 1000 + ipillar or residue.resnum == 1000 + end:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerEName:
                            linkMeE.append(atom.index)
        if residue.resnum == 1000 + ipillar or residue.resnum == 1000 + start:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerSName:
                            linkMeS.append(atom.index)
        ## link atoms from second linker to columns
        if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + end:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerEName:
                            linkMeE.append(atom.index)
        if residue.resnum == 2000 + ipillar or residue.resnum == 2000 + start:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerSName:
                            linkMeS.append(atom.index)
        ## link atoms from third linker to columns
        if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + end:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerEName:
                            linkMeE.append(atom.index)
        if residue.resnum == 3000 + ipillar or residue.resnum == 3000 + start:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerSName:
                            linkMeS.append(atom.index)
        ## link atoms from fourth linker to columns
        if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + end:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerEName:
                            linkMeE.append(atom.index)
        if residue.resnum == 4000 + ipillar or residue.resnum == 4000 + start:
                for atom in residue.atom:
                        if atom.pdbname.strip() == linkerSName:
                            linkMeS.append(atom.index)

    ## add bonds dont use linkMe 1/6
    ## should be c8 to o5s
    wall.addBonds([(linkMeS[2], connectMe[4], 1), # c12 o5
                 (linkMeS[1], connectMe[2], 1), # c8 o3
                 (linkMeE[3], connectMe[5], 1), # c12 o5
                 (linkMeS[0], connectMe[0], 1), # c8 o3
                 (linkMeE[1], connectMe[1], 1), # c12 o5
                 (linkMeE[2], connectMe[3], 1)
                 ]) # c8 o3
    wall.deleteAtoms(deleteMe)
    # return wall structure
    wall.write("models/cages/wall.pdb")
    return wall
