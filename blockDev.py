from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from seq_writer import gen_DNA
import build_dev as dev
import wallDev as wd
import numpy as np
import build_dna2b
import toml
import sys
import os

''' Classes used to build junction blocks '''

class block:
    ##------------------------------------------------------------------------------------------------------------
    # constructor
    def __init__(self, iblock, nblock, cname, width, ires, sequence, segment_length, wombat):
        self.iblock = iblock
        self.nblock = nblock
        self.cname = cname
        self.width = width
        self.ires = ires
        self.sequence = sequence
        self.segment_length = segment_length
        self.wombat = wombat

        # call private methods to build strands
        self.top = self._buildTop()
        self.mid = self._buildMid()
        self.bottom = self._buildBottom()
        self.st = self.top.merge(self.mid, copy_props=True)
        self.st = self.st.merge(self.bottom, copy_props=True)

    ##------------------------------------------------------------------------------------------------------------
    # method to build top DNA strand
    def _buildTop(self):
        ## create and move top pillar
        top = build_dna2b.process_sequence(self.sequence[0][:(self.segment_length*2 + 1)])

        ## rotate pillar
        transform_matrix = transform.get_alignment_matrix(np.array([0, 0, 1]), np.array([1, 0, 0]))
        atom_index_list = range(1, top.atom_total+1)
        transform.transform_structure(top, transform_matrix, atom_index_list)

        ## get centroid of top pillar
        centroid = transform.get_centroid(top, atom_list=[atom.index for atom in top.atom if atom.element=="P"])
        newpos = np.array([0.0, -self.width, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(top, x, y, z, atom_index_list)

        ## renumber residues
        i = 0
        for res in top.residue:
            res.resnum = 1000 + self.ires + i
            i += 1

        ## rename chain
        ## capital for outer chain / lower for inner
        for chain in top.chain:
            if chain.name == 'A':
                if self.iblock == (self.nblock-1):
                    chain.name = self.cname[0][self.nblock-1]
                else:
                    chain.name = self.cname[0][self.iblock]

            else:
                if self.iblock == (self.nblock-1):
                    chain.name = self.cname[1][self.nblock-1]
                else:
                    chain.name = self.cname[1][self.iblock]

        return top

    ##------------------------------------------------------------------------------------------------------------
    # method to build middle DNA strand
    def _buildMid(self):
        mid = build_dna2b.process_sequence(gen_DNA(self.segment_length))

        ## delete second strand
        deleteme = []
        for chain in mid.chain:
            if chain.name == 'A':
                for atom in chain.atom:
                    deleteme.append(atom.index)
        mid.deleteAtoms(deleteme)

        ## rotate
        transform_matrix = transform.get_alignment_matrix(np.array([0, 0, 1]), np.array([0, -1, 0]))
        atom_index_list = range(1, mid.atom_total+1)
        transform.transform_structure(mid, transform_matrix, atom_index_list)

        ## get centroid
        centroid = transform.get_centroid(mid, atom_list=[atom.index for atom in mid.atom if atom.element=="P"])
        atom_index_list = range(1, mid.atom_total+1)
        newpos = np.array([0.0, 0.0, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(mid, x, y, z, atom_index_list)

        ## rename residues
        self.lastres = 0
        for res in mid.residue:
            self.lastres = res.resnum

        ## rename chain
        for chain in mid.chain:
            chain.name = chr(109+self.wombat)
        self.wombat+=1

        return mid

    ##------------------------------------------------------------------------------------------------------------
    # method to build bottom DNA strand
    def _buildBottom(self):
        bottom = build_dna2b.process_sequence(self.sequence[0][(self.segment_length*2 + 1):])
        transform_matrix = transform.get_alignment_matrix(np.array([0, 0, 1]), np.array([-1, 0, 0]))
        atom_index_list = range(1, bottom.atom_total+1)
        transform.transform_structure(bottom, transform_matrix, atom_index_list)

        ## get centroid
        centroid = transform.get_centroid(bottom, atom_list=[atom.index for atom in bottom.atom if atom.element=="P"])
        newpos = np.array([0.0, self.width, 0.0, 0.0])

        ## move to refined position
        (x, y, z, q) = newpos - centroid
        transform.translate_structure(bottom, x, y, z, atom_index_list)

        ## rename residues

        i = 0
        for res in bottom.residue:
            res.resnum = 2000 + self.ires + i
            i += 1
        ## rename chains
        for chain in bottom.chain:
            if chain.name == 'A':
                if self.iblock == 0:
                    chain.name = self.cname[0][self.nblock-1]
                else:
                    chain.name = self.cname[0][self.iblock-1]

            else:
                if self.iblock == 0:
                    chain.name = self.cname[1][self.nblock-1]
                else:
                    chain.name = self.cname[1][self.iblock-1]

        return bottom

    ##------------------------------------------------------------------------------------------------------------
    # method to bond all blocks together
    def bond(self):
        # delete hydrogens
        delete_hydrogens(self.st)

        # initialize bonding liself.sts
        deleteme = []
        deletebond_top = []
        deletebond_bottom = []
        convertdouble = []
        connectme_mid = []
        connectme_out = []
        check1 = []
        check2 = []

        ## make resnums variables(?)
        for chain in self.st.chain:
            ## middle
            if chain.name == chr(109+self.wombat-1):
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
                    if atom.resnum == self.lastres:
                        if atom.pdbname.strip() == 'O5T' or atom.pdbname.strip() == 'O1P':
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'O2P':
                            deleteme.append(atom.index)
                    if atom.resnum == self.lastres - 1:
                        if atom.pdbname.strip() == "O3'":
                            connectme_mid.append(atom.index)

            ## top
            if chain.name == self.cname[0][self.iblock] or chain.name == self.cname[0][self.iblock-1]:
                for atom in chain.atom:
                    if atom.resnum == 1000 + self.ires + self.segment_length:
                        if atom.pdbname.strip() == 'P':
                            connectme_out.append(atom.index)
                            deletebond_top.append(atom.index)
                            check1.append(atom.index)
                    if atom.resnum == 1000 + self.ires + self.segment_length - 1:
                        if atom.pdbname.strip() == "O3'":
                            deletebond_top.append(atom.index)
                            check2.append(atom.index)
                    if atom.resnum == 1000 + self.ires + 1:
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'OP2' or atom.pdbname.strip() == 'OP1' or atom.pdbname.strip() == 'O3T':
                            deleteme.append(atom.index)
                    if atom.resnum == 1000 + self.ires:
                        if atom.pdbname.strip() == 'O3T':
                            deleteme.append(atom.index)


            ## bottom
            if chain.name == self.cname[0][self.iblock] or chain.name == self.cname[0][self.iblock-1]:
                for atom in chain.atom:
                    if atom.resnum == 2000 + self.ires + self.segment_length:
                        if atom.pdbname.strip() == 'P':
                            connectme_out.append(atom.index)
                            deletebond_bottom.append(atom.index)
                    if atom.resnum == 2000 + self.ires + self.segment_length - 1:
                        if atom.pdbname.strip() == "O3'":
                            deletebond_bottom.append(atom.index)
                    if atom.resnum == 2000 + self.ires + len(self.sequence[0])//2+1:
                        if atom.pdbname.strip() == 'O5T':
                            # print(atom.index)
                            deleteme.append(atom.index)
                        if atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'O2P':
                            # print('hit' + atom.pdbname.strip())
                            convertdouble.append(atom.index)
                    # if atom.resnum == 2000 + self.ires + len(DNA_seq[0][0]//2+1:
                    #     if atom.pdbname.strip == "O5'" or atom.pdbname.strip() == 'P' or atom.pdbname.strip() == 'OP2' or atom.pdbname.strip() == 'OP1' or atom.pdbname.strip() == 'O3T':


        ## delete double bond and create single bond between P and O2P

        # print('checking bonds')
        # print(check1, check2)
        # print(self.st.atom[check1[0]].pdbname)
        # print(self.st.atom[check2[0]].pdbname)

        # for i in range(len(check1)):
        #     print(self.st.areBound(check1[i], check2[i]))

        # print(convertdouble)
        try:
            self.st.deleteBond(convertdouble[0], convertdouble[1])
            self.st.addBonds([(convertdouble[0], convertdouble[1], 1)])
        except:
            print("could not delete double")

        ## bond pillars to form block
        if self.iblock != 0:
            self.st.deleteBond(deletebond_top[0], deletebond_top[1])
            self.st.deleteBond(deletebond_bottom[0], deletebond_bottom[1])
            self.st.addBonds([(connectme_out[0], connectme_mid[1], 1),
                         (connectme_out[1], connectme_mid[0], 1)])
        else:
            self.st.deleteBond(deletebond_top[0], deletebond_top[1])
            self.st.deleteBond(deletebond_bottom[0], deletebond_bottom[1])
            self.st.addBonds([(connectme_out[0], connectme_mid[0], 1),
                         (connectme_out[1], connectme_mid[1], 1)])

        self.st.deleteAtoms(deleteme)

        # print('checking bonds')
        # print(check1, check2)

        # for i in range(len(check1)):
            # print(self.st.areBound(check1[i], check)

    def getStructure(self):
        return self.st
