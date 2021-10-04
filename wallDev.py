from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize, measure
from schrodinger import structure
from math import sin, tan, pi, atan
import build_dna2b
import schrodinger
import sys, numpy
import random
import toml
import os
import copy

''' Tools to build wall using input object and pillar number '''

# dictionary to help give a unique number to each residue that is created
resnum_counter = {"A":0,"B":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"J":0,"K":0,"L":0}
# dictionary to help find the 3' oxygen on the end of each strand to connect to the linker
internal_atom_names = {"sense":{"A":"O53","C":"O55","G":"O55","T":"O53"},"antisense":{"A":"O74","C":"O80","G":"O76","T":"O76"}}

# Class to hold information on columns (i.e. DNA strands)
class column:
    # constructor
    def __init__(self, dnaSeqs, ipillar, npillar, columnNum):
        self._dnaSeqs = dnaSeqs
        self._ipillar = ipillar
        self._npillar = npillar
        self._chainName = chr(65+ipillar)
        self._columnNum = columnNum

        # Print functions for debugging
        # print("DNA: " + str(self._dnaSeqs))
        # print("Ipillar: " + str(self._ipillar))
        # print("Column: " + str(self._columnNum))

    ##------------------------------------------------------------------------------------------------------------
    # method to make column and move it
    def makeColumn(self, transform_matrix, newpos, rotation=None, keep=None, comp=True):
        ## build column
        ct = build_dna2b.process_sequence(self._dnaSeqs[self._columnNum], keep=keep, comp=comp)
        ## calculate atom index
        atom_index_list = range(1, ct.atom_total+1)
        ## rotate into posn
        transform.transform_structure(ct, transform_matrix, atom_index_list)
        ## rotate if necessary
        if rotation is not None:
            schrodinger.structutils.transform.rotate_structure(ct, x_angle=rotation[0], y_angle=rotation[1], z_angle=rotation[2], rot_center=rotation[3])
        # translate into position, using the centroid of the main (sense) DNA sequence
        # this improves consistency by eliminating complementary DNA, which may be
        # distributed unevenly along the main strand, from the calculation
        temp = copy.deepcopy(ct)
        deleteme = []
        for chain in temp.chain:
            if chain.name == "B":
                for atom in chain.atom:
                    deleteme.append(atom.index)
        temp.deleteAtoms(deleteme)
        centroid = transform.get_centroid(temp, atom_list=[atom.index for atom in temp.atom if atom.element=="P"])
        (x, y, z, q) = newpos - centroid
        del temp
        transform.translate_structure(ct, x, y, z, atom_index_list)

        ## rename chains
        self._chainsResnums(ct)
        return ct

    ##------------------------------------------------------------------------------------------------------------
    # method to rename chains and renumber resnums according to columnNum
    def _chainsResnums(self, ct):

        # if first column
        # OLD CODE if self._columnNum == 1:
        if self._columnNum == 0 or self._columnNum == 2:
            # delete unwanted HXL/POT atoms
            deleteme = []
            for chain in ct.chain:
                for residue in chain.residue:
                    if residue.pdbres.strip() == "POT":
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    elif residue.pdbres.strip() == "HXL":
                        for atom in residue.atom:
                            if atom.element == "H": # doesn't seem to be capable of permanently deleting H
                                deleteme.append(atom.index) # (see comments below)
            ct.deleteAtoms(deleteme)
            # rename chains and renumber residues
            for chain in ct.chain:
                if len(ct.chain) == 1:
                    #chain.name = chr(65+(self._npillar*2)+self._ipillar)
                    chain.name = self._chainName
                elif len(ct.chain) == 2:
                    if chain.name == "A":
                        chain.name = self._chainName
                    elif chain.name == "B":
                        chain.name = chr(65+(self._npillar*2)+self._ipillar)
                    else:
                        print("Unexpected else block: expected chain.name == 'A' or chain.name == 'B'")
                        raise SystemExit()
                else:
                    print("Unexpected else block: expected len(ct.chain)--basically the ct returned by build_dna2b--to have either 1 or 2 chains")
                    raise SystemExit()
                for residue in chain.residue:
                    #global resnum_counter
                    resnum_counter[chain.name]+=1
                    residue.resnum = resnum_counter[chain.name]
            # marking atoms
            for chain in ct.chain:
                for residue in chain.residue:
                    if residue.resnum == resnum_counter[chain.name]-len(chain.residue)+1 and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]: # should be HXL, which has just 1 atom
                        for atom in residue.atom:
                            if atom.element == "O":
                                atom.element = "Kr"
                            else: # sometimes _chainIterator adds deleted hydrogens back
                                pass # maybe try to delete the hydrogen if it's causing problems
                    elif residue.resnum == resnum_counter[chain.name] and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]:
                        for atom in residue.atom:
                            if atom.name == internal_atom_names["sense"][atom.getResidue().pdbres[2]]:
                                atom.element = "Kr"

        elif self._columnNum == 1:
            # deleting HXL and POT
            deleteme = []
            for chain in ct.chain:
                for residue in chain.residue:
                    if residue.pdbres.strip() == "POT":
                        for atom in residue.atom:
                            deleteme.append(atom.index)
                    elif residue.pdbres.strip() == "HXL":
                        for atom in residue.atom:
                            if atom.element == "H": # doesn't seem to be capable of permanently deleting H
                                deleteme.append(atom.index) # (see comments below)
            ct.deleteAtoms(deleteme)
            # renaming chains, renumbering residues
            first = True
            for chain in ct.chain:
                if first:
                    chain.name = self._chainName
                    first = False
                else:
                    #chain.name = chr(65+self._npillar-1)
                    chain.name = chr(65+self._npillar)

                for residue in chain.residue:
                    resnum_counter[chain.name]+=1
                    residue.resnum = resnum_counter[chain.name]
            # marking atoms
            for chain in ct.chain:
                if chain.name == chr(65+self._npillar): # what used to be "B"
                    for residue in chain.residue:
                        if residue.resnum == resnum_counter[chain.name]-len(chain.residue)+1 and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]:  # only 1 atom left after HXL hydrogen is deleted
                            for atom in residue.atom:
                                atom.element = "Xe"
                        elif residue.resnum == resnum_counter[chain.name] and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]:
                            for atom in residue.atom:
                                if atom.name == internal_atom_names["antisense"][atom.getResidue().pdbres[2]]: # "O74"
                                    atom.element = "Xe"
                                else:
                                    pass
                        else:
                            pass # we only want to mark atoms in the terminal residues
                elif chain.name == self._chainName:
                    for residue in chain.residue:
                        if residue.resnum == resnum_counter[chain.name]-len(chain.residue)+1 and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]:
                            for atom in residue.atom:
                                if atom.element == "O":
                                    atom.element = "Kr"
                                else: # sometimes _chainIterator adds deleted hydrogens back
                                    pass # maybe try to delete the hydrogen if it's causing problems
                        elif residue.resnum == resnum_counter[chain.name] and residue.pdbres not in [" Xa "," Xc "," Xg "," Xt "]:
                            for atom in residue.atom:
                                if atom.name == internal_atom_names["sense"][atom.getResidue().pdbres[2]]:
                                    atom.element = "Kr"
                                else:
                                    pass
                        else:
                            pass # we only want to mark atoms in the terminal residues
                else:
                    print(f"Unexpected else block: expected chain.name to be either {self._chainName} (self._chainName) or {chr(65+self._npillar)} (chr(65+self._npillar))")

        else:
            print("Unexpected else block: expected columnNum to be either 0, 1, or 2")
            raise SystemExit()


##------------------------------------------------------------------------------------------------------------
# class to hold linker information
class linker:
    # constructor
    def __init__(self, end, start, center, linkerfile, linkerNum):
        self._endResnum = end
        self._startResnum = start
        self._centerAtom = center
        self._linkerNum = linkerNum
        self._link = next(structure.StructureReader(linkerfile))

    ##------------------------------------------------------------------------------------------------------------
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

    ##------------------------------------------------------------------------------------------------------------
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
def buildColumns(seqs, ipillar, npillar, width, width2, keep, next_c1seq, next_c2seq, next_c3seq, num_walls):

    # first column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([1, 0, 0]))
    newpos = numpy.array([0.0, 0.5*width2, 0.0, 0.0])
    c1 = column(seqs, ipillar, npillar, 0)
    c1 = c1.makeColumn(transform_matrix, newpos, keep=keep, comp=False)

    # third column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([-1, 0, 0]))
    newpos = numpy.array([0.0, -0.5*width2, 0.0, 0.0])
    c3 = column(seqs, ipillar, npillar, 2)
    c3 = c3.makeColumn(transform_matrix, newpos, keep=keep, comp=False)

    # second (middle) column
    transform_matrix = transform.get_alignment_matrix(numpy.array([0, 0, 1]), numpy.array([0, -1, 0]))
    newpos = numpy.array([0.5*width, 0.0, 0.0, 0.0])
    rotation = [0, 0.5*pi,0, None]
    c2 = column(seqs, ipillar, npillar, 1)
    c2 = c2.makeColumn(transform_matrix, newpos, rotation)
    initial_top = copy.deepcopy(c2.atom[16].xyz)
    initial_bottom = copy.deepcopy(c2.atom[len(c2.atom)-16].xyz)

    # add in some additional transformations to account for the variety of
    # wall geometries needed to accomodate irregular cages
    # this is an imprecise heuristic which leaves the finer details to
    # schrodinger's minimize function and multisim/desmond

    # first, rotate c2 to account for different top (c1) and bottom (c3) sequence lengths
    # avoid counting complementary residues
    comp_list = ["Xa","Xc","Xg","Xt"]
    first = [residue.pdbres for residue in c1.residue if residue.pdbres.strip() not in comp_list]
    third = [residue.pdbres for residue in c3.residue if residue.pdbres.strip() not in comp_list]
    if len(third) > len(first):
        angle = atan(len(third)/len(c2.residue))
    elif len(third) < len(first):
        overhang = (len(first) - len(third))/2
        angle = -1 * atan(overhang/len(c2.residue))
    else: # len(third) == len(first)
        angle = 0
    old_bottom_coords = copy.deepcopy(c2.atom[len(c2.atom)-16].xyz) # recording for later
    transform.rotate_structure(c2,-angle,0,angle,c2.atom[16].xyz)
    rotation_movement = ((len(c2.residue)-1)/2) * sin(angle)
    new_bottom_coords = c2.atom[len(c2.atom)-16].xyz # recording for later

    # next, make sure that the c2s are spaced correctly
    if len(first) < len(third): # we want to translate out, not in
        coords_change = [new_bottom_coords[0]-old_bottom_coords[0],new_bottom_coords[1]-old_bottom_coords[1],new_bottom_coords[2]-old_bottom_coords[2]]
        const = 1/2
        transform.translate_structure(c2,const*coords_change[0],const*coords_change[1],const*coords_change[2])

    # third, make sure that the top and bottom segments of the wall
    # are angled correctly to connect the double-stranded (c2) segments

    # this requires some information about the next c2 in the cage, which
    # is part of a different wall
    # setting up a temporary copy of c2 helps with this
    temp = copy.deepcopy(c2)
    if len(next_c3seq) > len(next_c1seq):
        angle = atan(len(next_c3seq)/len(next_c2seq))  # + 1/6 just a constant to give some extra clearance
    elif len(next_c3seq) < len(next_c1seq):
        overhang = (len(next_c1seq) - len(next_c3seq))/2
        angle = -1 * atan(overhang/len(next_c2seq))
    else: # len(next_c3seq) == len(next_c1seq)
        angle = 0
    transform.rotate_structure(temp,-angle,0,angle,temp.atom[16].xyz)
    # now we get a vector aligned with the next c2 in the cage
    v_temp = [temp.atom[len(temp.atom)-16].x-temp.atom[16].x,temp.atom[len(temp.atom)-16].y-temp.atom[16].y,temp.atom[len(temp.atom)-16].z-temp.atom[16].z]
    # and now a vector aligned with the c2 of the current wall
    v_c2 = [c2.atom[len(c2.atom)-16].x-c2.atom[16].x,c2.atom[len(c2.atom)-16].y-c2.atom[16].y,c2.atom[len(c2.atom)-16].z-c2.atom[16].z]
    # now we scale v_temp to the length of the ACTUAL next c2 (not the copy of the current one)
    scalar = len(next_c2seq)/((len(c2.residue)/2)-1)
    if scalar <= 1:
        scalar = 1/scalar
    else:
        scalar = -scalar
    v_temp = [scalar*v_temp[0],scalar*v_temp[1],scalar*v_temp[2]]
    # this allows us to get a vector representing the movement along the c2 axis
    # that is required for c1 and c3
    if scalar == 1:
        v_diff = [0,0,0] # no adjustment necessary
    else: # scalar > 0, scalar != 1
        # divide by 4 to dampen the effects of different c2 angles on c1 and c3 angles
        v_diff = [(v_temp[0]-v_c2[0])/4,(v_temp[1]-v_c2[1])/4,(v_temp[2]-v_c2[2])/4]
        # rotate c3 to point toward the bottom of the next c2, which may have been rotated
        transform.rotate_structure(c3, 0, angle, 0, c3.atom[16].xyz)
    # adding v_diff to vectors representing c1 and c3 gives vectors representing
    # the new orientations of c1 and c3
    v_c1 = [c1.atom[16].x - c1.atom[len(c1.atom)-16].x, c1.atom[16].y - c1.atom[len(c1.atom)-16].y, c1.atom[16].z - c1.atom[len(c1.atom)-16].z]
    v_sum1 = [v_c1[0]+v_diff[0],v_c1[1]+v_diff[1],v_c1[2]+v_diff[2]]
    v_c3 = [c3.atom[len(c3.atom)-16].x - c3.atom[16].x, c3.atom[len(c3.atom)-16].y - c3.atom[16].y , c3.atom[len(c3.atom)-16].z - c3.atom[16].z]
    v_sum3 = [v_c3[0]-v_diff[0],v_c3[1]-v_diff[1],v_c3[2]-v_diff[2]]
    # schrodinger.structutils.transform provides functions for doing the rest:
    matrix1 = transform.get_alignment_matrix(v_c1,v_sum1)
    matrix3 = transform.get_alignment_matrix(v_c3,v_sum3)
    transform.transform_structure(c1,matrix1)
    transform.transform_structure(c3,matrix3)

    # fourth, we center c1 and c3 in the space between the current c2 and
    # the c2 of the next wall by translating
    c_m1 = (len(first)-1)/2
    c_m3 = (len(first)-1)/2
    midpoint1 = [c2.atom[16].x-4.96*(c_m1),c2.atom[16].y,c2.atom[16].z] # 4.96 is the width of a nucleotide
    midpoint3 = [c2.atom[len(c2.atom)-16].x-4.96*(rotation_movement+c_m3),c2.atom[len(c2.atom)-16].y,c2.atom[len(c2.atom)-16].z]
    list1 = [atom.index for atom in c1.atom if type(atom.chain)==type("string")]
    list3 = [atom.index for atom in c3.atom if type(atom.chain)==type("string")]
    transform.translate_to_origin(c1,list1)
    transform.translate_to_origin(c3,list3)
    transform.translate_structure(c1,midpoint1[0],midpoint1[1],midpoint1[2])
    transform.translate_structure(c3,midpoint3[0],midpoint3[1],midpoint3[2])

    # the transformations are now finished

    # but the coordinates where the linkers need to be placed should be
    # recorded and returned to buildWall
    v_l12 = [c2.atom[16].x-c1.atom[len(c1.atom)-16].x,c2.atom[16].y-c1.atom[len(c1.atom)-16].y,c2.atom[16].z-c1.atom[len(c1.atom)-16].z]
    v_l34 = [c3.atom[16].x-c2.atom[len(c2.atom)-16].x,c3.atom[16].y-c2.atom[len(c2.atom)-16].y,c3.atom[16].z-c2.atom[len(c2.atom)-16].z]
    coords_l2 = [c1.atom[len(c1.atom)-16].x+(v_l12[0]*1/2),c1.atom[len(c1.atom)-16].y+(v_l12[1]*1/2),c1.atom[len(c1.atom)-16].z+(v_l12[2]*1/2)]
    for c, number in enumerate(coords_l2): # just giving the linkers some extra clearance
        if number > 0:
            coords_l2[c] = number + 10
        elif number < 0:
            coords_l2[c] = number - 10
        else: # number == 0; unexpected case
            raise Exception("Unexpected Else Block: Origin is not the center of the cage")
    coords_l3 = [c2.atom[len(c2.atom)-16].x+(v_l34[0]/2),c2.atom[len(c2.atom)-16].y+(v_l34[1]/2),c2.atom[len(c2.atom)-16].z+(v_l34[2]/2)]
    for c, number in enumerate(coords_l3):
        if number > 0:
            coords_l3[c] = number + 10
        elif number < 0:
            coords_l3[c] = number - 10
        else: # number == 0; unexpected case
            raise Exception("Unexpected Else Block: Origin is not the center of the cage")
    ratio = len(first)/len(third)
    if ratio >= 1:
        v_pos = copy.deepcopy(c1.atom[16].xyz)
        for c, number in enumerate(v_pos):
            if number > 0:
                v_pos[c] = number + 10
            elif number < 0:
                v_pos[c] = number - 10
            else: # number == 0 ; unexpected case
                raise Exception("Unexpected Else Block: Origin is not the center of the cage")
        coords_l1 = copy.deepcopy(v_pos)
        scalar = (3.3*(len(first)+len(third))/2)/((v_c3[0]**2 + v_c3[1]**2 + v_c3[2]**2)**0.5)
        coords_l4 = [scalar*v_c3[0], scalar*v_c3[1], scalar*v_c3[2]]
        coords_l4 = [c3.atom[16].x + coords_l4[0], c3.atom[16].y + coords_l4[1], c3.atom[16].z + coords_l4[2]]
        for c, number in enumerate(coords_l4):
            if number > 0:
                coords_l4[c] = number + 10
            elif number < 0:
                coords_l4[c] = number - 10
            else: # number == 0 ; unexpected case
                raise Exception("Unexpected Else Block: Origin is not the center of the cage")
    else: # ratio < 1
        v_pos = copy.deepcopy(c3.atom[len(c3.atom)-16].xyz)
        for c, number in enumerate(v_pos):
            if number > 0:
                v_pos[c] = number + 10
            elif number < 0:
                v_pos[c] = number - 10
            else: # c == 0 ; unexpected case
                raise Exception("Unexpected Else Block: Origin is not the center of the cage")
        coords_l4 = copy.deepcopy(v_pos)
        scalar = (3.3*(len(first)+len(third))/2)/((v_c1[0]**2 + v_c1[1]**2 + v_c1[2]**2)**0.5)
        coords_l1 = [scalar*v_c1[0],scalar*v_c1[1],scalar*v_c1[2]]
        coords_l1 = [c1.atom[len(c1.atom)-16].x + coords_l1[0], c1.atom[len(c1.atom)-16].y + coords_l1[1], c1.atom[len(c1.atom)-16].z + coords_l1[2]]
        for c, number in enumerate(coords_l1):
            if number > 0:
                coords_l1[c] = number + 10
            elif number < 0:
                coords_l1[c] = number - 10
            else: # number == 0 ; unexpected case
                raise Exception("Unexpected Else Block: Origin is not the center of the cage")

    big_list = [coords_l1,coords_l2,coords_l3,coords_l4]

    # and the columns need to be merged into one st
    columns = c1.merge(c2, copy_props=True)
    columns = columns.merge(c3, copy_props=True)

    return columns, big_list

##------------------------------------------------------------------------------------------------------------
# function to build linkers of the wall
def buildLinkers(endResnum, startResnum, centerAtom, linkerFile, ipillar, width, width2, linker_coords):
    # first linker (top left)
    l1 = linker(endResnum, startResnum, centerAtom, linkerFile, 1)
    if linkerFile == "linkers/4T.pdb":
        l1 = l1.makeLinker(numpy.array([-0.5*width, 0.5*width2, 0.25*width]), [0, 0, 0], ipillar) #pi/2
    else:
        l1 = l1.makeLinker(numpy.array([-0.45*width, 0.5*width2, 0.15*width]), [0, 0, pi], ipillar)

    # second linker (top right)
    l2 = linker(endResnum, startResnum, centerAtom, linkerFile, 2)
    if linkerFile == "linkers/4T.pdb":
        l2 = l2.makeLinker(numpy.array([0.45*width, 0.52*width2, 0.2*width]), [pi, pi, 0], ipillar)
    else:
        l2 = l2.makeLinker(numpy.array([0.45*width, 0.52*width2, 0.2*width]), [0, 0, 0], ipillar)

    # third linker (bottom right)
    l3 = linker(endResnum, startResnum, centerAtom, linkerFile, 3)
    if linkerFile =="linkers/4T.pdb":
        l3 = l3.makeLinker(numpy.array([0.45*width, -0.5*width2, 0.15*width]), [pi, 0, 0], ipillar)
    else:
        l3 = l3.makeLinker(numpy.array([0.45*width, -0.5*width2, 0.15*width]), [0.25*pi, 0, 0], ipillar)

    # fourth linker (bottom left)
    l4 = linker(endResnum, startResnum, centerAtom, linkerFile, 4)
    if linkerFile == "linkers/4T.pdb":
        l4 = l4.makeLinker(numpy.array([-0.4*width, -0.5*width2, 0.25*width]), [0, pi, 0], ipillar)
    else:
        l4 = l4.makeLinker(numpy.array([-0.4*width, -0.5*width2, 0.15*width]), [0, 0, 0], ipillar)

    # translate linkers into position
    transform.translate_structure(l1, linker_coords[0][0]-l1.atom[len(l1.atom)//2].x, linker_coords[0][1]-l1.atom[len(l1.atom)//2].y, linker_coords[0][2]-l1.atom[len(l1.atom)//2].z)
    transform.translate_structure(l2, linker_coords[1][0]-l2.atom[len(l2.atom)//2].x, linker_coords[1][1]-l2.atom[len(l2.atom)//2].y, linker_coords[1][2]-l2.atom[len(l2.atom)//2].z)
    transform.translate_structure(l3, linker_coords[2][0]-l3.atom[len(l3.atom)//2].x, linker_coords[2][1]-l3.atom[len(l3.atom)//2].y, linker_coords[2][2]-l3.atom[len(l3.atom)//2].z)
    transform.translate_structure(l4, linker_coords[3][0]-l4.atom[len(l4.atom)//2].x, linker_coords[3][1]-l4.atom[len(l4.atom)//2].y, linker_coords[3][2]-l4.atom[len(l4.atom)//2].z)

    # merge linkers
    linkers = l1.merge(l2, copy_props=True)
    linkers = linkers.merge(l3, copy_props=True)
    linkers = linkers.merge(l4, copy_props=True)

    return linkers

##------------------------------------------------------------------------------------------------------------
# function to take build wall out of linkers and columns
def buildWall(struct, ipillar, npillar, width, width2, keep):
    # import linker information
    # create linker object to measure size
    linkerInfo = toml.load("linkers/linkerInformation.toml")
    linker = linkerInfo[struct["linker"]][0]

    # extract info from wall input
    seqs = struct["sequence"][ipillar]
    if ipillar + 1 == len(struct["sequence"]):
        next_c1seq = struct["sequence"][0][0]
        next_c2seq = struct["sequence"][0][1]
        next_c3seq = struct["sequence"][0][2]
    else:
        next_c1seq = struct["sequence"][ipillar+1][0]
        next_c2seq = struct["sequence"][ipillar+1][1]
        next_c3seq = struct["sequence"][ipillar+1][2]
    endResnum = linker["linkerEResnum"]
    startResnum = linker["linkerSResnum"]
    centerAtom = linker["linkerCenter"]
    linkerFile = linker["file"]
    linkerEName = linker["linkerEName"]
    linkerSName = linker["linkerSName"]
    linkerEResnum = linker["linkerEResnum"]
    linkerSResnum = linker["linkerSResnum"]

    # measurement variables
    chainName = chr(65+ipillar)
    chainLength = len(seqs[0])
    chainLength2 = len(seqs[1])
    chainLength3 = len(seqs[2])
    end = linkerEResnum * 10 + ipillar
    start = linkerSResnum * 10 + ipillar

    # generate columns and linkers
    (columns, linker_coords) = buildColumns(seqs, ipillar, npillar, width, width2, keep, next_c1seq, next_c2seq, next_c3seq, len(struct["sequence"]))
    linkers = buildLinkers(endResnum, startResnum, centerAtom, linkerFile, ipillar, width, width2, linker_coords)
    wall = columns.merge(linkers, copy_props=True)

    # bond columns to linkers
    ## initialize bonding lists
    deleteMe = []
    connectMe = []
    linkMeS = []
    linkMeE = []
    i = 1


    ## bond or delete residues
    for residue in wall.residue:

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

    for atom in wall.atom:
        if atom.element == "Kr":
            atom.element = "O"
            connectMe.append(atom.index)


    ## add bonds dont use linkMe 1/6
    ## should be c8 to o5s
    wall.addBonds([(linkMeS[2], connectMe[4], 1), # c12 o5
                 (linkMeS[1], connectMe[2], 1), # c8 o3
                 (linkMeE[3], connectMe[5], 1), # c12 o5
                 (linkMeS[0], connectMe[0], 1), # c8 o3
                 (linkMeE[1], connectMe[1], 1), # c12 o5
                 (linkMeE[2], connectMe[3], 1)
                 ]) # c8 o3

    wall.write(f"models/cages/wall{ipillar}.pdb")

    # return wall structure
    return wall
