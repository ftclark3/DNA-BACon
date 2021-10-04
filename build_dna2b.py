'''
NOTE TO THE USER: This module uses functions from the schrodinger.infra.mm
module, a wrapper which is not documented well. Comments explaining
the functions get_fragment_structure and grow_fragment can be found
in process_sequence.
'''

from schrodinger.structutils import analyze, build
from schrodinger import structure
from schrodinger.infra import mm
import random

# set the base pairing rules
_pdb_codes = [
  ("A",    "A:T"),
  ("C",    "C:G"),
  ("G",    "G:C"),
  ("T",    "T:A")
]
_pdb_lookup = dict(_pdb_codes)

# place the fragment
def get_fragment_structure(fragname):

    ## creates fragment
    frag_handle = mm.mmfrag_new("dna2b")

    first = True
    st = None
    while True:
        ## check if first atom(?)
        if first:
            which = mm.MM_FIRST
            first = False
        else:
            which = mm.MM_NEXT

        try:
            name, ct_handle = mm.mmfrag_get_fragment(frag_handle, which)
        except mm.MmException as e:
            if e.rc == mm.MMFRAG_DONE:
                break

        if name == fragname:
            copy_ct = mm.mmct_ct_duplicate(ct_handle)
            st = structure.Structure(copy_ct, True)
    mm.mmfrag_delete(frag_handle)

    if st is None:
        raise ValueError("Unknown fragment: %s" % fragname)

    return st

# grow
def grow_fragment(st, grow_atoms1, grow_atoms2, fraggroup, fragname):

    mmfrag_handle = mm.mmfrag_new(fraggroup)

    mm.mmfrag_set_fragment_name(mmfrag_handle, fragname)

    directions=[]

    d = mm.mmfrag_get_direction(mmfrag_handle, mm.MM_FIRST)
    directions.append(d)

    while True:
        try:
            d = mm.mmfrag_get_direction(mmfrag_handle, mm.MM_NEXT)
        except mm.MmException as e:
            if e.rc == mm.MMFRAG_DONE:
                break  # exit the loop (end of list)
        else:
            directions.append(d)

    conn_to, conn_from, new_from, new_to = mm.mmfrag_get_connections(
        mmfrag_handle, st, mm.MM_FIRST)

    #------
    conn_to2, conn_from2, new_from2, new_to2 = mm.mmfrag_get_connections(
        mmfrag_handle, st, mm.MM_NEXT)


    #------

    grow_atoms1=[new_from, new_from2]
    grow_atoms2=[new_to, new_to2]

    (new_grow_atom1, new_grow_atom2, renumber_dict) = build.grow_mult(st, mmfrag_handle, grow_atoms1, grow_atoms2)  # problem solved!!!

    return new_grow_atom1, new_grow_atom2


# read in one sequence of a chain

def process_sequence(line,keep=None,comp=True):

    st = None
    grow_atoms1=[]
    grow_atoms2=[]

    for i, c in enumerate(line):
        # get the complementary base for the fragment
        try:
            fragname = _pdb_lookup[c]
        except KeyError:
            try:
                c = c.upper()
                fragname = _pdb_lookup[c]
            except KeyError:
                print(f"Error: Unrecognized character {c} in sequence")
        if i==0:
            # get a structure representing the first base and its complementary nucleotide
            st = get_fragment_structure(fragname)
        else:
            # extend the two strands by simultaneously adding main sequence and complementary nucleotides
            (grow_atoms1, grow_atoms2) = grow_fragment(st, grow_atoms1, grow_atoms2, "dna2b", fragname)

    '''
    Everything below here should be clear
    '''

    # if there is a particular complementary sequence to keep, mark these bases
    if keep:
        # reverse the order of the string keep, because Schrodinger's iterator
        # for structure._Residue objects, _ResidueIterator, will always return
        # the complementary strand residues in the opposite order of the main strand
        reverse = []
        for i in range(len(keep)):
            reverse.append(keep[-1-i])

        # make a list of complementary residues so that they can be accessed later in the script
        # this is done because _ResidueIterator does not support indexing of residues
        residues = []
        for residue in st.residue:
            if residue.chain == "B" and residue.pdbres.strip() in ["DA","DC","DG","DT"]:
                residues.append(residue)

        # initialize a list to store resnums for the residues on the complementary strand
        # that should be saved
        numbers = []

        # set up counters and loops to go through complementary residues and mark the resnums
        # of residues that should be kept
        c = 0
        i = 0
        k = 0
        # the idea behind this loop is that we start at the first residue (set by k) in
        # the complementary chain and see if we can find keep anywhere in the strand,
        # and we record the residues that are part of keep
        # next, we start at the second residue (k is changed) and see if we can find the sequence
        # anywhere in the complementary strand, and so on
        while True:
            # break if there are fewer residues left in the strand than the length of keep_sequence
            if len(residues) - k < len(reverse):
                break
            # c and i are used to track where we are in the complementary strand and append
            # each residue.resnum to numbers as necessary
            c = k
            for residue in residues[k:]:
                c+=1
                if reverse[i] in residue.pdbres:
                    i+=1
                    if i == len(reverse):
                        break
                else:
                    i = 0
            if i < len(reverse):
                i = 0
            while i > 0:
                numbers.append(residues[c-i].resnum)
                # rename complementary residues so that they can be identified later
                residues[c-i].pdbres = f"X{residues[c-i].pdbres[2].lower()}"
                i-=1
            k+=1

        # delete all complementary DNA bases not designated to be saved
        deleteme = []
        for residue in st.residue:
            if residue.chain == "B" and residue.resnum not in numbers:
                for atom in residue.atom:
                    deleteme.append(atom.index)
        st.deleteAtoms(deleteme)

    # if no particular sequence was supposed to be kept, then delete the entire
    # complementary strand, or delete none of the complementary strand, depending
    # on the value of comp passed to the function
    elif not comp:
        deleteme = []
        for residue in st.residue:
            if residue.chain == "B":
                for atom in residue.atom:
                    deleteme.append(atom.index)
        st.deleteAtoms(deleteme)

    return st

# main here
def main():
    ct = process_sequence("AAAA") 
    ct.write("4A_full.pdb")
    return True

if __name__ == "__main__":
    main()
