import random
from schrodinger import structure
from schrodinger.infra import mm
from schrodinger.structutils import analyze, build

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

def process_sequence(line):

    st = None
    grow_atoms1=[]
    grow_atoms2=[]

    for i, c in enumerate(line):
        fragname = _pdb_lookup[c]
        if i==0:
            st = get_fragment_structure(fragname)
        else:
            (grow_atoms1, grow_atoms2) = grow_fragment(st, grow_atoms1, grow_atoms2, "dna2b", fragname)

    return st

'''
# main here
def main():
    ct = process_sequence("ACGTACAAATTT")
    ct.write("DNA_model.pdb")
    return True


if __name__ == "__main__":
    main()
'''
