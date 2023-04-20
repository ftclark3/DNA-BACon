## DNA-Bacon
Scripts to build, assemble, and construct DNA nanocages and junctions. DNA-BACon also includes tools to automate the process of running MD jobs with Desmond.

### Prerequisites
- Schrodinger API
- Numpy

### File descriptions
- build.py -- accepts arguments and runs functions to build cage
- cage_dev.py -- contains function to assemble cage from walls
- wallDev.py -- tools and classes for wall creation
- inputDev.py -- tools and classes for input
- mdSim.py -- tools and classes used to run an MD job
- sysBuilder.py -- tools and classes used to run a System Builder job
- run_example.py -- example MD job script
- junction_readin.py -- constructs junction using input file
- seq_writer.py -- generates an input file for cages


## Getting started
### Cage
Input files allow for customized cages. Here is an example input file.

```
# !Name trial
# !Linker linker.pdb
# !Ion Mg     
# !Add true   
# !Min
>
TCGCTGAGTA TCAACTGCTCTCAACTGCTC GCAAGTGTGGGCACGCACAC TCAACTGCTCTCAACTGCTC CACAAATCTG
CTATCGGTAG TCAACTGCTCTCAACTGCTC TACTCAGCGACAGATTTGTG TCAACTGCTCTCAACTGCTC CAACTAGCGG
CACTGGTCAG TCAACTGCTCTCAACTGCTC CTACCGATAGCCGCTAGTTG TCAACTGCTCTCAACTGCTC GGTTTGCTGA
CACTGGTCAG TCAACTGCTCTCAACTGCTC CTACCGATAGCCGCTAGTTG TCAACTGCTCTCAACTGCTC GGTTTGCTGA
<
```

- !Name -- (required) provides name of final file
- !Linker -- (required) HEG.pdb / 4T.pdb / custom linker 
- !Ion Mg -- if present Mg will be added to final file in order to balance charge
- !Add -- (cage-specific) if present, user can specify molecule to attach to top and bottom edges (BETA)
- !Min -- if present script will minimize final file

To run, execute: $SCHRODINGER/run build.py (input file)

The input file must be in the "input/" directory
The linker file must be in the "linker/" directory

## Junction
Junction input files are analogous:
```
# !Name test
# !Ion Mg
# !Min
>
TGATCTCTTTGCGAGAACCATCCGTTAAGACCAGAA
AAATGAAGAGACTACGTGAACTAATGAGAAACAGTG
TGATCTCTTTGCGAGAACCATCCGTTAAGACCAGAA
<
```
To run, execute: $SCHRODINGER/run junction_readin.py (input file)
Each line is built into one side

## Running MD jobs
DNA-BACon can automatically run customized Desmond MD jobs. See 'run_example.py' for an example.   

### Authors
* Jianing Li
* Jonathon Ferrell
* Marlo Zorman

### License
I should ask Jonathon about this
