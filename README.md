## DNA-Bacon
Scripts to (B)uild, (A)ssemble, and (C)onstruct DNA nanocages, junctions and rings. DNA-BACon also includes tools to automate the process of running MD jobs with Desmond.

### Prerequisites
- Schrodinger 2020-2
- Desmond 2020-2 (only needed for MD jobs)
- Numpy 1.19.1
- TOML 0.10.1

### File descriptions
- blockDev.py -- classes used to build blocks of a junction
- build_dev.py -- functions to build cages and additions
- build_dna2b.py -- functions to grow DNA strands
- build.py -- accepts input filename as argument and builds cage or junction accordingly
- mdSim.py -- tools and classes used to run an MD job
- run_example.py -- example MD job script
- seq_writer.py -- generates an input file for cages
- sysBuilder.py -- tools and classes used to run a System Builder job
- wallDev.py -- tools and classes for wall creation


## Getting started
DNA cages and junctions are built according to '.toml' input files, which should be placed in the 'input/' directory. Input files contain DNA sequences, linker instructions, and structure modification instructions. See 'mockup.toml' and the next sections for examples and instructions. In order to build a DNA Structure, construct an input file and use the command:

```
$SCHRODINGER/run build.py input_file.toml
```

### Cage
Below is an example input file for a generator of a four-walled DNA Nanocage. This example includes comments that describe each line. By default, complementary strands on top and bottom strands are deleted for the final, minimized structure. However, an initial version with full double strands is also saved.  

```
[[cage]]
	linker = "benzene" # Filename of pdb to use for linker (for custom linkers, place pdb in linkers/ directory)
	name = "long_seq_benzene_noComp_test" # Filename of final structure
	min = false # Optimize bond lengths using the Schrodinger API (recommended before simulations)
	sequence = [
	["TCGCTGAGTATCGCTGAGTA", "TCCCCGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTATCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA"],
	["TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTATCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA"],
	["TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTATCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA"],
	["TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTATCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA", "TCGCTGAGTATCGCTGAGTA"],
	] # Each line represents a cage wall, at least 3 walls must be present. See .doc for further information.
	keep_sequence = ["TCCCC"] # Unique sequence to specify complementary sequence used to bind cage-decorations
	add_name = ["new.pdb"] # Filename of cage decoration, place in input/ directory
	add_num = [26] # Atom index of the atom that should be attached to the cage
```

### Junction
```
[[junction]]
	name = "junction_test" # name of the final file
	min = false # optimize bond lengths using the Schrodinger API (recommended before simulations)
	sequence = [
	["TCTTCTATACTGGCAAAAAAAAAACAGGATTAGCAGAG"],
	["TCTTCTATACTGGCAAAAAAAAAACAGGATTAGCAGAG"],
	["TCTTCTATACTGGCAAAAAAAAAACAGGATTAGCAGAG"]
	] # each sequence represents one strand of a DNA junction. A minimum of three strands are required
```

## Running MD jobs
DNA-BACon can automatically run customized Desmond MD jobs. See 'run_example.py' for an example.   

### Authors
* Jonathon Ferrell
* Marlo Zorman
* Finn Clark
* Jianing Li
