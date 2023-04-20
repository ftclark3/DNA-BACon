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

##------------------------------------------------------------------------------------------------------------
## DNA structure class
class dnaStruct:
    def __init__(self):
        self.dnaSeqs = []
        self.mg = False
        self.name = "wombat"
        self.min = False
        self.add = False
        self.addSeqs = []
        self.addStructs = []
        self.linker = "tabmow"

    def setLinker(self):
        ## set linker variables
        if self.linker == 'linkers/HEG.pdb': # HEG LINKER
            self.HEG = True
            self.linkerS = 11
            self.linkerSName = "C8"
            self.linkerSResnum = 1
            self.linkerE = 17
            self.linkerEName = "C12"
            self.linkerEResnum = 1
            self.linkerCenter = "O1"
        elif self.linker == 'linkers/4T.pdb': # 4T linker
            self.HEG = False
            self.linkerS = 5
            self.linkerSName = "C5*"
            self.linkerSResnum = 3
            self.linkerE = 104
            self.linkerEName = "C3*"
            self.linkerEResnum = 6
            self.linkerCenter = "P*"
        elif self.linker == 'linkers/benzene.pdb': # Benzene linker
            self.HEG = False
            self.linkerS = 7
            self.linkerSName = "C7"
            self.linkerSResnum = 1
            self.linkerE = 13
            self.linkerEName = "C13"
            self.linkerEResnum = 1
            self.linkerCenter = "C5"
        else:
            if self.linker == "tabmow":
                self.linker = input("Please enter the linker filename: ")
            self.linkerS = input('linker start atom number: ')
            self.linkerSName = str(input('linker start name: '))
            self.linkerSResnum = input('linker start resnum: ')
            self.linkerE = input('linker end atom number: ')
            self.linkerEName = str(input('linker end name: '))
            self.linkerEResnum = input('linker end resnum :')
            self.linkerCenter = str(input('linker center atom: '))


##------------------------------------------------------------------------------------------------------------
## function to read DNA input file and extract settings/sequence
def fileReader(inputFile):
    # initialize variable
    inputFile = 'input/' + inputFile
    dna = dnaStruct()
    flag = False

    # read input file
    fh = open(inputFile)
    for line in fh:
        temp = line.split()
    #    print(temp)
        if temp[0] == "":
            continue
        if temp[0] == "#":
            if temp[1] == "!Ion":
                dna.mg=True
            if temp[1] == "!Min":
                dna.min=True
            if temp[1] == "!Name":
                try:
                    dna.name = temp[2]
                except:
                    print("You didn't include a name!")
            if temp[1] == "!Add":
                dna.add = True
                for i in temp[2:]:
                    if ".pdb" not in i:
                        dna.addSeqs.append(i)
                    else:
                        dna.addStructs.append(i)
            if temp[1] == "!Linker":
                dna.linker = "linkers/" + temp[2]
        if temp[0] == ">":
            flag = True
            continue
        if temp[0] == "<":
            flag = False
            continue
        if flag:
            dna.dnaSeqs.append(temp)
    fh.close()

    # update linker variables
    dna.setLinker()

    # return dna structure
    return dna
