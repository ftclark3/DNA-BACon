import random
import sys, numpy
import build_dna2b
import schrodinger
import os
import dev
from math import tan, pi
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
from schrodinger.structutils import transform, analyze, minimize
from sklearn.externals import joblib

class dnaStruct:
    def __init__(self, inputfile):
        self.DNA_seq = []
        self.AddSeqs = []
        self.AddStructs = []
        self.Mg = False
        self.name = "wombat"
        self.Min = False
        self.Add = False
        self.inputFile

        # size calculations
    ##------------------------------------------------------------------------------------------------------------
    ## method to read DNA input file and extract settings/sequence
    def sizeCalc(self):
        self.npillar = len(self.DNA_seq)


    ##------------------------------------------------------------------------------------------------------------
    ## method to read DNA input file and extract settings/sequence
    def fileRreader(self):
        inputfile = 'input/' + self.inputFile
        flag = False

        # open and read input file
        fh = open(inputfile)
        for line in fh:
            temp = line.split()
            if temp[0] == "":
                continue
            if temp[0] == "#":
                if temp[1] == "!Ion":
                    self.Mg=True
                if temp[1] == "!Min":
                    self.Min=True
                if temp[1] == "!Name":
                    try:
                        self.name = temp[2]
                    except:
                        print("You didn't include a name!")
                if temp[1] == "!Add":
                    self.Add = True
                    for i in temp[2:]:
                        if ".pdb" not in i:
                            self.AddSeqs.append(i)
                        else:
                            self.AddStructs.append(i)
            if temp[0] == ">":
                flag = True
                continue
            if temp[0] == "<":
                flag = False
            if flag:
                self.DNA_seq.append(temp)
        fh.close()
