from shutil import copyfile, rmtree
from os import mkdir, chdir, getcwd
import subprocess
import sys

''' Class to run MD system builders on cage models '''

##------------------------------------------------------------------------------------------------------------
class sysBuilder:
    # st is a list that contains a structure and its name
    # path should be the desmond working directory
    # constructor
    def __init__(self, st, path):
        self.st = st
        self.path = path
        self.jobName = self.st[1] + "_setup"
        self.bashCommand = '"${SCHRODINGER}/utilities/multisim" -JOBNAME ' + self.jobName + ' -m ' \
                            + self.jobName + '.msj ' + self.jobName + '.mae -o ' + self.jobName + \
                            "-out.cms -HOST localhost -TMPLAUNCHDIR -ATTACHED"

        # generate new folder
        self.dst = self.path + self.st[1] + "_setup/" + self.st[1] + "_setup" + ".mae"
        try:
            mkdir(self.path + self.st[1] + "_setup/")
        except:
            rmtree(self.path + self.st[1] + "_setup/")
            mkdir(self.path + self.st[1] + "_setup/")

    ##------------------------------------------------------------------------------------------------------------
    # method to move mae to job folder
    def moveMAE(self):
        copyfile(self.st[0], self.dst)

    ##------------------------------------------------------------------------------------------------------------
    # method to write .msj
    def writeMSJ(self):
        msg = '''task {
	task = "desmond:auto"
}

build_geometry {
	add_counterion = {
		 ion = Na
		 number = neutralize_system
	}
	box = {
		 shape = cubic
		 size = 15.0
		 size_type = buffer
	}
	override_forcefield = OPLS3e
	rezero_system = false
	solvent = SPC
}

assign_forcefield {
	forcefield = OPLS3e
}'''
        msj = open(self.path + self.st[1] + "_setup/" + self.st[1] + "_setup.msj", "w")
        msj.write(msg)
        msj.close()

    ##------------------------------------------------------------------------------------------------------------
    # method to run system builder
    def run(self):
        subprocess.Popen(self.bashCommand, shell=True, cwd=self.path + self.st[1] + "_setup")
