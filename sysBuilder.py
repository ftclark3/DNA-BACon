import sys
import subprocess
from os import mkdir, chdir, getcwd
from shutil import copyfile, rmtree

## syself.stem builder class
# allows user to run syself.stem builder on generated model

class sysBuilder:
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

    # method to move mae to job folder
    def moveMAE(self):
        copyfile(self.st[0], self.dst)

    # method to write .msj
    def writeMSJ(self):
        # write msj to job folder
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

    # method to run system builder
    def run(self):
        subprocess.Popen(self.bashCommand, shell=True, cwd=self.path + self.st[1] + "_setup")
