import build
from random import randint
from mdSim import mdSim
from sysBuilder import sysBuilder

# this is an example script that shows how to create a cage model and simulate it
##------------------------------------------------------------------------------------------------------------
# set job variables
#path must exist
path = "/home/marlo/HDD/Maestro/jobs/"
# job settings
# change these to customize job
seed = str(randint(0000, 9999)).rjust(4, "0") # generate random 4 digit seed
time = "100000"
interval = "9.6"
job_settings = {"seed":seed, "time":time, "interval":interval}


# build cage
# main() function returns a list with the savepath and struct name
cage = build.main("bacon") # feed in input file

# run sys builder
sysBuild = sysBuilder(cage, path)
sysBuild.moveMAE()
sysBuild.writeMSJ() # buffer: cubic/15 angstroms, edit MSJ string in sysBuilder to change
sysBuild.run()

# run MD
md = mdJob(cage path)
# moveCMS will wait for the sysBuild job to complete before continuing
md.moveCMS()
md.writeMSJ()
md.writeCFG(job_settings)
md.run()
