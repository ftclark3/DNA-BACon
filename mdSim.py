from shutil import copyfile, rmtree
from os import mkdir, chdir, getcwd
import subprocess
import os.path
import time
import sys

''' Class to run MD simulation with Desmond using built cage '''

class mdJob:
    ##------------------------------------------------------------------------------------------------------------
    # constructor
    def __init__(self, st, path):
        self.st = st
        self.path = path
        self.jobName = self.st[1]
        self.bashCommand = "$SCHRODINGER/utilities/multisim -VIEWNAME desmond_molecular_dynamics_gui.MDApp -JOBNAME " + \
                            self.jobName + " -HOST localhost -maxjob 1 -cpu 1 -m " + self.jobName + ".msj -c " + self.jobName + \
                            '.cfg -description "Molecular Dynamics" ' + self.jobName + '.cms -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] ' + \
                            "-PROJ /home/marlo/.schrodinger/tmp/tproj62725a4923 -DISP append -o " + self.jobName+ "-out.cms -ATTACHED"

        # generate new folder
        self.dir = self.path + self.st[1] + "/"
        self.location = self.path + self.st[1] + "_setup/" + self.st[1] + "_setup-out.cms"
        self.dst = self.dir + self.st[1] + ".cms"
        try:
            mkdir(self.dir + "/")
        except:
            rmtree(self.dir + "/")
            mkdir(self.dir + "/")

    ##------------------------------------------------------------------------------------------------------------
    # method to move out.cms to job folder
    def moveCMS(self):
        while not os.path.exists(self.location):
            time.sleep(20)
        copyfile(self.location, self.dst)

    ##------------------------------------------------------------------------------------------------------------
    # method to build msj
    def writeMSJ(self):
        msg = '''
# Desmond standard NPT relaxation protocol
# All times are in the unit of ps.
# Energy is in the unit of kcal/mol.
task {
   task = "desmond:auto"
   set_family = {
      desmond = {
         checkpt.write_last_step = no
      }
   }
}

simulate {
   title       = "Brownian Dynamics NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 100ps"
   annealing   = off
   time        = 100
   timestep    = [0.001 0.001 0.003 ]
   temperature = 10.0
   ensemble = {
      class = "NVT"
      method = "Brownie"
      brownie = {
         delta_max = 0.1
      }
   }
   restrain = {
      atom = "solute_heavy_atom"
      force_constant = 50.0
   }
}

simulate {
   effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   title       = "NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 12ps"
   annealing   = off
   time        = 12
   timestep    = [0.001 0.001 0.003]
   temperature = 10.0
   restrain    = { atom = solute_heavy_atom force_constant = 50.0 }
   ensemble    = {
      class  = NVT
      method = Berendsen
      thermostat.tau = 0.1
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

simulate {
   title       = "NPT, T = 10 K, and restraints on solute heavy atoms, 12ps"
   effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   annealing   = off
   time        = 12
   temperature = 10.0
   restrain    = retain
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 50.0
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

solvate_pocket {
   should_skip = true
   ligand_file = ?
}

simulate {
   title       = "NPT and restraints on solute heavy atoms, 12ps"
   effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                  ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   time        = 12
   restrain    = retain
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 50.0
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

simulate {
   title       = "NPT and no restraints, 24ps"
   effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                  ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   time        = 24
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 2.0
   }

   eneseq.interval   = 0.3
   trajectory.center = solute
}

simulate {
   cfg_file = "'''
        msg += self.jobName
        msg += '''.cfg"
   jobname  = "$MASTERJOBNAME"
   dir      = "."
   compress = ""
}

# Job launching command:
# $SCHRODINGER/utilities/multisim -VIEWNAME desmond_molecular_dynamics_gui.MDApp -JOBNAME g_rec_prism_1 -HOST localhost -maxjob 1 -cpu 1 -m g_rec_prism_1.msj -c g_rec_prism_1.cfg -description "Molecular Dynamics" g_rec_prism_1.cms -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] -PROJ /home/marlo/.schrodinger/tmp/tproj62725a4923 -DISP append -o g_rec_prism_1-out.cms -ATTACHED
        '''

        msj = open(self.path + self.st[1] + "/" + self.st[1] + ".msj", "w")
        msj.write(msg)
        msj.close()

    ##------------------------------------------------------------------------------------------------------------
    # method to build cfg, accepts dictionary of settings
    def writeCFG(self, jobSettings):
        seed = jobSettings["seed"]
        time = jobSettings["time"]
        interval = jobSettings["interval"]

        msg = '''annealing = false
backend = {
}
bigger_rclone = false
checkpt = {
   first = 0.0
   interval = 240.06
   name = "$JOBNAME.cpt"
   write_last_step = true
}
cpu = 1
cutoff_radius = 9.0
elapsed_time = 0.0
energy_group = false
eneseq = {
   first = 0.0
   interval = 1.2
   name = "$JOBNAME$[_replica$REPLICA$].ene"
}
ensemble = {
   barostat = {
      tau = 2.0
   }
   class = NPT
   method = MTK
   thermostat = {
      tau = 1.0
   }
}
glue = solute
maeff_output = {
   first = 0.0
   interval = 120.0
   name = "$JOBNAME$[_replica$REPLICA$]-out.cms"
   periodicfix = true
   trjdir = "$JOBNAME$[_replica$REPLICA$]_trj"
}
meta = false
meta_file = ?
pressure = [1.01325 isotropic ]
randomize_velocity = {
   first = 0.0
   interval = inf
   seed = '''
        msg += seed + '\n'
        msg += '''   temperature = "@*.temperature"
}
restrain = none
simbox = {
   first = 0.0
   interval = 1.2
   name = "$JOBNAME$[_replica$REPLICA$]_simbox.dat"
}
surface_tension = 0.0
taper = false
temperature = [
   [300.0 0 ]
]
time = '''
        msg += time + '\n'
        msg += '''timestep = [0.002 0.002 0.006 ]
trajectory = {
   center = []
   first = 0.0
   format = dtr
   frames_per_file = 250
   interval = '''
        msg += interval + '\n'
        msg += '''   name = "$JOBNAME$[_replica$REPLICA$]_trj"
   periodicfix = true
   write_velocity = false
}
        '''
        msj = open(self.path + self.st[1] + "/" + self.st[1] + ".cfg", "w")
        msj.write(msg)
        msj.close()

    ##------------------------------------------------------------------------------------------------------------
    # method to run job
    def run(self):
        # print(self.path + self.st[1])
        subprocess.Popen(self.bashCommand, shell=True, cwd=self.path + self.st[1] + "/")
