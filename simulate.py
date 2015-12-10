import os
import sys
import json
import subprocess
import shlex
import time
import random
from multiprocessing import Pool

NAMD_COMMAND = '/opt/namd/namd2'

eq_conf = """#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          model1_protein_autopsf.psf
coordinates        model1_protein_autopsf.pdb

set temperature    %(temperature)s
set outputname     %(output)s
set inputname      %(input)s

firsttimestep      0

# Continuing a job from the restart files
if {%(continue)s} {
 binCoordinates     $inputname.restart.coor
# binVelocities      $inputname.restart.vel  ;# remove the "temperature" entry if you use this!
 extendedSystem     $inputname.restart.xsc
}



#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          ../../../../common/par_all27_prot_na.prm
temperature         $temperature


# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              24.0
switching           on
switchdist          20.0
pairlistdist        28.0


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2 
stepspercycle       10


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     1     ;# damping coefficient (gamma) of 1/ps
langevinTemp        $temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens

dielectric	%(dielectric)s


# Output
outputName          $outputname
restartfreq         5000    ;# 5000steps = every 10ps
dcdfreq             5000
xstFreq             1000
outputEnergies      500
outputPressure      500


#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################
if {0} {
  Constraints                    yes
  ConsExp                        2
  ConsRef                        constraints.pdb
  ConsKFile                      constraints.pdb
  ConskCol                       B
}

if {1} {
	fixedAtoms	on
	fixedAtomsFile	constraints.pdb
	fixedAtomsCol	B
}

#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Minimization
%(minimize)s
reinitvels          $temperature

run %(time)s"""

analyze_tcl = """
mol new model1_protein_autopsf.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile eq1.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile eq2.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set mol 0
set file [open "rmsd.dat" w]   

# use frame 0 for the reference
set reference [atomselect $mol "protein and resid %(residues)s" frame 0]
# the frame being compared
set compare [atomselect $mol "protein and resid %(residues)s"]

set num_steps [molinfo $mol get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
        # get the correct frame
        $compare frame $frame

        # compute the transformation
        set trans_mat [measure fit $compare $reference]
        # do the alignment
        $compare move $trans_mat
        # compute the RMSD
        set rmsd [measure rmsd $compare $reference]
        # print the RMSD
        puts $file "$frame $rmsd"
}
exit
"""

def runSimulation(parameters):
	time.sleep(float( random.randint(1000,30000)/1000.0))
	os.chdir(parameters['currentDir'])
	if not os.path.exists(parameters['folders'][0]):
		try:
			os.makedirs(parameters['folders'][0])
		except:
			pass
	os.chdir(parameters['folders'][0])
	if not os.path.exists(parameters['folders'][1]):
		try:
			os.makedirs(parameters['folders'][1])
		except:
			pass
	os.chdir(parameters['folders'][1])
	if not os.path.exists(parameters['folders'][2]):
		try:
			os.makedirs(parameters['folders'][2])
		except:
			pass
	os.chdir(parameters['folders'][2])
	iteration = "1"
	for i in range(1000):
		iteration = str(i)
		if not os.path.exists(iteration):
			try:
				os.makedirs(iteration) 
				break
			except:
				pass
	os.chdir(iteration)

	with open('model1_protein.pdb','w') as f2:
		with open('../../../../model1_protein.pdb','r') as f:
			for line in f:
				newline = line
				if 'ATOM' in line:
					residue = int(line.split()[5])
					if residue <= parameters['proteinFull'][1] and residue >= parameters['proteinFull'][0]:
						newline = newline[:56] + "0.00" + newline[60:]
						newline = newline[:62] + "0.00" + newline[66:]
						f2.write(newline)
				else:
					f2.write(line)
	with open('eq1.conf','w') as f:
		f.write(eq_conf % {'dielectric':"100.0",
			'minimize':'minimize            50',
			'continue':'0',
			'input':'eq1',
			'output':'eq1',
			'temperature':'1000',
			'time':str(int(0.03*500000))})
	with open('eq2.conf','w') as f:
		f.write(eq_conf % {'dielectric':str(parameters['dielectric']),
			'minimize':'minimize            50',
			'continue':'1',
			'input':'eq1',
			'output':'eq2',
			'temperature':'310',
			'time':str(int(0.2*500000))})
	with open('analyze.tcl','w') as f:
		f.write(analyze_tcl % {'residues': str(parameters['proteinFree'][0]) + ' to ' + str(parameters['proteinFree'][1])})


	cmd = ['vmd','-dispdev','text','-e','../../../../build.tcl'] 
	print(cmd)
	tstart = time.time()
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "Building PSF"
	out, err = p.communicate()
	with open('build.log','w') as f:
		f.write(out)
	print "finished " + str(time.time()-tstart)
	time.sleep(1)
	with open('constraints.pdb','w') as f2:
		with open('model1_protein_autopsf.pdb','r') as f:
			for line in f:
				newline = line
				if 'ATOM' in line:
					residue = int(line.split()[5])
					if residue <= parameters['proteinFree'][1] and residue >= parameters['proteinFree'][0]:
						f2.write(newline)
					else:
						newline = newline[:56] + "0.00" + newline[60:]
						newline = newline[:62] + "1.00" + newline[66:]
						f2.write(newline)
				else:
					f2.write(line)
	cmd = [NAMD_COMMAND,'+p1',os.getcwd() + '/eq1.conf']
	print(cmd)
	tstart = time.time()
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "Running thermal denaturation"
	out, err = p.communicate()
	with open('eq1.log','w') as f:
		f.write(out)
	print "finished " + str(time.time()-tstart)

	cmd = [NAMD_COMMAND,'+p1',os.getcwd() + '/eq2.conf']
	print(cmd)
	tstart = time.time()
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "Running folding simulation"
	out, err = p.communicate()
	with open('eq2.log','w') as f:
		f.write(out)
	print "finished " + str(time.time()-tstart)


	cmd = ['vmd','-dispdev','text','-e','analyze.tcl']
	print(cmd)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	print "Running analysis"
	out, err = p.communicate()
	with open('analysis.log','w') as f:
		f.write(out)
	print "finished"


if __name__ == "__main__":
	allParameters = []
	for i in range(80,540,5):
		parameters = {}
		parameters['proteinFull'] = [1,i]
		parameters['proteinFree'] = [i-15,i]
		parameters['dielectric'] = 0.75
		parameters['folders'] = ["1-"+ str(i)]
		parameters['folders'].append(str(i-15) + "-" + str(i))
		parameters['folders'].append(str(float(parameters['dielectric'])))
		parameters['currentDir'] = os.getcwd()
		for j in range(2):
			allParameters.append(parameters)
	#runSimulation(parameters)
	p = Pool(62)
	p.map(runSimulation,allParameters)


