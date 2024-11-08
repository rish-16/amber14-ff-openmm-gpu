import sys
import MDAnalysis as mda

from openmm.app import *
from openmm import *
from openmm.unit import *

ANG_TO_NM = 0.1

# convert atomic coordinates from trajectory to OpenMM `Quantity` objects (in nm, not angstroms)
uni = mda.Universe('16pk/16pk_A_R3.xtc', format="XTC")
frame_positions = []
for frame in uni.trajectory:
    frame_positions.append(quantity.Quantity(uni.atoms.positions * ANG_TO_NM, nanometer))
    
# instantiate initial crystal structure
pdb = PDBFile('16pk/16pk_A.pdb')
N_residues = pdb.topology.getNumResidues()
N_atoms = pdb.topology.getNumAtoms()

# set up AMBER14 FF with sampling method (downloads potential function automatically)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
unit_cell = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME, # use particle-mesh ewald summation for long-range interactions
    nonbondedCutoff=1.0 * nanometer,
)

# perform at room temperature in implicit waters with 0.004ps step size with minimal friction
integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 0.004*picoseconds)	
simulation = Simulation(
					topology=pdb.topology, 
					system=unit_cell, 
					integrator=integrator
				)

frame_energies = [] # shape: [T, N, 1]
frame_forces = [] # shape: [T, N, 3]

# iterate trajectory to get frame-wise energies and forces
for frame_pos in frame_positions:
    simulation.context.setPositions(frame_pos) # set new coordinates for current frame
    state = simulation.context.getState(
									getEnergy=True,
									getForces=True
								) # computes AMBER14 frame-wise metadata
	
    frame_forces.append(state.getForces())
    frame_energies.append(state.getPotentialEnergy())