import openmm, os
import numpy as np
from openmm import unit, app, Platform
from openmm.app import PDBFile

# from openfold.np.relax.amber_minimize import _add_restraints
# from openfold.np.relax.cleanup import fix_pdb as _fix_pdb
# from src.utils.protein.faspr import faspr_pack
# from src.utils.misc.process import mp_imap_unordered

CPU_COUNT = os.cpu_count()
KJ_PER_MOL = unit.kilojoules_per_mole
KCAL_PER_MOL = unit.kilocalories_per_mole
ANG = unit.angstroms

# Force
KCAL_PER_MOL_NM = unit.kilocalories_per_mole / (unit.nano * unit.meter)
STIFFNESS_UNIT = KCAL_PER_MOL / (ANG ** 2)

def get_forces_for_frame(pdb_path, add_H=False):
    def get_energy_profile(system, simulation, ret={}):
        # Get all forces
        state = simulation.context.getState(getEnergy=True, getForces=True)
        ret['potential_energy'] = state.getPotentialEnergy() / KCAL_PER_MOL

        # Get energy profile
        for i, f in enumerate(system.getForces()):
            state = simulation.context.getState(getEnergy=True, groups={i})
            ret[f.getName()] = state.getPotentialEnergy() / KCAL_PER_MOL

        if 'CustomExternalForce' in ret.keys():
            ret['PotentialEnergy'] = ret['potential_energy'] - ret['CustomExternalForce']
        
        return ret

    def get_atom_forces(simulation, topology):
        atom_forces = simulation.context.getState(getForces=True).getForces() / KCAL_PER_MOL_NM
        atom_forces = np.array(atom_forces)

        ca_forces = []
        atom_idx = 0
        for j, res in enumerate(topology._chains[0].residues()):
            for atom in res.atoms():
                if atom.name =='CA':
                    ca_forces.append(atom_forces[atom_idx])
                atom_idx += 1
        return atom_forces, ca_forces
    platform = openmm.Platform.getPlatformByName('CUDA') # set to GPU mode
    
    forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')
    pdbfile = PDBFile(str(pdb_path))
    modeller = openmm.app.Modeller(pdbfile.topology, pdbfile.positions)

    if add_H:
        modeller.addHydrogens(forcefield, pH=7.0)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        soluteDielectric=1.0, solventDielectric=78.5
    )

    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds) # dummy intergrator. Not used in minimizeEnergy
    simulation = app.Simulation(
        modeller.topology, system, integrator, platform = platform, 
    )
    simulation.context.setPositions(modeller.positions) # update frame coordinates
    total_system_energy = get_energy_profile(system=system, simulation=simulation)

    atom_forces, ca_forces = get_atom_forces(simulation, modeller.topology)
    atom_forces = np.array(atom_forces)
    ca_forces = np.array(ca_forces)

    return atom_forces, ca_forces, total_system_energy['potential_energy']

# usage
atom_forces, ca_forces, system_energy = get_forces_for_frame("16pk/16pk_A.pdb")
print (atom_forces.shape, ca_forces.shape, system_energy)