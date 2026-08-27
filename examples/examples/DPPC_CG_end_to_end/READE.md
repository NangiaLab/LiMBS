#LiMBS1.0 End-to-End Example: DPPC Coarse-grained Bilayer

##Purpose
This directory provides a complete example showing how a LiMBS1.0 membrane specification can be used in an end-to-end coarse-grained membrane workflow.
The example follows the sequence:
LiMBS specification:
-LiMBS validation
-LiMBS-INSANE interface
-INSANE embrane construction
-MARTINI3/GROMACS 
-Energy minimization
-Equilibration
-Production molecular dynamics
LiMBS itself provides the standardized membrane-system representation, validation, and translation of the specification into inputs for the external INSANE membrane-building program. Coordinates generation, force-field assignment, equilibration, and molecular-dynamics simulation are performed by external software.

##Example System
The example is a symmetric coarse-grained DPPC planar bilayer with DPPC in both leaflets, solvated with MARTINI water and NaCl.
The complete LiMBS specification is provided in: 'DPPC_CG.txt'
##Files
1. DPPC_CG.txt
2. membrane.gro
3. membrane.top
4. em.mdp
5. nvt.mdp
6. npt.mdp
7. md.mdp
LiMBS does not generate the additional force-field files required by MARTINI3; obtain them separately from the appropriate MARTINI distribution.
## Software requirements
The LiMBS1.0 reference implementation used for this example was tested with:
Python 3.11.15
GROMACS 2023.2
INSANE (insane.py executed as a Python script)
MARTINI3 coarse-grained force field
The LiMBS parser and LiMBS-INSANE interface require Python 3. INSANE and GROMACS are external programs and are not distributed or implemented by LiMBS.
