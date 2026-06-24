# LiMBS
**LipidMem-BigSMILES (LiMBS): A Standardized Line Notation for Lipid Membrane Simulations**

LiMBS employs a dual-block architecture consisting of a system-level descriptor block and a chemical structure block, enabling independent specification of membrane geometry, composition, environmental conditions, and molecular identity.

## Current Features

- Human-readable membrane notation
- Machine-parseable grammar
- Planar membrane representation
- Vesicle membrane representation
- Atomistic lipid notation
- Coarse-grained lipid notation
- Asymmetric membrane systems
  
## Repository Structure

```text
parser/
   LiMBS INSANE and vesicle parser scripts

```
## Included Tools

-  LiMBSv1_parser.py: Core parser and validator for LiMBS notation
-  LiMBS_insane.py: Converts LiMBS CG planar membrane notation into INSANE-compatible input to produce .gro and .top files.
- LiMBS_vesicle.py:  Converts LiMBS CG vesicle notation into TS2CG input (input.str) and runs PLM + PCG steps.
  
##  No Python dependencies required (Python 3.6+ standard library only).

## External tools:
  - insane.py or insane_peptoid_py3.py  (for planar membranes)
  - TS2CG executable in PATH            (for vesicles)




