# LiMBS
**LipidMem-BigSMILES (LiMBS): A Standardized line notation for Lipid Membrane Simulations**
LiMBS employs a dual-block architecture consisting of a system-level descriptor block and a chemical structure block, enabling independent specification of membrane geometry, composition, environmental conditions, and molecular identity.

## Current Features

- Human-readable membrane notation
- Machine-parseable grammar
- Planar membrane representation
- Vesicle membrane representation
- Atomistic lipid notation
- Coarse-grained lipid notation
- Protein-embedded membrane systems
- Asymmetric membrane systems
- Prototype web interface
- Prototype parser workflows for INSANE and TS2CG
## Repository Structure

```text
parser/
    Parser-related scripts

web/
    Flask-based web interface for generating LiMBS notation
```
## Included Tools

- Limbs_insane.py: converts LiMBS notation into INSANE-compatible membrane builder input.
- LiMBS_vesicle.py: converts LiMBS notation into TS2CG vesicle input.
- Flask web interface for generating LiMBS notation interactively.
