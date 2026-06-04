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
- Protein-embedded membrane systems
- Asymmetric membrane systems
- Prototype web interface
- Prototype parser workflows for INSANE and TS2CG
## Repository Structure

```text
parser/
   LiMBS INSANE and vesicle parser scripts

web/
    app.py
    lipid_registry.py
    cg_templates.py

web/templates/
    form.html
```
## Included Tools

- LiMBS_insane.py: Converts LiMBS notation into INSANE-compatible membrane builder input.
- LiMBS_vesicle.py: Converts LiMBS notation into TS2CG vesicle input.
- Web interface: Generates LiMBS notation interactively through a Flask-based web application.

  ## Installation

### Clone the repository

```bash
git clone https://github.com/NangiaLab/LiMBS.git
cd LiMBS
```

### Install dependencies

```bash
pip install -r requirements.txt
```

### Launch the web interface

```bash
cd web
python app.py
```

Open your browser and navigate to:

```text
http://127.0.0.1:5000
```

