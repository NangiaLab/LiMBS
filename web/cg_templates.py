from pathlib import Path
import re


ITP_DIR = Path("~/MARTINI3/ff").expanduser()


HEADGROUPS = {
    "PC": {
        "atoms": "[N(C)(C)(C)CC:^NC3] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^NC3] [^PO4] [^GL1] [^GL2]",
    },
    "PE": {
        "atoms": "[N(C)(C):^NH3] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^NH3] [^PO4] [^GL1] [^GL2]",
    },
    "PS": {
        "atoms": "[OC(=O)CNC:^CNO] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^CNO] [^PO4] [^GL1] [^GL2]",
    },
    "PG": {
        "atoms": "[OCCOCO:^GL0] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^GL0] [^PO4] [^GL1] [^GL2]",
    },
    "PA": {
        "atoms": "[OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^PO4] [^GL1] [^GL2]",
    },
    "PI": {
        "atoms": "[CCOCO:^C1][CO:^C2][OC:^C3][OC:^C4] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^C1][^C2][^C3][^C4] [^PO4] [^GL1] [^GL2]",
    },
    "DG": {
        "atoms": "[CO:^COH][COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^COH][^GL1] [^GL2]",
    },
    "SM": {
        "atoms": "[N(C)(C)(C)CC:^NC3] [OP(=O)(O)(O):^PO4][NC:^OH1] [C(=O):^AM2]",
        "refs":  "[^NC3] [^PO4] [^OH1][^AM2]",
    },
    "CE": {
        "atoms": "[CO:^COH] [NC:^AM2] [CO:^OH1]",
        "refs":  "[^COH] [^AM2] [^OH1]",
    },
}


def read_itp_beads_for_lipid(lipid_name):
    """
    Read the real Martini bead order from ;@BEADS for a lipid.
    This is safer than guessing bead names from tail length and unsaturation.
    """
    lipid_name = lipid_name.upper()

    if not ITP_DIR.exists():
        return None

    for itp in ITP_DIR.glob("*.itp"):
        text = itp.read_text(errors="ignore")

        insane_matches = list(re.finditer(r";@INSANE[^\n]*", text))

        for m in insane_matches:
            line = m.group(0)

            alname_match = re.search(r"alname=([A-Za-z0-9]+)", line)
            if not alname_match:
                continue

            alname = alname_match.group(1).upper()
            if alname != lipid_name:
                continue

            # Look after the ;@INSANE line for the matching ;@BEADS line
            chunk = text[m.start():m.start() + 2000]
            beads_match = re.search(r";@BEADS\s+(.+)", chunk)

            if beads_match:
                return beads_match.group(1).split()

    return None


def make_tail(length, double_bonds, suffix):
    """
    Fallback tail generator.
    Used only when ;@BEADS is not found in the Martini .itp file.
    """
    bead_count = {
        16: 4,
        18: 4,
        20: 5,
        22: 6,
        24: 6,
    }

    if length not in bead_count:
        raise ValueError(f"Unsupported tail length: {length}")

    n_beads = bead_count[length]

    double_bond_positions = {
        (16, 0): [],
        (16, 1): [3],

        (18, 0): [],
        (18, 1): [3],
        (18, 2): [3, 4],

        (20, 0): [],
        (20, 1): [3],
        (20, 2): [3, 4],
        (20, 3): [3, 4, 5],
        (20, 4): [2, 3, 4, 5],

        (22, 0): [],
        (22, 1): [4],
        (22, 2): [4, 5],
        (22, 4): [3, 4, 5, 6],
        (22, 5): [2, 3, 4, 5, 6],
        (22, 6): [2, 3, 4, 5, 6],

        (24, 0): [],
        (24, 1): [4],
    }

    key = (length, double_bonds)

    if key not in double_bond_positions:
        raise ValueError(
            f"No fallback Martini tail rule defined for {length}:{double_bonds}. "
            f"Please check the ;@BEADS line in the .itp file."
        )

    d_positions = double_bond_positions[key]

    beads = []
    for i in range(1, n_beads + 1):
        if i in d_positions:
            beads.append(f"D{i}{suffix}")
        else:
            beads.append(f"C{i}{suffix}")

    return beads


def bead_smiles(bead):
    bead_smiles_map = {
        "NC3": "N(C)(C)(C)CC",
        "NH3": "N(C)(C)",
        "PO4": "OP(=O)(O)(O)",
        "GL1": "COC(=O)",
        "GL2": "COC(=O)",
        "CNO": "OC(=O)CNC",
        "GL0": "OCCOCO",
        "C1": "C1(O)C",
        "C2": "CO",
        "C3": "CO",
        "C4": "CO",
        "INO": "C1(O)C(O)C(O)C(O)C(O)O",
        "AM2": "NC(=O)",
        "OH1": "CO",
        "COH": "CO",
        "GM1":  "C1(O)C(O)C(O)C(O)O",
        "GM2":  "C1(O)C(O)C(O)C(O)O",
        "GM3":  "C1(O)C(O)C(O)C(O)O",
        "GM4":  "C1(O)C(O)C(O)C(O)O",
        "GM5":  "C1(O)C(O)C(O)C(O)O",

        "GM6":  "NC(CO)C(O)O",
        "GM7":  "NC(CO)C(O)O",
        "GM8":  "NC(CO)C(O)O",
        "GM9":  "NC(CO)C(O)O",
        "GM10": "NC(CO)C(O)O",

        "GM11": "OC(=O)C(O)C(O)",
        "GM12": "OC(=O)C(O)C(O)",
        "GM13": "OC(=O)C(O)C(O)",

        "GM14": "NC(CO)C(O)C(=O)O",
        "GM15": "NC(CO)C(O)C(=O)O",
        "GM16": "NC(CO)C(O)C(=O)O",
        "OH1" : "NC",
        "AM1": "NC(=O)",
        "AM2": "C(=O)",
        "T1A": "CCC",
        "ROH": "CCC(O)",
    }

    if bead in bead_smiles_map:
        return bead_smiles_map[bead]

    if bead.startswith("D"):
        return "CC(C=C)"

    return "CCCC"

def normalize_bead_name(bead):
    """
    Normalize bead names from Martini .itp files to the LiMBS convention.
    Example: some PC lipids may use NC4 in .itp, but LiMBS uses NC3.
    """
    if bead == "NC4":
        return "NC3"
    return bead


def generate_chol_block(name):
    return f"""{name}=CG|[
[CCC(O):^ROH] [#TailA] [#TailB].
{{#TailA=[>][(C=C)C:^R1][CCC:^R2][CCC:^R3][CC:^R4][CC:^R5][CCC:^R6]}}.
{{#TailB=[>][CCC:^C1][CCCCC:^C2]}}
|
[^ROH] [#TailA] [#TailB].
{{#TailA=[>][^R1][^R2][^R3][^R4][^R5][^R6]}}.
{{#TailB=[>][^C1][^C2]}}
]"""


def generate_explicit_block(name, bead_list):
    """
    Generate a grouped LiMBS CG block using exact Martini bead names from .itp.
    Head beads stay on the main line.
    Tail beads ending in A go into #TailA.
    Tail beads ending in B go into #TailB.
    """
    bead_list = [normalize_bead_name(b) for b in bead_list]

    tail_a = [b for b in bead_list if b.endswith("A")]
    tail_b = [b for b in bead_list if b.endswith("B")]
    head = [b for b in bead_list if b not in tail_a and b not in tail_b]

    head_atoms = " ".join(f"[{bead_smiles(b)}:^{b}]" for b in head)
    head_refs = " ".join(f"[^{b}]" for b in head)

    tail_a_atoms = "".join(f"[{bead_smiles(b)}:^{b}]" for b in tail_a)
    tail_b_atoms = "".join(f"[{bead_smiles(b)}:^{b}]" for b in tail_b)

    tail_a_refs = "".join(f"[^{b}]" for b in tail_a)
    tail_b_refs = "".join(f"[^{b}]" for b in tail_b)

    return f"""{name}=CG|[
{head_atoms} [#TailA] [#TailB].
{{#TailA=[>]{tail_a_atoms}}}.
{{#TailB=[>]{tail_b_atoms}}}
|
{head_refs} [#TailA] [#TailB].
{{#TailA=[>]{tail_a_refs}}}.
{{#TailB=[>]{tail_b_refs}}}
]"""


def normalize_bead_name(bead):
    if bead == "NC4":
        return "NC3"
    return bead
def generate_cg_block(name, lipid_class, tail_a, tail_b):
    if lipid_class in ("STEROL", "CHOL"):
        return generate_chol_block(name)

    itp_beads = read_itp_beads_for_lipid(name)
    if itp_beads is not None:
        return generate_explicit_block(name, itp_beads)

    if lipid_class in ("GM1", "GM3"):
        bead_list = tail_a + tail_b
        return generate_explicit_block(name, bead_list)

    if lipid_class not in HEADGROUPS:
        raise ValueError(f"Unsupported lipid class: {lipid_class}")

    head = HEADGROUPS[lipid_class] 
def generate_cg_block(name, lipid_class, tail_a, tail_b):
    lipid_class = lipid_class.upper()

    if lipid_class in ("STEROL", "CHOL"):
        return generate_chol_block(name)

    # First try to use exact Martini beads from .itp
    itp_beads = read_itp_beads_for_lipid(name)
    if itp_beads is not None:
        return generate_explicit_block(name, itp_beads)

    # Glycolipids such as PNG1, PNG3, DPG1, DPG3, XNG1, XNG3
    if lipid_class in ("GM1", "GM3"):
        bead_list = tail_a + tail_b
        return generate_explicit_block(name, bead_list)

    if lipid_class not in HEADGROUPS:
        raise ValueError(f"Unsupported lipid class: {lipid_class}")

    head = HEADGROUPS[lipid_class]

    tail_a_beads = make_tail(tail_a[0], tail_a[1], "A")
    tail_b_beads = make_tail(tail_b[0], tail_b[1], "B")

    tail_a_atoms = "".join(f"[{bead_smiles(b)}:^{b}]" for b in tail_a_beads)
    tail_b_atoms = "".join(f"[{bead_smiles(b)}:^{b}]" for b in tail_b_beads)

    tail_a_refs = "".join(f"[^{b}]" for b in tail_a_beads)
    tail_b_refs = "".join(f"[^{b}]" for b in tail_b_beads)

    return f"""{name}=CG|[
{head["atoms"]} [#TailA] [#TailB].
{{#TailA=[>]{tail_a_atoms}}}.
{{#TailB=[>]{tail_b_atoms}}}
|
{head["refs"]} [#TailA] [#TailB].
{{#TailA=[>]{tail_a_refs}}}.
{{#TailB=[>]{tail_b_refs}}}
]"""
