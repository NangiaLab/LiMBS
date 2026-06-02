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
        "atoms": "[C1(O)C(O)C(O)C(O)C(O)O:^INO] [OP(=O)(O)(O):^PO4] [COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^INO] [^PO4] [^GL1] [^GL2]",
    },
    "DG": {
        "atoms": "[COC(=O):^GL1] [COC(=O):^GL2]",
        "refs":  "[^GL1] [^GL2]",
    },
    "SM": {
        "atoms": "[N(C)(C)(C)CC:^NC3] [OP(=O)(O)(O):^PO4] [NC(=O):^AM2]",
        "refs":  "[^NC3] [^PO4] [^AM2]",
    },
}


def make_tail(length, double_bonds, suffix):
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
    beads = []

    for i in range(1, n_beads + 1):
        if i >= 2 and i < 2 + double_bonds:
            beads.append(f"D{i}{suffix}")
        else:
            beads.append(f"C{i}{suffix}")

    return beads


def bead_smiles(bead):
    if bead.startswith("D"):
        return "CC(C=C)"
    return "CCCC"


def generate_chol_block(name):
    return f"""{name}=CG|[
[CCC(O):^ROH] [(C=C)C:^R1] [CCC:^R2] [CCC:^R3] [CC:^R4] [CC:^R5] [CCC:^R6] [CCC:^C1] [CCCCC:^C2]
|
[^ROH] [^R1] [^R2] [^R3] [^R4] [^R5] [^R6] [^C1] [^C2]
]"""


def generate_cg_block(name, lipid_class, tail_a, tail_b):
    if lipid_class in ("STEROL", "CHOL"):
        return generate_chol_block(name)

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
