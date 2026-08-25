"""
cg_templates.py

Reusable coarse-grained structural template generation for LiMBS1.0.

This module generates LiMBS CG species blocks. It does not generate
coordinates, force-field parameters, or simulation-ready membrane systems.

Important conventions
1. When an annotated Martini .itp file is available, the observed Martini
   TailA/TailB bead patterns are used to determine which supplied tail
   composition belongs to TailA and which belongs to TailB.

2. Verified LiMBS fragment overrides are used for special tails that cannot
   be represented correctly by the generic bead-fragment rules.

"""

from pathlib import Path
import os
import re


# Optional Martini 3 topology directory

ITP_DIR = Path(
    os.environ.get(
        "LIMBS_MARTINI3_ITP_DIR",
        "~/MARTINI3/ff",
    )
).expanduser()


# Headgroup templates used by the fallback generator

HEADGROUPS = {
    "PC": {
        "atoms": (
            "[N(C)(C)(C)CC:^NC3] "
            "[OP(=O)(O)(O):^PO4] "
            "[CCOC(=O):^GL1] "
            "[COC(=O):^GL2]"
        ),
        "refs": "[^NC3] [^PO4] [^GL1] [^GL2]",
    },

    "PE": {
        "atoms": (
            "[NCC:^NH3] "
            "[OP(=O)(O)(O):^PO4] "
            "[CCOC(=O):^GL1] "
            "[COC(=O):^GL2]"
        ),
        "refs": "[^NH3] [^PO4] [^GL1] [^GL2]",
    },

    "PS": {
        "atoms": (
            "[OC(=O)CNC:^CNO] "
            "[OP(=O)(O)(O):^PO4] "
            "[CCOC(=O):^GL1] "
            "[COC(=O):^GL2]"
        ),
        "refs": "[^CNO] [^PO4] [^GL1] [^GL2]",
    },

    "PG": {
        "atoms": (
            "[OCCOCO:^GL0] "
            "[OP(=O)(O)(O):^PO4] "
            "[CCOC(=O):^GL1] "
            "[COC(=O):^GL2]"
        ),
        "refs": "[^GL0] [^PO4] [^GL1] [^GL2]",
    },

    "PA": {
        "atoms": (
            "[OP(=O)(O)(O):^PO4] "
            "[COC(=O):^GL1] "
            "[CCOC(=O):^GL2]"
        ),
        "refs": "[^PO4] [^GL1] [^GL2]",
    },

    "PI": {
        "atoms": (
            "[CCOCO:^C1][CO:^C2][OC:^C3][OC:^C4] "
            "[OP(=O)(O)(O):^PO4] "
            "[CCOC(=O):^GL1] "
            "[COC(=O):^GL2]"
        ),
        "refs": "[^C1][^C2][^C3][^C4] [^PO4] [^GL1] [^GL2]",
    },

    "DG": {
        "atoms": "[CO:^COH] [CCOC(=O):^GL1] [COC(=O):^GL2]",
        "refs": "[^COH] [^GL1] [^GL2]",
    },

    "SM": {
        "atoms": (
            "[N(C)(C)(C)CC:^NC3] "
            "[OP(=O)(O)(O):^PO4] "
            "[NC:^OH1] "
            "[C(=O):^AM2]"
        ),
        "refs": "[^NC3] [^PO4] [^OH1] [^AM2]",
    },

    "CE": {
        "atoms": "[CO:^COH] [NC:^AM2] [CO:^OH1]",
        "refs": "[^COH] [^AM2] [^OH1]",
    },
}


# Head  bead-to-fragment definitions

BEAD_SMILES = {
    "NC3": "N(C)(C)(C)CC",
    "NH3": "NCC",
    "PO4": "OP(=O)(O)(O)",
    "GL1": "CCOC(=O)",
    "GL2": "COC(=O)",
    "CNO": "OC(=O)CNC",
    "GL0": "OCCOCO",
    "C1": "C1(O)C",
    "C2": "CO",
    "C3": "CO",
    "C4": "CO",
    "INO": "C1(O)C(O)C(O)C(O)C(O)O",
    "COH": "CO",
    "AM1": "NC(=O)",

    # Sphingomyelin sphingosine delta-4 double bond.
    "T1A": "(C=C)C",

    "ROH": "CCC(O)",
}


CLASS_SPECIFIC_BEAD_SMILES = {
    "SM": {
        "OH1": "NC",
        "AM2": "C(=O)",
    },

    "CE": {
        "COH": "CO",
        "AM2": "NC",
        "OH1": "CO",
    },
}


# Verified tail-specific LiMBS fragment definitions

TAIL_FRAGMENT_OVERRIDES = {
    # 20:3, e.g. PQPE
    (20, 3): {
        "C1": "CCC",
    },

    # 20:4, e.g. PAPE / DAPE / PAPA
    (20, 4): {
        "C1": "CCC",
        "D2": "(C=C)C",
        "D3": "(C=C)C",
        "D4": "(C=C)C(C=C)",
        "C5": "CCCCC",
    },

    # 20:5
    # D1 is one carbon shorter than D2-D5.
    (20, 5): {
        "D1": "C(C=C)",
        "D2": "C(C=C)C",
        "D3": "C(C=C)C",
        "D4": "C(C=C)C",
        "D5": "C(C=C)C",
    },

    # 22:1, e.g. PEPC
    # Only C1 needs the five-carbon representation.
    (22, 1): {
        "C1": "CCCCC",
    },

    # 22:6, e.g. PUPC / PUPE
    (22, 6): {
        "D1": "C(C=C)C",
        "D2": "C(C=C)C",
        "D3": "C(C=C)C",
        "D4": "C(C=C)C",
        "D5": "C(C=C)CC",
    },
}


# Lipid-class-specific tail overrides

CLASS_TAIL_FRAGMENT_OVERRIDES = {
    # DPSM and other SM 18:0 acyl chains:
    # C1 must contain five carbons.
    ("SM", (18, 0)): {
        "C1": "CCCCC",
    },
}


# Verified special tail bead patterns

EXPLICIT_TAIL_BEAD_PATTERNS = {
    # 20:3
    (20, 3): [
        "C1",
        "D2",
        "D3",
        "D4",
        "C5",
    ],

    # 20:4
    (20, 4): [
        "C1",
        "D2",
        "D3",
        "D4",
        "C5",
    ],

    # 20:5
    (20, 5): [
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
    ],

    # 22:6
    (22, 6): [
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
    ],
}


# Generic fallback bead counts and D-bead positions

TAIL_BEAD_COUNTS = {
    12: 3,
    16: 4,
    18: 4,
    20: 5,
    22: 6,
    24: 6,
}


DOUBLE_BOND_POSITIONS = {
    # 12 carbon
    (12, 0): [],

    # 16 carbon
    (16, 0): [],
    (16, 1): [2],

    # 18 carbon
    (18, 0): [],
    (18, 1): [2],
    (18, 2): [2, 3],

    # 20 carbon
    (20, 0): [],
    (20, 1): [3],
    (20, 2): [3, 4],
    # 20:3, 20:4, 20:5 are explicit above.

    # 22 carbon
    (22, 0): [],
    (22, 1): [4],
    (22, 2): [4, 5],
    (22, 4): [3, 4, 5, 6],
    (22, 5): [2, 3, 4, 5, 6],
    # 22:6 is explicit above.

    # 24 carbon
    (24, 0): [],
    (24, 1): [4],
}


# Annotated Martini .itp support

def read_itp_beads_for_lipid(lipid_name):
    """
    Read bead ordering from an annotated Martini .itp file.

    Returns a list of bead names when a compatible ;@INSANE / ;@BEADS
    annotation is found. Otherwise returns None.
    """
    lipid_name = lipid_name.upper()

    if not ITP_DIR.exists():
        return None

    for itp in ITP_DIR.glob("*.itp"):
        text = itp.read_text(errors="ignore")

        insane_matches = list(
            re.finditer(
                r";@INSANE[^\n]*",
                text,
            )
        )

        for index, match in enumerate(insane_matches):
            line = match.group(0)

            alname_match = re.search(
                r"alname=([A-Za-z0-9_.-]+)",
                line,
            )

            if not alname_match:
                continue

            if alname_match.group(1).upper() != lipid_name:
                continue

            block_end = (
                insane_matches[index + 1].start()
                if index + 1 < len(insane_matches)
                else len(text)
            )

            chunk = text[
                match.start():block_end
            ]

            beads_match = re.search(
                r";@BEADS\s+(.+)",
                chunk,
            )

            if beads_match:
                return beads_match.group(1).split()

    return None


# Tail helpers

def normalize_tail_spec(tail_spec):
    """
    Normalize a tail specification to:

        (carbon_count, double_bond_count)

    Example:
        [18, 0] -> (18, 0)
    """
    if tail_spec is None:
        return None

    if len(tail_spec) != 2:
        raise ValueError(
            f"Invalid tail specification {tail_spec!r}; "
            "expected (carbon_count, double_bond_count)."
        )

    return (
        int(tail_spec[0]),
        int(tail_spec[1]),
    )


def _is_standard_tail_bead(bead):
    """Return True for C1A, D2A, C3B, etc."""
    return re.fullmatch(
        r"[CD]\d+[AB]",
        bead,
    ) is not None


def _is_tail_bead(bead, lipid_class=None):
    """
    Return True for a LiMBS/Martini tail bead.

    T1A is treated as part of TailA only for sphingomyelin.
    """
    if _is_standard_tail_bead(bead):
        return True

    if (
        lipid_class is not None
        and lipid_class.upper() == "SM"
        and bead == "T1A"
    ):
        return True

    return False


def _base_tail_labels(beads):
    """
    Convert standard bead labels to unsuffixed labels.

    T1A remains T1 because it is a special sphingomyelin bead.
    """
    result = []

    for bead in beads:
        if bead == "T1A":
            result.append("T1")
        elif _is_standard_tail_bead(bead):
            result.append(bead[:-1])

    return result


def _tail_unsaturation_count(beads):
    """
    Approximate the number of unsaturated Martini tail beads.

    T1A counts as one unsaturated bead because it contains the sphingosine delta-4 double bond.
    """
    count = 0

    for bead in beads:
        if bead == "T1A":
            count += 1
        elif re.fullmatch(r"D\d+[AB]", bead):
            count += 1

    return count


def tail_contains_d2(tail_beads):
    """Return True when the specified tail contains D2A or D2B."""
    return any(
        re.fullmatch(
            r"D2[AB]",
            bead,
        )
        for bead in tail_beads
    )


# Fallback tail bead-name generation

def make_tail(length, double_bonds, suffix):
    """
    Generate fallback Martini-style bead labels for a standard lipid tail.

    suffix must be "A" or "B".
    """
    if suffix not in ("A", "B"):
        raise ValueError(
            f"Tail suffix must be 'A' or 'B', received {suffix!r}."
        )

    length = int(length)
    double_bonds = int(double_bonds)

    key = (
        length,
        double_bonds,
    )

    # Use verified special patterns first.
    if key in EXPLICIT_TAIL_BEAD_PATTERNS:
        return [
            f"{bead}{suffix}"
            for bead in EXPLICIT_TAIL_BEAD_PATTERNS[key]
        ]

    if length not in TAIL_BEAD_COUNTS:
        raise ValueError(
            f"Unsupported tail length: {length}"
        )

    if key not in DOUBLE_BOND_POSITIONS:
        raise ValueError(
            f"No fallback Martini tail rule defined for "
            f"{length}:{double_bonds}. "
            "Use a compatible annotated .itp file or add a verified "
            "fallback rule."
        )

    d_positions = set(
        DOUBLE_BOND_POSITIONS[key]
    )

    return [
        (
            f"{'D' if index in d_positions else 'C'}"
            f"{index}{suffix}"
        )
        for index in range(
            1,
            TAIL_BEAD_COUNTS[length] + 1,
        )
    ]


def make_sm_tail_a(tail_spec):
    """
    Generate fallback TailA for a sphingomyelin sphingoid chain.

    For a monounsaturated SM sphingoid chain, the first tail bead is T1A.
    Example for 18:1:

        T1A C2A C3A C4A

    Example for 24:1:

        T1A C2A C3A C4A C5A C6A
    """
    tail_spec = normalize_tail_spec(
        tail_spec
    )

    length, double_bonds = tail_spec

    if (
        double_bonds == 1
        and length in TAIL_BEAD_COUNTS
    ):
        bead_count = TAIL_BEAD_COUNTS[
            length
        ]

        return (
            ["T1A"]
            + [
                f"C{index}A"
                for index in range(
                    2,
                    bead_count + 1,
                )
            ]
        )

    # For an SM tail not covered by the special sphingoid rule,
    # use the standard fallback rather than guessing.
    return make_tail(
        length,
        double_bonds,
        "A",
    )


# Match supplied tail metadata to annotated Martini TailA / TailB

def _expected_tail_beads_for_scoring(
    tail_spec,
):
    """
    Return a generic expected bead pattern for scoring an annotated tail.
    """
    tail_spec = normalize_tail_spec(
        tail_spec
    )

    return make_tail(
        tail_spec[0],
        tail_spec[1],
        "A",
    )


def _tail_pattern_score(
    observed_beads,
    tail_spec,
):
    """
    Score how well an observed Martini tail matches a supplied tail spec.

    Lower is better.

    The score considers:
    - number of beads
    - number of unsaturated beads
    - exact C/D bead labels when possible

    T1A is handled specially because it is the sphingomyelin Delta-4
    unsaturated bead and does not follow the ordinary C1/D2 naming.
    """
    tail_spec = normalize_tail_spec(
        tail_spec
    )

    expected_beads = (
        _expected_tail_beads_for_scoring(
            tail_spec
        )
    )

    observed_count = len(
        observed_beads
    )

    expected_count = len(
        expected_beads
    )

    observed_unsat = (
        _tail_unsaturation_count(
            observed_beads
        )
    )

    expected_unsat = (
        _tail_unsaturation_count(
            expected_beads
        )
    )

    # Bead-count disagreement is the strongest penalty.
    score = (
        20
        * abs(
            observed_count
            - expected_count
        )
    )

    # Unsaturation disagreement is also important.
    score += (
        6
        * abs(
            observed_unsat
            - expected_unsat
        )
    )

    # If the observed tail contains T1A, do not force exact C/D positional
    # matching because T1A is a special sphingoid representation.
    if "T1A" in observed_beads:
        return score

    observed_base = (
        _base_tail_labels(
            observed_beads
        )
    )

    expected_base = (
        _base_tail_labels(
            expected_beads
        )
    )

    score += sum(
        observed_bead != expected_bead
        for observed_bead, expected_bead
        in zip(
            observed_base,
            expected_base,
        )
    )

    return score


def resolve_itp_tail_specs(
    tail_a_beads,
    tail_b_beads,
    first_tail_spec,
    second_tail_spec,
):
    """
    Match two supplied tail compositions to annotated Martini TailA/TailB.

    Both possible assignments are scored. If the scores tie, the supplied
    order is retained.
    """
    first_tail_spec = normalize_tail_spec(
        first_tail_spec
    )

    second_tail_spec = normalize_tail_spec(
        second_tail_spec
    )

    direct_score = (
        _tail_pattern_score(
            tail_a_beads,
            first_tail_spec,
        )
        + _tail_pattern_score(
            tail_b_beads,
            second_tail_spec,
        )
    )

    swapped_score = (
        _tail_pattern_score(
            tail_a_beads,
            second_tail_spec,
        )
        + _tail_pattern_score(
            tail_b_beads,
            first_tail_spec,
        )
    )

    if swapped_score < direct_score:
        return (
            second_tail_spec,
            first_tail_spec,
        )

    return (
        first_tail_spec,
        second_tail_spec,
    )


def apply_explicit_tail_pattern_if_needed(
    tail_beads,
    tail_spec,
):
    """
    Replace an annotated standard tail pattern only when a verified LiMBS
    pattern is defined for that exact tail composition.

    SM T1A-containing tails are left unchanged.
    """
    tail_spec = normalize_tail_spec(
        tail_spec
    )

    if "T1A" in tail_beads:
        return tail_beads

    if tail_spec not in EXPLICIT_TAIL_BEAD_PATTERNS:
        return tail_beads

    if not tail_beads:
        return tail_beads

    suffixes = {
        bead[-1]
        for bead in tail_beads
        if _is_standard_tail_bead(
            bead
        )
    }

    if len(suffixes) != 1:
        raise ValueError(
            "Could not determine a unique A/B suffix for "
            f"tail beads: {tail_beads}"
        )

    suffix = next(
        iter(suffixes)
    )

    return [
        f"{bead}{suffix}"
        for bead in EXPLICIT_TAIL_BEAD_PATTERNS[
            tail_spec
        ]
    ]


# Fragment generation

def bead_smiles(
    bead,
    lipid_class=None,
):
    """
    Return the fragment representation for a non-tail bead.
    """
    lipid_class = (
        lipid_class.upper()
        if lipid_class
        else None
    )

    if (
        lipid_class
        in CLASS_SPECIFIC_BEAD_SMILES
    ):
        class_map = (
            CLASS_SPECIFIC_BEAD_SMILES[
                lipid_class
            ]
        )

        if bead in class_map:
            return class_map[
                bead
            ]

    if bead in BEAD_SMILES:
        return BEAD_SMILES[
            bead
        ]

    # Legacy direct-call fallback for standard tail labels.
    if re.fullmatch(
        r"C\d+[AB]",
        bead,
    ):
        return "CCCC"

    if re.fullmatch(
        r"D\d+[AB]",
        bead,
    ):
        return "CC(C=C)"

    raise ValueError(
        f"No LiMBS fragment mapping defined for bead {bead!r}"
        + (
            f" in lipid class {lipid_class!r}."
            if lipid_class
            else "."
        )
    )


def tail_bead_smiles(
    bead,
    tail_beads,
    tail_spec=None,
    lipid_class=None,
):
    """
    Return the LiMBS structural fragment for one tail bead.

    Priority:
        1. SM T1A special fragment
        2. lipid-class-specific tail override
        3. general tail-specific override
        4. generic fallback
    """

    # Sphingomyelin sphingosine Delta-4 double bond
    if bead == "T1A":
        return "(C=C)C"

    match = re.fullmatch(
        r"([CD])(\d+)([AB])",
        bead,
    )

    if not match:
        raise ValueError(
            f"Invalid tail bead label: {bead!r}"
        )

    bead_type = match.group(1)
    bead_index = int(
        match.group(2)
    )

    base_bead = (
        f"{bead_type}{bead_index}"
    )

    tail_spec = normalize_tail_spec(
        tail_spec
    )

    lipid_class = (
        lipid_class.upper()
        if lipid_class
        else None
    )

    # Lipid-class-specific override
    if (
        lipid_class is not None
        and tail_spec is not None
    ):
        class_key = (
            lipid_class,
            tail_spec,
        )

        class_override = (
            CLASS_TAIL_FRAGMENT_OVERRIDES.get(
                class_key
            )
        )

        if (
            class_override is not None
            and base_bead
            in class_override
        ):
            return class_override[
                base_bead
            ]

    # General tail-specific override
    if tail_spec is not None:
        override = (
            TAIL_FRAGMENT_OVERRIDES.get(
                tail_spec
            )
        )

        if (
            override is not None
            and base_bead
            in override
        ):
            return override[
                base_bead
            ]

    # Generic fallback
    if bead_type == "D":
        return "CC(C=C)"

    if (
        bead_type == "C"
        and bead_index == 1
    ):
        if tail_contains_d2(
            tail_beads
        ):
            return "CCCCC"

        return "CCC"

    return "CCCC"


# Display-label normalization

def normalize_bead_name(
    bead,
):
    """
    Apply LiMBS compatibility aliases without modifying source topology.
    """
    if bead == "NC4":
        return "NC3"

    return bead


# Formatting helpers

def _join_atoms(
    beads,
    lipid_class,
):
    return " ".join(
        (
            f"["
            f"{bead_smiles(bead, lipid_class)}"
            f":^{bead}]"
        )
        for bead in beads
    )


def _join_refs(
    beads,
):
    return " ".join(
        f"[^{bead}]"
        for bead in beads
    )


def _join_tail_atoms(
    beads,
    tail_spec=None,
    lipid_class=None,
):
    return "".join(
        (
            f"["
            f"{tail_bead_smiles(bead, beads, tail_spec, lipid_class)}"
            f":^{bead}]"
        )
        for bead in beads
    )


def _join_tail_refs(
    beads,
):
    return "".join(
        f"[^{bead}]"
        for bead in beads
    )


# Cholesterol

def generate_chol_block(
    name,
):
    """
    Generate the LiMBS CG block used for cholesterol.
    """
    return f"""{name}=CG|[
[CCC(O):^ROH] [#TailA] [#TailB]
{{#TailA=[>][(C=C)C:^R1][CCC:^R2][CCC:^R3][CC:^R4][CC:^R5][CCC:^R6]}}
{{#TailB=[>][CCC:^C1][CCCCC:^C2]}}
|
[^ROH] [#TailA] [#TailB]
{{#TailA=[>][^R1][^R2][^R3][^R4][^R5][^R6]}}
{{#TailB=[>][^C1][^C2]}}
]"""


# Explicit block generation from annotated .itp bead ordering

def generate_explicit_block(
    name,
    bead_list,
    lipid_class,
    tail_a_spec=None,
    tail_b_spec=None,
):
    """
    Generate a LiMBS CG block from annotated Martini .itp bead ordering.

    The two supplied tail specs are automatically matched to observed TailA
    and TailB before fragment generation.

    For SM:
        T1A is included in TailA rather than in the headgroup.
    """
    lipid_class = (
        lipid_class.upper()
    )

    bead_list = [
        normalize_bead_name(
            bead
        )
        for bead in bead_list
    ]

    # Observed TailA
    observed_tail_a = [
        bead
        for bead in bead_list
        if (
            re.fullmatch(
                r"[CD]\d+A",
                bead,
            )
            or (
                lipid_class == "SM"
                and bead == "T1A"
            )
        )
    ]

    # Observed TailB
    observed_tail_b = [
        bead
        for bead in bead_list
        if re.fullmatch(
            r"[CD]\d+B",
            bead,
        )
    ]

    # Everything else belongs to the non-tail part
    head_beads = [
        bead
        for bead in bead_list
        if (
            bead not in observed_tail_a
            and bead not in observed_tail_b
        )
    ]

    # Resolve supplied tail metadata against actual .itp A/B bead patterns
    (
        resolved_tail_a_spec,
        resolved_tail_b_spec,
    ) = resolve_itp_tail_specs(
        observed_tail_a,
        observed_tail_b,
        tail_a_spec,
        tail_b_spec,
    )

    # Apply verified special bead patterns
    tail_a_beads = (
        apply_explicit_tail_pattern_if_needed(
            observed_tail_a,
            resolved_tail_a_spec,
        )
    )

    tail_b_beads = (
        apply_explicit_tail_pattern_if_needed(
            observed_tail_b,
            resolved_tail_b_spec,
        )
    )
    #Build head
    head_atoms = _join_atoms(
        head_beads,
        lipid_class,
    )

    head_refs = _join_refs(
        head_beads
    )

    tail_a_atoms = _join_tail_atoms(
        tail_a_beads,
        resolved_tail_a_spec,
        lipid_class,
    )

    tail_b_atoms = _join_tail_atoms(
        tail_b_beads,
        resolved_tail_b_spec,
        lipid_class,
    )

    tail_a_refs = _join_tail_refs(
        tail_a_beads
    )

    tail_b_refs = _join_tail_refs(
        tail_b_beads
    )

    return f"""{name}=CG|[
{head_atoms} [#TailA] [#TailB]
{{#TailA=[>]{tail_a_atoms}}}
{{#TailB=[>]{tail_b_atoms}}}
|
{head_refs} [#TailA] [#TailB]
{{#TailA=[>]{tail_a_refs}}}
{{#TailB=[>]{tail_b_refs}}}
]"""


# Main CG generator

def generate_cg_block(
    name,
    lipid_class,
    tail_a,
    tail_b,
):
    """
    Generate a LiMBS coarse-grained species block.

    Existing website interface:

        generate_cg_block(
            name,
            lipid_class,
            tail_a,
            tail_b,
        )
    """
    lipid_class = (
        lipid_class.upper()
    )

    # Cholesterol / sterol
    if lipid_class in (
        "STEROL",
        "CHOL",
    ):
        return generate_chol_block(
            name
        )

    # Validate class
    if lipid_class not in HEADGROUPS:
        raise ValueError(
            f"Unsupported lipid class: "
            f"{lipid_class}"
        )

    if (
        tail_a is None
        or tail_b is None
    ):
        raise ValueError(
            f"{name}: tail_a and tail_b are required for "
            "non-sterol CG template generation."
        )

    tail_a = normalize_tail_spec(
        tail_a
    )

    tail_b = normalize_tail_spec(
        tail_b
    )

    # Preferred path: annotated Martini .itp bead ordering
    itp_beads = (
        read_itp_beads_for_lipid(
            name
        )
    )

    if itp_beads is not None:
        return generate_explicit_block(
            name=name,
            bead_list=itp_beads,
            lipid_class=lipid_class,
            tail_a_spec=tail_a,
            tail_b_spec=tail_b,
        )

    # Fallback path
    head = HEADGROUPS[
        lipid_class
    ]

    if lipid_class == "SM":
        # For SM the first supplied tail is the sphingoid chain.
        tail_a_beads = make_sm_tail_a(
            tail_a
        )
    else:
        tail_a_beads = make_tail(
            tail_a[0],
            tail_a[1],
            "A",
        )

    tail_b_beads = make_tail(
        tail_b[0],
        tail_b[1],
        "B",
    )

    tail_a_atoms = _join_tail_atoms(
        tail_a_beads,
        tail_a,
        lipid_class,
    )

    tail_b_atoms = _join_tail_atoms(
        tail_b_beads,
        tail_b,
        lipid_class,
    )

    tail_a_refs = _join_tail_refs(
        tail_a_beads
    )

    tail_b_refs = _join_tail_refs(
        tail_b_beads
    )

    return f"""{name}=CG|[
{head["atoms"]} [#TailA] [#TailB]
{{#TailA=[>]{tail_a_atoms}}}
{{#TailB=[>]{tail_b_atoms}}}
|
{head["refs"]} [#TailA] [#TailB]
{{#TailA=[>]{tail_a_refs}}}
{{#TailB=[>]{tail_b_refs}}}
]"""


# Internal self-tests

def _run_self_tests():
    """
    Run small internal tests that do not depend on the user's Martini .itp
    directory.
    """

    # 12:0 support
    assert make_tail(
        12,
        0,
        "A",
    ) == [
        "C1A",
        "C2A",
        "C3A",
    ]

    # PQPE 20:3 first bead
    assert tail_bead_smiles(
        "C1A",
        [
            "C1A",
            "D2A",
            "D3A",
            "D4A",
            "C5A",
        ],
        (20, 3),
        "PE",
    ) == "CCC"

    # PAPE/DAPE 20:4 fragments
    assert tail_bead_smiles(
        "D4A",
        [
            "C1A",
            "D2A",
            "D3A",
            "D4A",
            "C5A",
        ],
        (20, 4),
        "PE",
    ) == "(C=C)C(C=C)"

    # 20:5 first bead
    assert tail_bead_smiles(
        "D1A",
        [
            "D1A",
            "D2A",
            "D3A",
            "D4A",
            "D5A",
        ],
        (20, 5),
        "PE",
    ) == "C(C=C)"

    # PEPC 22:1 C1
    assert tail_bead_smiles(
        "C1A",
        [
            "C1A",
            "C2A",
            "C3A",
            "D4A",
            "C5A",
            "C6A",
        ],
        (22, 1),
        "PC",
    ) == "CCCCC"

    # SM sphingosine T1A
    assert tail_bead_smiles(
        "T1A",
        [
            "T1A",
            "C2A",
            "C3A",
            "C4A",
        ],
        (18, 1),
        "SM",
    ) == "(C=C)C"

    # DPSM 18:0 acyl C1B
    assert tail_bead_smiles(
        "C1B",
        [
            "C1B",
            "C2B",
            "C3B",
            "C4B",
        ],
        (18, 0),
        "SM",
    ) == "CCCCC"

    # Synthetic PSM-like annotated topology.
    psm_test = generate_explicit_block(
        name="PSM_TEST",
        bead_list=[
            "NC3",
            "PO4",
            "OH1",
            "AM2",
            "T1A",
            "C2A",
            "C3A",
            "C4A",
            "C1B",
            "C2B",
            "C3B",
            "C4B",
        ],
        lipid_class="SM",
        tail_a_spec=(18, 1),
        tail_b_spec=(16, 0),
    )

    assert "[(C=C)C:^T1A]" in psm_test
    assert "[^T1A]" in psm_test

    # Synthetic DPSM-like annotated topology.
    dpsm_test = generate_explicit_block(
        name="DPSM_TEST",
        bead_list=[
            "NC3",
            "PO4",
            "OH1",
            "AM2",
            "T1A",
            "C2A",
            "C3A",
            "C4A",
            "C1B",
            "C2B",
            "C3B",
            "C4B",
        ],
        lipid_class="SM",
        tail_a_spec=(18, 1),
        tail_b_spec=(18, 0),
    )

    assert "[CCCCC:^C1B]" in dpsm_test

    print("cg_templates.py self-tests passed.")


if __name__ == "__main__":
    _run_self_tests()
