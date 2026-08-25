"""
aa_templates.py

Reusable atomistic (AA) species templates for the LiMBS1.0 website.

This module generates LiMBS AA species blocks only. It does not generate
coordinates, topology/parameter files, or simulation-ready AA membranes.
"""

from typing import Dict, List, Any


AA_LIPID_TEMPLATES: Dict[str, str] = {
    "CHOL": (
        "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2"
        "[C@@]1(CC[C@H]3[C@H]2CC=C4"
        "[C@@]3(CC[C@@H](C4)O)C)C"
    ),
    "DOPE": (
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "[N]CCOP(=O)(O)OCCO(=O)[#T18_1]CO(=O)[#T18_1]"
    ),
    "DOPG": (
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "OCCOCOP(=O)(O)OCCO(=O)[#T18_1]CO(=O)[#T18_1]"
    ),
    "DOPS": (
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "[N]C(C(=O)O)COP(=O)(O)OCCO(=O)[#T18_1]CO(=O)[#T18_1]"
    ),
    "DPPC": (
        "{#T16=CCCCCCCCCCCCCCCC}\n\n"
        "[N](C)(C)(C)CCOP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T16]"
    ),
    "POPA": (
        "{#T16=CCCCCCCCCCCCCCCC}\n"
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "OP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T18_1]"
    ),
    "POPC": (
        "{#T16=CCCCCCCCCCCCCCCC}\n"
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "[N](C)(C)(C)CCOP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T18_1]"
    ),
    "POPE": (
        "{#T16=CCCCCCCCCCCCCCCC}\n"
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "[N]CCOP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T18_1]"
    ),
    "POPG": (
        "{#T16=CCCCCCCCCCCCCCCC}\n"
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "OCCOCOP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T18_1]"
    ),
    "POPS": (
        "{#T16=CCCCCCCCCCCCCCCC}\n"
        "{#T18_1=CCCCCCCCC=CCCCCCCCC}\n\n"
        "[N]C(C(=O)O)COP(=O)(O)OCCO(=O)[#T16]CO(=O)[#T18_1]"
    ),
}


AA_ENVIRONMENT_TEMPLATES: Dict[str, str] = {
    "NaCl": "[Na+][Cl-]",
    "KCl": "[K+][Cl-]",
    "CaCl2": "[Ca+2][Cl-][Cl-]",
    "MgCl2": "[Mg+2][Cl-][Cl-]",
    "TIP3P": "[O]",
    "TIP4P": "[O]",
    "SPC": "[O]",
    "SPCE": "[O]",
}


def generate_aa_block(
    lipid: str,
    lipid_data: Dict[str, Any] | None = None,
) -> str:
    """
    Generate one AA lipid species block.

    Signature intentionally matches app.py:
        generate_aa_block(lipid, LIPIDS[lipid])

    lipid_data is accepted for compatibility with the web application.
    The structural definition itself comes only from the verified template
    table and is never inferred from metadata.
    """
    if lipid not in AA_LIPID_TEMPLATES:
        available = ", ".join(sorted(AA_LIPID_TEMPLATES))
        raise ValueError(
            f"No verified AA structural template is available for '{lipid}'. "
            f"Available AA lipid templates: {available}"
        )

    return f"{lipid}=AA|[\n{AA_LIPID_TEMPLATES[lipid]}\n]"


def generate_aa_environment_block(species: str) -> str:
    """Generate one AA solvent or salt species block."""
    if species not in AA_ENVIRONMENT_TEMPLATES:
        raise ValueError(
            f"No AA environment template is available for '{species}'."
        )

    return f"{species}=AA|[\n{AA_ENVIRONMENT_TEMPLATES[species]}\n]"


def build_aa_environment_blocks(
    solvent: str,
    salt_species: str,
) -> List[str]:
    """Return salt and solvent AA species blocks in LiMBS block order."""
    blocks: List[str] = []

    if salt_species:
        blocks.append(generate_aa_environment_block(salt_species))

    if solvent:
        blocks.append(generate_aa_environment_block(solvent))

    return blocks


def available_aa_lipids() -> List[str]:
    return sorted(AA_LIPID_TEMPLATES)


if __name__ == "__main__":
    print(generate_aa_block("DPPC", {}))
    print("|+|")
    print(generate_aa_environment_block("NaCl"))
    print("|+|")
    print(generate_aa_environment_block("TIP3P"))
