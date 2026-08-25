#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LiMBS1.0 to TS2CG exporter.

Converts a validated LiMBS CG vesicle specification into a TS2CG input.str
file and can optionally execute TS2CG PLM and PCG to build a vesicle.

The current TS2CG v2 input.str lipid-line format is:

    LipidName   RatioUp   RatioDown   Area/Lipid

Therefore, LiMBS leaflet ratios can be transferred directly, but the
builder-specific area-per-lipid (APL) value for every lipid must be supplied
explicitly by the user. LiMBS does not invent or maintain universal APL values.

Workflow

    LiMBS .txt
         LiMBS parser / validation
         JSON-compatible system dictionary
         TS2CG compatibility check
         user-supplied APL values
         input.str
         TS2CG PLM
         TS2CG PCG
         vesicle coordinate/topology output

Usage examples

# Parse and validate LiMBS only:
python limbs_ts2cg.py mySystem_CG.txt --validate-only

# Write input.str only:
python limbs_ts2cg.py mySystem_CG.txt --apl POPC=0.63 --apl POPE=0.64 --write-str-only

# Dry-run:
python limbs_ts2cg.py mySystem_CG.txt --apl POPC=0.6 --apl POPE=0.64 --dry-run

# Full run with wrapper defaults:
python limbs_ts2cg.py mySystem_CG.txt --apl POPC=0.63 --apl POPE=0.64

Requirements
1. LiMBSv1_parser.py in the same directory
2. TS2CG executable reachable in PATH for non-dry-run execution
3. A compatible triangulated-surface file, e.g. sphere.tsi
4. A compatible Martini lipid library, e.g. Martini3.LIB
5. Explicit APL values for all lipids used in the TS2CG [Lipids List]
6. Python 3.6+; no extra Python dependencies
"""

import argparse
import math
import os
import subprocess
import sys
from typing import Dict, List, Tuple


# Import the LiMBS parser

try:
    from LiMBSv1_parser import (
        LiMBSParser,
        LiMBSError,
        system_to_json,
    )
except ImportError as _e:
    sys.exit(
        "ERROR: Could not import the LiMBS parser.\n"
        "Make sure 'LiMBSv1_parser.py' is in the same directory as this script.\n"
        f"Detail: {_e}"
    )


# TS2CG wrapper settings

TS2CG_COMMAND = "TS2CG"
TSI_FILE_DEFAULT = "sphere.tsi"
MARTINI_LIB_DEFAULT = "Martini3.LIB"
STR_FILE_DEFAULT = "input.str"
OUT_NAME_DEFAULT = "vesicle"

# These are retained as LiMBS wrapper defaults from the existing exporter.
# They remain user-configurable via command-line flags.
BILAYER_THICKNESS_DEFAULT = 3.0
RESCALE_DEFAULT = (4.8, 4.8, 4.8)
MASHNO_DEFAULT = 4
BOND_LENGTH_DEFAULT = 0.2


# File I/O helpers

def read_text_file(path: str) -> str:
    """Read a UTF-8 text file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


# LiMBS -> parsed system dictionary

def limbs_to_data(path: str) -> dict:
    """Parse a LiMBS text file and return the JSON-compatible dictionary."""
    raw = read_text_file(path)
    parser = LiMBSParser()

    try:
        system = parser.parse(raw)
    except LiMBSError as exc:
        raise RuntimeError(f"LiMBS parse error: {exc}") from exc

    data = system_to_json(system)

    if data.get("errors"):
        msgs = "\n  ".join(
            f"[{entry['code']}] {entry['message']}"
            for entry in data["errors"]
        )
        raise RuntimeError(f"LiMBS validation errors:\n  {msgs}")

    return data


# TS2CG compatibility validation

def validate_for_ts2cg(data: dict) -> None:
    """
    Check compatibility with the current LiMBS -> TS2CG adapter.

    Adapter requirements:
    - CG resolution
    - vesicle geometry
    - leaflet:ratio composition mode
    - non-empty upper/outer and lower/inner compositions
    - positive three-dimensional LiMBS box
    - solvent field present
    - salt field present

    The leaflet:ratio restriction belongs to this TS2CG adapter. It is not
    intended to redefine the general LiMBS grammar.
    """
    errors = []

    if data.get("resolution") != "CG":
        errors.append(
            f"Resolution is '{data.get('resolution')}'. "
            "The current LiMBS-TS2CG adapter supports CG systems."
        )

    if data.get("membrane_type") != "vesicle":
        errors.append(
            f"Membrane type is '{data.get('membrane_type')}'. "
            "The current LiMBS-TS2CG adapter supports type:vesicle. "
            "For compatible planar CG systems use the LiMBS-INSANE interface."
        )

    leaflets = data.get("leaflets")
    if not leaflets:
        errors.append("Missing leaflet information in parsed system.")
    else:
        if leaflets.get("mode") != "ratio":
            errors.append(
                f"Leaflet mode is '{leaflets.get('mode')}'. "
                "The current LiMBS-TS2CG adapter requires leaflet:ratio mode."
            )

        if not leaflets.get("upper"):
            errors.append("Upper/outer leaflet composition is empty.")

        if not leaflets.get("lower"):
            errors.append("Lower/inner leaflet composition is empty.")

        for leaflet_name in ("upper", "lower"):
            composition = leaflets.get(leaflet_name, {})
            for lipid, value in composition.items():
                try:
                    numeric = float(value)
                except (TypeError, ValueError):
                    errors.append(
                        f"{leaflet_name.capitalize()} leaflet value for "
                        f"{lipid} is not numeric: {value!r}"
                    )
                    continue

                if not math.isfinite(numeric) or numeric < 0:
                    errors.append(
                        f"{leaflet_name.capitalize()} leaflet value for "
                        f"{lipid} must be finite and non-negative. Found {value!r}."
                    )

    box = data.get("box")
    if not box or "values" not in box:
        errors.append("Missing box dimensions.")
    elif len(box["values"]) != 3:
        errors.append(
            f"Box has {len(box['values'])} values; expected 3 (x, y, z)."
        )
    elif any(v <= 0 for v in box["values"]):
        errors.append(f"Box dimensions must be positive: {box['values']}")

    if not data.get("solvent"):
        errors.append("Missing solvent field.")

    if data.get("salt") is None:
        errors.append("Missing salt field.")

    if errors:
        raise RuntimeError(
            "System is not compatible with the current LiMBS-TS2CG adapter:\n  "
            + "\n  ".join(errors)
        )


def print_ratio_sum_notes(data: dict) -> None:
    """
    Print advisory notes when leaflet values do not sum to approximately 1.

    LiMBS preserves the user's ratio values. This function does not renormalize
    them and does not convert them into lipid counts.
    """
    leaflets = data.get("leaflets", {})

    for leaflet_name, display_name in (
        ("upper", "Upper/outer"),
        ("lower", "Lower/inner"),
    ):
        composition = leaflets.get(leaflet_name, {})
        if not composition:
            continue

        total = sum(float(v) for v in composition.values())

        if abs(total - 1.0) > 1e-6:
            print(
                f"\n  Note: {display_name} leaflet ratio values sum to "
                f"{total:g}, not 1.0."
            )
            print(
                "        Values will be written exactly as supplied; "
                "this exporter does not renormalize them."
            )


# APL handling

def parse_apl_assignments(assignments: List[str]) -> Dict[str, float]:
    """
    Parse repeated --apl LIPID=VALUE assignments.

    Example:
        --apl POPC=0.63 --apl POPE=0.64

    Returns:
        {"POPC": 0.63, "POPE": 0.64}
    """
    apl_map: Dict[str, float] = {}

    for item in assignments:
        if "=" not in item:
            raise ValueError(
                f"Invalid --apl value '{item}'. Expected format LIPID=VALUE, "
                "for example --apl POPC=0.63."
            )

        lipid, raw_value = item.split("=", 1)
        lipid = lipid.strip().upper()
        raw_value = raw_value.strip()

        if not lipid:
            raise ValueError(
                f"Invalid --apl value '{item}': lipid name is empty."
            )

        try:
            value = float(raw_value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid APL value for {lipid}: '{raw_value}' is not numeric."
            ) from exc

        if not math.isfinite(value) or value <= 0:
            raise ValueError(
                f"APL for {lipid} must be a positive finite number. "
                f"Found {raw_value!r}."
            )

        if lipid in apl_map:
            raise ValueError(
                f"APL for lipid '{lipid}' was supplied more than once."
            )

        apl_map[lipid] = value

    return apl_map


def lipid_names_for_ts2cg(data: dict) -> List[str]:
    """Return deterministic union of upper and lower leaflet lipid names."""
    leaflets = data["leaflets"]
    upper = leaflets.get("upper", {})
    lower = leaflets.get("lower", {})
    return sorted(set(upper) | set(lower))


def validate_apl_map(data: dict, apl_map: Dict[str, float]) -> None:
    """Require an explicit positive APL value for every lipid used by TS2CG."""
    required = lipid_names_for_ts2cg(data)

    missing = [lipid for lipid in required if lipid not in apl_map]
    extra = sorted(set(apl_map) - set(required))

    errors = []

    if missing:
        errors.append(
            "Missing APL value(s) for: " + ", ".join(missing)
        )

    if extra:
        errors.append(
            "APL value(s) supplied for lipid(s) not present in this system: "
            + ", ".join(extra)
        )

    if errors:
        raise RuntimeError(
            "TS2CG area-per-lipid input is incomplete or inconsistent:\n  "
            + "\n  ".join(errors)
            + "\n\nSupply one --apl LIPID=VALUE option for every lipid, "
              "for example:\n"
              "  --apl POPC=0.63 --apl POPE=0.64"
        )


# Write TS2CG input.str

def write_input_str(
    data: dict,
    path: str,
    apl_map: Dict[str, float],
) -> None:
    """
    Write a single-domain TS2CG input.str.

    TS2CG v2 lipid lines use:
        LipidName   RatioUp   RatioDown   Area/Lipid

    LiMBS upper leaflet is mapped to TS2CG RatioUp.
    LiMBS lower leaflet is mapped to TS2CG RatioDown.

    APL values are builder-specific and must be supplied explicitly by the
    user through --apl options. No lipid APL values are invented here.
    """
    validate_apl_map(data, apl_map)

    upper = data["leaflets"]["upper"]
    lower = data["leaflets"]["lower"]
    lipid_names = lipid_names_for_ts2cg(data)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("[Lipids List]\n")
        fh.write(";LipidName  RatioUp  RatioDown  Area/Lipid\n")
        fh.write("Domain 0\n")

        for lipid in lipid_names:
            ratio_up = float(upper.get(lipid, 0.0))
            ratio_down = float(lower.get(lipid, 0.0))
            apl = apl_map[lipid]

            fh.write(
                f"{lipid}  {ratio_up:g}  {ratio_down:g}  {apl:g}\n"
            )

        fh.write("End\n")


# Run TS2CG commands

def run_command(command: List[str], dry_run: bool) -> bool:
    """Print and optionally execute a shell command. Return True on success."""
    print("  " + " ".join(command))

    if dry_run:
        print("  dry-run: command not executed.")
        return True

    try:
        completed = subprocess.run(command, check=False)
    except FileNotFoundError:
        print(f"\nError: command not found: {command[0]}")
        return False

    if completed.returncode != 0:
        print(
            f"\nCommand failed with return code {completed.returncode}."
        )
        return False

    return True


def build_ts2cg(
    data: dict,
    apl_map: Dict[str, float],
    dry_run: bool = False,
    tsi_file: str = TSI_FILE_DEFAULT,
    martini_lib: str = MARTINI_LIB_DEFAULT,
    str_file: str = STR_FILE_DEFAULT,
    out_name: str = OUT_NAME_DEFAULT,
    bilayer_thickness: float = BILAYER_THICKNESS_DEFAULT,
    rescale: Tuple[float, float, float] = RESCALE_DEFAULT,
    mashno: int = MASHNO_DEFAULT,
    bond_length: float = BOND_LENGTH_DEFAULT,
) -> bool:
    """
    Write input.str, run TS2CG PLM, then run TS2CG PCG.

    In dry-run mode, commands are printed but external programs are not run.
    """
    validate_apl_map(data, apl_map)

    if not dry_run:
        if not os.path.isfile(tsi_file):
            print(f"\nError: TS2CG geometry file not found: {tsi_file}")
            return False

        if not os.path.isfile(martini_lib):
            print(f"\nError: Martini lipid library not found: {martini_lib}")
            return False

    write_input_str(data, str_file, apl_map)
    print(f"\n[TS2CG input.str]\n  Wrote: {str_file}")

    plm_cmd = [
        TS2CG_COMMAND,
        "PLM",
        "-TSfile",
        tsi_file,
        "-bilayerThickness",
        str(bilayer_thickness),
        "-rescalefactor",
        str(rescale[0]),
        str(rescale[1]),
        str(rescale[2]),
        "-Mashno",
        str(mashno),
    ]

    pcg_cmd = [
        TS2CG_COMMAND,
        "PCG",
        "-dts",
        "point",
        "-str",
        str_file,
        "-Bondlength",
        str(bond_length),
        "-LLIB",
        martini_lib,
        "-defout",
        out_name,
    ]

    print("\n[Step 1] TS2CG PLM")
    if not run_command(plm_cmd, dry_run):
        return False

    print("\n[Step 2] TS2CG PCG")
    if not run_command(pcg_cmd, dry_run):
        return False

    if dry_run:
        print("\nDry-run completed successfully.")
    else:
        print("\nTS2CG completed successfully.")
        print(f"  Output prefix : {out_name}")

    return True


# Summary printer

def print_summary(data: dict) -> None:
    """Print a concise parsed-system summary."""
    print("\nParsed LiMBS vesicle system:")
    print(f"  Resolution    : {data.get('resolution')}")
    print(f"  Geometry      : {data.get('membrane_type')}")
    print(f"  Box           : {data.get('box', {}).get('values')} nm")
    print(f"  Solvent       : {data.get('solvent')}")

    salt = data.get("salt", {})
    if salt:
        print(
            f"  Salt          : {salt.get('concentration')} "
            f"{salt.get('unit', 'M')} {salt.get('species_text', '')}"
        )
    else:
        print("  Salt          : None")

    leaflets = data.get("leaflets", {})
    print(f"  Leaflet mode  : {leaflets.get('mode')}")
    print(f"  Outer (upper) : {leaflets.get('upper')}")
    print(f"  Inner (lower) : {leaflets.get('lower')}")

    warnings = data.get("warnings", [])
    if warnings:
        print("\nParser warnings:")
        for warning in warnings:
            if isinstance(warning, dict):
                code = warning.get("code", "warning")
                message = warning.get("message", warning)
                print(f"  - [{code}] {message}")
            else:
                print(f"  - {warning}")


# Main

def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Convert a validated LiMBS CG vesicle specification into a "
            "TS2CG input.str and optionally run TS2CG PLM + PCG."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    ap.add_argument(
        "input",
        help="LiMBS text file (.txt).",
    )

    ap.add_argument(
        "--validate-only",
        action="store_true",
        help=(
            "Parse and validate LiMBS and check TS2CG adapter compatibility; "
            "write nothing and do not require APL values."
        ),
    )

    ap.add_argument(
        "--write-str-only",
        action="store_true",
        help="Write input.str only; do not run TS2CG.",
    )

    ap.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Write input.str and print TS2CG commands without executing them."
        ),
    )

    ap.add_argument(
        "--apl",
        action="append",
        default=[],
        metavar="LIPID=VALUE",
        help=(
            "Builder-specific TS2CG area per lipid in nm^2. "
            "Repeat once for every lipid, e.g. "
            "--apl POPC=0.63 --apl POPE=0.64."
        ),
    )

    ap.add_argument(
        "--tsi",
        default=TSI_FILE_DEFAULT,
        metavar="FILE",
        help=(
            "TS2CG triangulated-surface geometry file "
            f"(LiMBS wrapper default: {TSI_FILE_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--lib",
        default=MARTINI_LIB_DEFAULT,
        metavar="FILE",
        help=(
            "TS2CG Martini lipid library file "
            f"(LiMBS wrapper default: {MARTINI_LIB_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--str-file",
        default=STR_FILE_DEFAULT,
        metavar="FILE",
        help=(
            "Output TS2CG composition file "
            f"(default: {STR_FILE_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--out-name",
        default=OUT_NAME_DEFAULT,
        metavar="NAME",
        help=(
            "TS2CG output-file prefix "
            f"(default: {OUT_NAME_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--bilayer-thickness",
        type=float,
        default=BILAYER_THICKNESS_DEFAULT,
        metavar="FLOAT",
        help=(
            "TS2CG PLM bilayer thickness in nm "
            f"(LiMBS wrapper default: {BILAYER_THICKNESS_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--rescale",
        type=float,
        nargs=3,
        default=RESCALE_DEFAULT,
        metavar=("X", "Y", "Z"),
        help=(
            "TS2CG PLM rescale factors "
            f"(LiMBS wrapper default: "
            f"{RESCALE_DEFAULT[0]} {RESCALE_DEFAULT[1]} {RESCALE_DEFAULT[2]})."
        ),
    )

    ap.add_argument(
        "--mashno",
        type=int,
        default=MASHNO_DEFAULT,
        metavar="INT",
        help=(
            "TS2CG PLM Mashno value "
            f"(LiMBS wrapper default: {MASHNO_DEFAULT})."
        ),
    )

    ap.add_argument(
        "--bond-length",
        type=float,
        default=BOND_LENGTH_DEFAULT,
        metavar="FLOAT",
        help=(
            "TS2CG PCG initial bond-length guess in nm "
            f"(LiMBS wrapper default: {BOND_LENGTH_DEFAULT})."
        ),
    )

    args = ap.parse_args()

    # Step 1: Parse and validate LiMBS

    try:
        print(f"Parsing LiMBS file: {args.input}")
        data = limbs_to_data(args.input)
        print("  Validation: PASSED")
    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    print_summary(data)

    # Step 2: Check TS2CG adapter compatibility

    try:
        validate_for_ts2cg(data)
        print("\n  TS2CG compatibility: OK")
        print_ratio_sum_notes(data)
    except RuntimeError as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.validate_only:
        print("\nValidation and TS2CG compatibility checks passed. Exiting.")
        sys.exit(0)

    # Step 3: Parse and validate explicit APL inputs

    try:
        apl_map = parse_apl_assignments(args.apl)
        validate_apl_map(data, apl_map)
    except (ValueError, RuntimeError) as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    print("\nTS2CG area-per-lipid values:")
    for lipid in lipid_names_for_ts2cg(data):
        print(f"  {lipid}: {apl_map[lipid]:g} nm^2")

    # Step 4: Write input.str only

    if args.write_str_only:
        try:
            write_input_str(data, args.str_file, apl_map)
        except Exception as exc:
            print(f"\nError: {exc}", file=sys.stderr)
            sys.exit(1)

        print(f"\nWrote {args.str_file}")
        sys.exit(0)

    # Step 5: Build / dry-run TS2CG

    success = build_ts2cg(
        data=data,
        apl_map=apl_map,
        dry_run=args.dry_run,
        tsi_file=args.tsi,
        martini_lib=args.lib,
        str_file=args.str_file,
        out_name=args.out_name,
        bilayer_thickness=args.bilayer_thickness,
        rescale=tuple(args.rescale),
        mashno=args.mashno,
        bond_length=args.bond_length,
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
