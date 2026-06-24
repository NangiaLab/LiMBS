#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LiMBS1.0 to TS2CG exporter.

Converts a validated LiMBS CG vesicle file into a TS2CG input.str,
then optionally executes TS2CG PLM and PCG to produce a vesicle
coordinate and topology file.

Workflow

LiMBS .txt to  parser to  system dict  to  input.str  to  TS2CG PLM to  TS2CG PCG to  vesicle output

Usage examples

# Parse and validate only:
python limbs_ts2cg.py mySystem_CG.txt --validate-only

# Write input.str only:
python limbs_ts2cg.py mySystem_CG.txt --write-str-only

# Dry-run:
python limbs_ts2cg.py mySystem_CG.txt --dry-run

# Full run with default settings:
python limbs_ts2cg.py mySystem_CG.txt

# Full run with custom geometry and library:
python limbs_ts2cg.py mySystem_CG.txt --tsi mySphere.tsi --lib MyMartini.LIB --out-name myvesicle

Requirements
1. LiMBSv1_parser.py  (must be in the same directory)
2. TS2CG executable reachable in PATH
3. sphere.tsi  (TS2CG triangulated surface geometry file)
4. Martini3.LIB  (TS2CG Martini lipid library)
5. Python 3.6+, no extra dependencies
"""

import argparse
import os
import subprocess
import sys
from typing import List, Tuple

#Import the LiMBS parser
try:
    from LiMBSv1_parser import (
        LiMBSParser,
        LiMBSError,
        system_to_json,
    )
except ImportError as _e:
    sys.exit(
        f"ERROR: Could not import the LiMBS parser.\n"
        f"Make sure 'LiMBSv1_parser.py' is in "
        f"the same directory as this script.\nDetail: {_e}"
    )

#TS2CG settings
TS2CG_COMMAND = "TS2CG"
TSI_FILE_DEFAULT = "sphere.tsi"
MARTINI_LIB_DEFAULT = "Martini3.LIB"
STR_FILE_DEFAULT = "input.str"
OUT_NAME_DEFAULT = "vesicle"


#File I/O helpers

def read_text_file(path: str) -> str:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


#LiMBS to  parsed system

def limbs_to_data(path: str) -> dict:
    """Parse a LiMBS text file and return the JSON-compatible dict."""
    raw = read_text_file(path)
    parser = LiMBSParser()

    try:
        system = parser.parse(raw)
    except LiMBSError as exc:
        raise RuntimeError(f"LiMBS parse error: {exc}") from exc

    data = system_to_json(system)

    if data.get("errors"):
        msgs = "\n  ".join(
            f"[{e['code']}] {e['message']}" for e in data["errors"]
        )
        raise RuntimeError(f"LiMBS validation errors:\n  {msgs}")

    return data


#Validation

def validate_for_ts2cg(data: dict) -> None:
    """
    Check that the parsed LiMBS system is compatible with TS2CG.
    TS2CG supports: CG resolution, vesicle geometry, ratio-mode leaflets.
    """
    errors = []

    if data.get("resolution") != "CG":
        errors.append(
            f"Resolution is '{data.get('resolution')}'. "
            "TS2CG only supports CG systems."
        )

    if data.get("membrane_type") != "vesicle":
        errors.append(
            f"Membrane type is '{data.get('membrane_type')}'. "
            "TS2CG only supports type:vesicle. "
            "For planar membranes use INSANE instead."
        )

    leaflets = data.get("leaflets")
    if not leaflets:
        errors.append("Missing leaflet information in parsed system.")
    else:
        if leaflets.get("mode") != "ratio":
            errors.append(
                f"Leaflet mode is '{leaflets.get('mode')}'. "
                "TS2CG requires leaflet:ratio mode "
                "(use leaflet:count only for planar membranes with INSANE)."
            )
        if not leaflets.get("upper"):
            errors.append("Upper leaflet (outer monolayer) composition is empty.")
        if not leaflets.get("lower"):
            errors.append("Lower leaflet (inner monolayer) composition is empty.")

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
            "System is not compatible with TS2CG:\n  "
            + "\n  ".join(errors)
        )


#Write TS2CG input.str

def write_input_str(data: dict, path: str) -> None:
    """
    Write a TS2CG input.str for a single-domain vesicle.

    Both monolayers are written under Domain 0 because the default
    sphere.tsi geometry contains only Domain 0 points.

    Format:
        [Lipids List]
        Domain 0
        DPPC 1 1 0.64   # upper/outer monolayer  (leaflet id 1)
        DPPC 2 1 0.64   # lower/inner monolayer  (leaflet id 2)
        End
    """
    upper = data["leaflets"]["upper"]
    lower = data["leaflets"]["lower"]

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("[Lipids List]\n")
        fh.write("Domain 0\n")

        for lipid, value in upper.items():
            fh.write(f"{lipid} 1 1 {value:g}\n")

        for lipid, value in lower.items():
            fh.write(f"{lipid} 2 1 {value:g}\n")

        fh.write("End\n")


#Run TS2CG commands

def run_command(command: List[str], dry_run: bool) -> bool:
    """Print and optionally execute a shell command. Returns True on success."""
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
        print(f"\nCommand failed with return code {completed.returncode}.")
        return False

    return True


def build_ts2cg(
    data: dict,
    dry_run: bool = False,
    tsi_file: str = TSI_FILE_DEFAULT,
    martini_lib: str = MARTINI_LIB_DEFAULT,
    str_file: str = STR_FILE_DEFAULT,
    out_name: str = OUT_NAME_DEFAULT,
    bilayer_thickness: float = 3.0,
    rescale: Tuple[float, float, float] = (4.8, 4.8, 4.8),
    mashno: int = 4,
    bond_length: float = 0.2,
) -> bool:
    """Write input.str and run TS2CG PLM then PCG."""

    if not dry_run:
        if not os.path.isfile(tsi_file):
            print(f"\nError: TS2CG geometry file not found: {tsi_file}")
            return False
        if not os.path.isfile(martini_lib):
            print(f"\nError: Martini lipid library not found: {martini_lib}")
            return False

    #Write input.str
    write_input_str(data, str_file)
    print(f"\n[TS2CG input.str]\n  Wrote: {str_file}")

    plm_cmd = [
        TS2CG_COMMAND, "PLM",
        "-TSfile", tsi_file,
        "-bilayerThickness", str(bilayer_thickness),
        "-rescalefactor", str(rescale[0]), str(rescale[1]), str(rescale[2]),
        "-Mashno", str(mashno),
    ]

    pcg_cmd = [
        TS2CG_COMMAND, "PCG",
        "-dts", "point",
        "-str", str_file,
        "-Bondlength", str(bond_length),
        "-LLIB", martini_lib,
        "-defout", out_name,
    ]

    print("\n[Step 1] TS2CG PLM")
    if not run_command(plm_cmd, dry_run):
        return False

    print("\n[Step 2] TS2CG PCG")
    if not run_command(pcg_cmd, dry_run):
        return False

    print(f"\nTS2CG completed successfully.")
    print(f"  Output prefix : {out_name}")

    return True


#Summary printer 

def print_summary(data: dict) -> None:
    print("\nParsed LiMBS vesicle system:")
    print(f"  Resolution    : {data.get('resolution')}")
    print(f"  Geometry      : {data.get('membrane_type')}")
    print(f"  Box           : {data.get('box', {}).get('values')} nm")
    print(f"  Solvent       : {data.get('solvent')}")
    salt = data.get("salt", {})
    print(f"  Salt          : {salt.get('concentration')} M  {salt.get('ions', '')}")
    leaflets = data.get("leaflets", {})
    print(f"  Outer (upper) : {leaflets.get('upper')}")
    print(f"  Inner (lower) : {leaflets.get('lower')}")

    warnings = data.get("warnings", [])
    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  - {w.get('message', w)}")


#Main

def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Convert a LiMBS CG vesicle file into a TS2CG input.str "
            "and optionally run TS2CG PLM + PCG to build the vesicle."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    ap.add_argument("input", help="LiMBS text file (.txt).")
    ap.add_argument("--validate-only", action="store_true",
                    help="Parse and validate only; write nothing.")
    ap.add_argument("--write-str-only", action="store_true",
                    help="Write input.str only; do not run TS2CG.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print TS2CG commands without executing them.")

    ap.add_argument("--tsi", default=TSI_FILE_DEFAULT, metavar="FILE",
                    help=f"TS2CG triangulated surface geometry file (default: {TSI_FILE_DEFAULT}).")
    ap.add_argument("--lib", default=MARTINI_LIB_DEFAULT, metavar="FILE",
                    help=f"TS2CG Martini lipid library file (default: {MARTINI_LIB_DEFAULT}).")
    ap.add_argument("--str-file", default=STR_FILE_DEFAULT, metavar="FILE",
                    help=f"Output TS2CG composition file (default: {STR_FILE_DEFAULT}).")
    ap.add_argument("--out-name", default=OUT_NAME_DEFAULT, metavar="NAME",
                    help=f"TS2CG output file prefix (default: {OUT_NAME_DEFAULT}).")

    ap.add_argument("--bilayer-thickness", type=float, default=3.0, metavar="FLOAT",
                    help="TS2CG PLM bilayer thickness in nm (default: 3.0).")
    ap.add_argument("--rescale", type=float, nargs=3, default=(4.8, 4.8, 4.8),
                    metavar=("X", "Y", "Z"),
                    help="TS2CG PLM rescale factors (default: 4.8 4.8 4.8).")
    ap.add_argument("--mashno", type=int, default=4, metavar="INT",
                    help="TS2CG PLM Mashno value (default: 4).")
    ap.add_argument("--bond-length", type=float, default=0.2, metavar="FLOAT",
                    help="TS2CG PCG bond length in nm (default: 0.2).")

    args = ap.parse_args()

    #Step 1: parse and validate LiMBS 
    try:
        print(f"Parsing LiMBS file: {args.input}")
        data = limbs_to_data(args.input)
        print("  Validation: PASSED")
    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    print_summary(data)

    if args.validate_only:
        print("\nValidation passed. Exiting.")
        sys.exit(0)

    #Step 2: check TS2CG compatibility
    try:
        validate_for_ts2cg(data)
        print("\n  TS2CG compatibility: OK")
    except RuntimeError as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    #Step 3: write input.str only
    if args.write_str_only:
        write_input_str(data, args.str_file)
        print(f"\nWrote {args.str_file}")
        sys.exit(0)

    #Step 4: build with TS2CG 
    success = build_ts2cg(
        data=data,
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
