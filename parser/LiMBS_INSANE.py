#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LiMBS1.0 to INSANE exporter.

Converts a validated LiMBS CG planar membrane file into an INSANE command line, then optionally executes it to produce membrane.gro and membrane.top.

Workflow
LiMBS .txt  to  parser  to  JSON dict  to  INSANE command  to  .gro / .top

Usage examples
# Print INSANE command only (dry-run):
python limbs_insane.py mySystem_CG.txt --dry-run

# Run INSANE with default output names:
python limbs_insane.py mySystem_CG.txt --insane-script insane.py

# Save intermediate JSON and run INSANE:
python limbs_insane.py mySystem_CG.txt --write-json out.json --insane-script insane.py

# Start from a previously saved LiMBS JSON:
python limbs_insane.py out.json --from-json --insane-script insane.py

Requirements
- limbs_insane.py  (must be in the same directory)
- insane.py reachable at --insane-script path
- Python 3.6+, no extra dependencies
"""

import argparse
import json
import os
import subprocess
import sys

#Import the parser
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

#INSANE script name 
INSANE_SCRIPT_DEFAULT = "insane.py"


#File I/O helpers 

def read_text_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def read_json_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"JSON file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


#LiMBS to JSON

def limbs_to_json(path):
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


#Validation before building INSANE command

def validate_for_insane(data):
    """
    Check that the parsed LiMBS system is compatible with INSANE.
    INSANE supports: CG resolution, planar geometry, count-mode leaflets.
    """
    errors = []

    if data.get("resolution") != "CG":
        errors.append(
            f"Resolution is '{data.get('resolution')}'. "
            "INSANE only supports CG systems."
        )

    if data.get("membrane_type") != "planar":
        errors.append(
            f"Membrane type is '{data.get('membrane_type')}'. "
            "INSANE only supports type:planar. "
            "For vesicles use TS2CG instead."
        )

    leaflets = data.get("leaflets")
    if not leaflets:
        errors.append("Missing leaflet information in parsed system.")
    else:
        if leaflets.get("mode") not in ( "count", "ratio"):
            errors.append(
                f"Leaflet mode is '{leaflets.get('mode')}'. "
                "INSANE supports leaflet:count or leaflet:ratio "
                "for planar CG membranes."
            )
        if not leaflets.get("upper"):
            errors.append("Upper leaflet composition is empty.")
        if not leaflets.get("lower"):
            errors.append("Lower leaflet composition is empty.")

    box = data.get("box")
    if not box or "values" not in box:
        errors.append("Missing box dimensions.")
    elif len(box["values"]) != 3:
        errors.append(
            f"Box has {len(box['values'])} values; expected 3 (x, y, z)."
        )

    if not data.get("solvent"):
        errors.append("Missing solvent field.")
    elif data["solvent"] != "W":
        errors.append(
            f"Solvent is '{data['solvent']}'. "
            "INSANE (MARTINI) requires sol:W."
        )

    if data.get("salt") is None:
        errors.append("Missing salt field.")

    if errors:
        raise RuntimeError(
            "System is not compatible with INSANE:\n  "
            + "\n  ".join(errors)
        )


#Build INSANE command

def build_insane_command(data, insane_script, out_gro, out_top,rand, solr, python_exe):
    """
    Translate the LiMBS JSON dict into an INSANE command-line list.

    INSANE flags used:
        -u  LIPID:count   upper leaflet lipid
        -l  LIPID:count   lower leaflet lipid
        -x / -y / -z      box dimensions (nm)
        -sol              solvent type  (W for MARTINI)
        -salt             salt concentration (M)
        -rand             randomisation seed
        -solr             solute rotational freedom
        -o                output .gro file
        -p                output .top file
    """
    cmd = [python_exe, insane_script]

    upper = data["leaflets"]["upper"]
    lower = data["leaflets"]["lower"]

    # Lipids in alphabetical order for reproducibility
    for lipid in sorted(upper):
        cmd += ["-u", f"{lipid}:{upper[lipid]}"]
    for lipid in sorted(lower):
        cmd += ["-l", f"{lipid}:{lower[lipid]}"]

    box = data["box"]["values"]
    cmd += ["-x", str(box[0])]
    cmd += ["-y", str(box[1])]
    cmd += ["-z", str(box[2])]

    cmd += ["-sol", data["solvent"]]

    salt_conc = data["salt"]["concentration"]
    cmd += ["-salt", str(salt_conc)]

    cmd += ["-rand", str(rand)]
    cmd += ["-solr", str(solr)]

    cmd += ["-o", out_gro]
    cmd += ["-p", out_top]

    return cmd


#Main

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convert a LiMBS CG planar membrane file into an INSANE command "
            "and optionally run it to produce .gro and .top files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    ap.add_argument(
        "input",
        help="LiMBS text file (.txt) or LiMBS JSON file (.json).",
    )
    ap.add_argument(
        "--from-json",
        action="store_true",
        help="Treat input as a LiMBS JSON file instead of a LiMBS text file.",
    )
    ap.add_argument(
        "--write-json",
        metavar="PATH",
        help="Save the intermediate LiMBS JSON to this path.",
    )
    ap.add_argument(
        "--insane-script",
        default=INSANE_SCRIPT_DEFAULT,
        metavar="PATH",
        help=f"Path to the INSANE Python script (default: {INSANE_SCRIPT_DEFAULT}).",
    )
    ap.add_argument(
        "--python-exe",
        default="python3",
        metavar="EXE",
        help="Python executable used to run INSANE (default: python3).",
    )
    ap.add_argument(
        "--out-gro",
        default="membrane.gro",
        metavar="FILE",
        help="Output GROMACS .gro coordinate file (default: membrane.gro).",
    )
    ap.add_argument(
        "--out-top",
        default="membrane.top",
        metavar="FILE",
        help="Output GROMACS .top topology file (default: membrane.top).",
    )
    ap.add_argument(
        "--rand",
        default="0.1",
        metavar="FLOAT",
        help="INSANE -rand value: randomisation seed (default: 0.1).",
    )
    ap.add_argument(
        "--solr",
        default="0.5",
        metavar="FLOAT",
        help="INSANE -solr value: solute rotational freedom (default: 0.5).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the INSANE command without executing it.",
    )

    args = ap.parse_args()

    #Step 1: parse LiMBS
    try:
        if args.from_json:
            print(f"Loading LiMBS JSON: {args.input}")
            data = read_json_file(args.input)
        else:
            print(f"Parsing LiMBS file: {args.input}")
            data = limbs_to_json(args.input)
            print("  Validation: PASSED")

        #Step 2: check INSANE compatibility 
        validate_for_insane(data)
        print("  INSANE compatibility: OK")

        #Step 3: optionally save JSON 
        if args.write_json:
            with open(args.write_json, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)
                fh.write("\n")
            print(f"  JSON saved: {args.write_json}")

        #Step 4: build command
        cmd = build_insane_command(
            data=data,
            insane_script=args.insane_script,
            out_gro=args.out_gro,
            out_top=args.out_top,
            rand=args.rand,
            solr=args.solr,
            python_exe=args.python_exe,
        )

    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)

    #Step 5: print command
    print("\nGenerated INSANE command:")
    print("  " + " ".join(cmd))

    if args.dry_run:
        print("\nDry-run mode: command was not executed.")
        return

    #Step 6: run INSANE
    print("\nRunning INSANE...")
    retcode = subprocess.call(cmd)

    if retcode != 0:
        print(f"\nINSANE failed with return code {retcode}.", file=sys.stderr)
        sys.exit(retcode)

    print("\nINSANE completed successfully.")
    print(f"  Coordinates : {args.out_gro}")
    print(f"  Topology    : {args.out_top}")


if __name__ == "__main__":
    main()
