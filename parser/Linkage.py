#!/usr/bin/env python3
"""
LiMBS Linkage Module
Enrich a validated LiMBS file with external lipid-database identifiers.

Outputs
Enriched JSON
Optional TSV and/or CSV tables 

Examples
python linkage_robust.py system.txt
python linkage_robust.py system.txt --parser LiMBSv1_parser.py
python linkage_robust.py system.txt --format tsv
python linkage_robust.py system.txt --format all
python linkage_robust.py system.txt --out results.json
python linkage_robust.py system.txt --summary

"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# Linkage-table metadata

TABLE_VERSION = "1.0"
TABLE_DATE = "2026-06"


# External identifier table
# Confidence levels:
#   full       = LMID + ChEBI ID + InChIKey are all present
#   partial    = LMID is present but one or more other identifiers are missing
#   class-only = only a lipid-class identifier is available
#   unknown    = lipid is not present in this static linkage table
# IMPORTANT:
# The values below are treated as curated static data. This script validates
# identifier *format* only; it does not verify database records over the web.

EXTERNAL_IDS: dict[str, dict[str, Any]] = {
    "POPA": {
        "lmid": "LMGP10010002",
        "chebi_id": "CHEBI:75966",
        "inchikey": "HQKQBZAQYXBDGN-NSHDSACASA-N",
        "confidence": "full",
    },
    "PIPA": {
        "lmid": "LMGP10010023",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PAPA": {
        "lmid": "LMGP10010043",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DPPC": {
        "lmid": "LMGP01010564",
        "chebi_id": "CHEBI:16700",
        "inchikey": "IIZPIQOKICCWNB-UHFFFAOYSA-N",
        "confidence": "full",
    },
    "POPC": {
    "lmid": "LMGP01010005",
    "chebi_id": "CHEBI:73001",
    "inchikey": "WTJKGGKOPKCXLL-VYOBOKEXSA-N",
    "pubchem_cid": 5497103,
    "confidence": "full",
    },
    "DOPC": {
    "lmid": "LMGP01010890",
    "chebi_id": "CHEBI:74669",
    "inchikey": "SNKAWJBJQDLSFF-NVKMUCNASA-N",
    "pubchem_cid": 10350317,
    "confidence": "full",
    },
    "PAPC": {
        "lmid": "LMGP01012254",
        "chebi_id": "CHEBI:77103",
        "inchikey": "RZRNAYUHWVFMIP-NSHDSACASA-N",
        "confidence": "full",
    },
    "PUPC": {
        "lmid": "LMGP01010004",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DSPC": {
        "lmid": "LMGP01010006",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "POPE": {
        "lmid": "LMGP02010009",
        "chebi_id": "CHEBI:74544",
        "inchikey": "WTJKGGKOPKCXLL-RYUDHWBXSA-N",
        "pubchem_cid": 5283496,
        "confidence": "full",
    },
    "DOPE": {
        "lmid": "LMGP02010052",
        "chebi_id": "CHEBI:74539",
        "inchikey": "JYIIAQXHSRFMRE-NSHDSACASA-N",
        "confidence": "full",
    },
    "DPPE": {
        "lmid": "LMGP02010037",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DAPE": {
        "lmid": "LMGP02010010",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PUPE": {
        "lmid": "LMGP02010004",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "POPS": {
        "lmid": "LMGP03010024",
        "chebi_id": "CHEBI:77517",
        "inchikey": "WTJKGGKOPKCXLL-NSHDSACASA-N",
        "confidence": "full",
    },
    "DOPS": {
        "lmid": "LMGP03010001",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PIPS": {
        "lmid": "LMGP03010006",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PQPS": {
        "lmid": "LMGP03010007",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "POPI": {
        "lmid": "LMGP06010005",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PIPI": {
        "lmid": "LMGP06010006",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PAPI": {
        "lmid": "LMGP06010003",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PUPI": {
        "lmid": "LMGP06010004",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "POPG": {
        "lmid": "LMGP04010005",
        "chebi_id": "CHEBI:75926",
        "inchikey": "JRMCHWNYKPMVIG-NSHDSACASA-N",
        "pubchem_cid": 5497103,
        "confidence": "full",
    },
    "DOPG": {
        "lmid": "LMGP04010001",
        "chebi_id": "CHEBI:75929",
        "inchikey": "XPKDNKHAZMTKAK-NSHDSACASA-N",
        "confidence": "full",
    },
    "PODG": {
        "lmid": "LMGP17010001",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DPSM": {
        "lmid": "LMSP01010001",
        "chebi_id": "CHEBI:17076",
        "inchikey": "TZCPCKNHTURCKP-NJSFGKDESA-N",
        "confidence": "full",
    },
    "CHOL": {
        "lmid": "LMST01010001",
        "chebi_id": "CHEBI:16113",
        "inchikey": "HVYWMOMLDIMFJA-DPAQBDIFSA-N",
        "pubchem_cid": 5997,
        "confidence": "full",
    },
    "DPCE": {
        "lmid": "LMSP02010001",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DXCE": {
        "lmid": "LMSP02010002",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PNCE": {
        "lmid": "LMSP02010003",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "XNCE": {
        "lmid": "LMSP02010004",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DPG1": {
        "lmid": "LMGL02010001",
        "chebi_id": "CHEBI:28009",
        "inchikey": None,
        "confidence": "partial",
    },
    "DXG1": {
        "lmid": "LMGL02010002",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "PNG1": {
        "lmid": "LMGL02010003",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "XNG1": {
        "lmid": "LMGL02010004",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
    },
    "DPG3": {
        "lmid": "LMGL02020001",
        "chebi_id": "CHEBI:28217",
        "inchikey": None,
        "confidence": "partial",
    },
    "PVPE": {
        "lmid": "LMGP02010003",
        "chebi_id": None,
        "inchikey": None,
        "confidence": "partial",
        "note": "Vaccenic acid PE (18:1 Delta11); found in bacterial membranes",
    },
    "PVCL2": {
        "lmid": "LMGP12010005",
        "chebi_id": "CHEBI:28494",
        "inchikey": "DSNRWDQKZIEDDB-UHFFFAOYSA-N",
        "confidence": "full",
        "note": "Cardiolipin with vaccenic acid tails",
    },
}

TABLE_COVERAGE = f"{len(EXTERNAL_IDS)} predefined lipid entries"


# Identifier normalization and format validation

LMID_PATTERN = re.compile(r"^LM[A-Z]{2}\d{8}$")
INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEBI_PATTERN = re.compile(r"^CHEBI:\d+$")


def add_lipidmaps_urls() -> None:
    """Add a LIPID MAPS URL to every entry containing a valid LMID."""
    for entry in EXTERNAL_IDS.values():
        lmid = entry.get("lmid")
        if lmid and LMID_PATTERN.fullmatch(str(lmid)):
            entry["lipidmaps_url"] = (
                f"https://www.lipidmaps.org/databases/lmsd/{lmid}"
            )
        elif not lmid and entry.get("lipidmaps_class_url"):
            entry["lipidmaps_url"] = entry["lipidmaps_class_url"]


add_lipidmaps_urls()


def validate_lmid(lmid: str | None) -> bool:
    """Return True when an LMID has the expected textual format."""
    return bool(lmid and LMID_PATTERN.fullmatch(lmid))


def validate_inchikey(inchikey: str | None) -> bool:
    """Return True when an InChIKey has the expected textual format."""
    return bool(inchikey and INCHIKEY_PATTERN.fullmatch(inchikey))


def validate_chebi(chebi: str | None) -> bool:
    """Return True when a ChEBI identifier has the expected textual format."""
    return bool(chebi and CHEBI_PATTERN.fullmatch(chebi))


def validate_entry(name: str, entry: dict[str, Any]) -> list[str]:
    """
    Validate identifier formatting for one linkage-table entry.

    This checks syntax/format only. It does not confirm that an identifier
    corresponds to the intended molecule in an external database.
    """
    warnings: list[str] = []

    lmid = entry.get("lmid")
    if lmid and not validate_lmid(str(lmid)):
        warnings.append(
            f"{name}: LMID '{lmid}' does not match expected format LMxx########"
        )

    inchikey = entry.get("inchikey")
    if inchikey and not validate_inchikey(str(inchikey)):
        warnings.append(
            f"{name}: InChIKey '{inchikey}' does not match expected format"
        )

    chebi = entry.get("chebi_id")
    if chebi and not validate_chebi(str(chebi)):
        warnings.append(
            f"{name}: ChEBI ID '{chebi}' does not match expected format CHEBI:#####"
        )

    return warnings


# Parser discovery and execution

PARSER_NAMES = ["LiMBSv1_parser.py"]


def find_parser(limbs_file: str, explicit_path: str | None = None) -> str:
    """Find the LiMBS reference parser."""
    if explicit_path:
        explicit = Path(explicit_path).expanduser().resolve()
        if explicit.is_file():
            return str(explicit)
        raise FileNotFoundError(f"Parser not found at: {explicit_path}")

    input_path = Path(limbs_file).expanduser().resolve()
    search_dirs = [input_path.parent, Path.cwd()]

    for directory in search_dirs:
        for parser_name in PARSER_NAMES:
            candidate = directory / parser_name
            if candidate.is_file():
                return str(candidate.resolve())

    raise FileNotFoundError(
        "LiMBS parser not found.\n"
        "Place LiMBSv1_parser.py in the same directory as the input file or "
        "current working directory,\n"
        "or provide it explicitly with --parser /path/to/LiMBSv1_parser.py"
    )


def parse_limbs_to_json(limbs_file: str, parser_path: str) -> dict[str, Any]:
    """
    Run the LiMBS parser with --json and return the parsed JSON object.
    """
    input_path = Path(limbs_file)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input LiMBS file not found: {limbs_file}")

    result = subprocess.run(
        [sys.executable, parser_path, str(input_path), "--json"],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print("Parser returned an error.", file=sys.stderr)
        if result.stdout.strip():
            print(result.stdout.strip(), file=sys.stderr)
        if result.stderr.strip():
            print(result.stderr.strip(), file=sys.stderr)
        raise SystemExit(result.returncode or 1)

    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        print(
            f"Failed to parse JSON output from parser: {exc}",
            file=sys.stderr,
        )
        print(
            "Raw parser output (first 500 characters):",
            result.stdout[:500],
            file=sys.stderr,
        )
        raise SystemExit(1)

    if not isinstance(data, dict):
        print("Parser JSON output is not a JSON object.", file=sys.stderr)
        raise SystemExit(1)

    if not data.get("valid", True):
        print(
            "LiMBS validation failed. Fix these errors before running linkage:",
            file=sys.stderr,
        )
        for err in data.get("errors", []):
            if isinstance(err, dict):
                print(
                    f"  {err.get('code', 'ERROR')}: {err.get('message', '')}",
                    file=sys.stderr,
                )
            else:
                print(f"  {err}", file=sys.stderr)
        raise SystemExit(1)

    return data


# Linkage enrichment

def build_case_insensitive_lookup() -> dict[str, str]:
    """Map uppercase lipid names to the canonical key used in EXTERNAL_IDS."""
    return {key.upper(): key for key in EXTERNAL_IDS}


def enrich_with_external_ids(limbs_json: dict[str, Any]) -> dict[str, Any]:
    """
    Add curated external identifiers and a linkage summary to parsed LiMBS JSON.
    """
    enriched = deepcopy(limbs_json)

    lipids = enriched.get("lipids", {})
    if not isinstance(lipids, dict):
        raise ValueError(
            "Parser JSON does not contain the expected 'lipids' dictionary."
        )

    counts = {
        "full": 0,
        "partial": 0,
        "class-only": 0,
        "unknown": 0,
    }
    format_warnings: list[str] = []
    lookup = build_case_insensitive_lookup()

    for name, lipid_data in lipids.items():
        if not isinstance(lipid_data, dict):
            lipid_data = {"raw_value": lipid_data}
            lipids[name] = lipid_data

        canonical_key = lookup.get(str(name).upper())

        if canonical_key is not None:
            entry = deepcopy(EXTERNAL_IDS[canonical_key])
            format_warnings.extend(validate_entry(str(name), entry))

            confidence = entry.get("confidence", "partial")
            if confidence not in counts:
                format_warnings.append(
                    f"{name}: unrecognized confidence value '{confidence}'; "
                    "counted as partial"
                )
                confidence = "partial"
                entry["confidence"] = confidence

            counts[confidence] += 1
            lipid_data["external_ids"] = entry

        else:
            lipid_data["external_ids"] = {
                "lmid": None,
                "chebi_id": None,
                "inchikey": None,
                "pubchem_cid": None,
                "confidence": "unknown",
                "note": (
                    f"'{name}' not found in linkage table "
                    f"(v{TABLE_VERSION}). Add or curate external identifiers "
                    "manually if required."
                ),
            }
            counts["unknown"] += 1

    enriched["lipids"] = lipids
    enriched["linkage_summary"] = {
        "total_lipids": len(lipids),
        "full": counts["full"],
        "partial": counts["partial"],
        "class_only": counts["class-only"],
        "unknown": counts["unknown"],
        "format_warnings": format_warnings,
        "strategy": "static_table",
        "table_version": TABLE_VERSION,
        "table_date": TABLE_DATE,
        "table_coverage": TABLE_COVERAGE,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    return enriched


# Reporting

def print_linkage_report(enriched: dict[str, Any], verbose: bool = True) -> None:
    """Print a human-readable linkage report."""
    print("\n" + "=" * 68)
    print("  LiMBS LINKAGE REPORT")
    print("=" * 68)

    print(f"  Resolution : {enriched.get('resolution')}")
    print(f"  Type       : {enriched.get('membrane_type')}")

    box = enriched.get("box", {})
    if isinstance(box, dict):
        print(f"  Box        : {box.get('values')} {box.get('unit')}")
    else:
        print(f"  Box        : {box}")

    print(f"  Valid      : {enriched.get('valid')}")

    summary = enriched.get("linkage_summary", {})
    total = summary.get("total_lipids", 0)

    print(
        f"\n  Linkage table : v{summary.get('table_version')} "
        f"({summary.get('table_date')})"
    )
    print(f"  Coverage      : {summary.get('table_coverage')}")

    print(f"\n  Lipid linkage summary ({total} total):")
    print(f"    Full       : {summary.get('full', 0)}")
    print(f"    Partial    : {summary.get('partial', 0)}")
    print(f"    Class-only : {summary.get('class_only', 0)}")
    print(f"    Unknown    : {summary.get('unknown', 0)}")

    warnings = summary.get("format_warnings") or []
    if warnings:
        print("\n  Format warnings:")
        for warning in warnings:
            print(f"    - {warning}")

    if not verbose:
        print("=" * 68)
        return

    print("\n  Per-lipid external IDs:")
    print("  " + "-" * 64)

    badge_map = {
        "full": "[OK]",
        "partial": "[PARTIAL]",
        "class-only": "[CLASS]",
        "unknown": "[UNKNOWN]",
    }

    for name, data in enriched.get("lipids", {}).items():
        ids = data.get("external_ids", {}) if isinstance(data, dict) else {}
        confidence = ids.get("confidence", "unknown")
        badge = badge_map.get(confidence, "[?]")

        print(f"\n  {badge} {name} ({confidence})")

        if ids.get("note"):
            print(f"    note          : {ids['note']}")

        for key in ("lmid", "chebi_id", "pubchem_cid", "inchikey"):
            value = ids.get(key)
            if value not in (None, ""):
                print(f"    {key:<14}: {value}")

        url = ids.get("lipidmaps_url")
        if url:
            print(f"    {'url':<14}: {url}")

    print("\n" + "=" * 68)


# Output helpers

def default_output_paths(input_file: str) -> dict[str, Path]:
    """Return default JSON/TSV/CSV output paths beside the input file."""
    input_path = Path(input_file)
    stem_path = input_path.with_suffix("")
    return {
        "json": Path(f"{stem_path}_linked.json"),
        "tsv": Path(f"{stem_path}_linked.tsv"),
        "csv": Path(f"{stem_path}_linked.csv"),
    }


def ensure_parent_directory(path: Path) -> None:
    """Create the parent directory when needed."""
    path.parent.mkdir(parents=True, exist_ok=True)


def check_no_overwrite(paths: list[Path]) -> None:
    """Exit if any requested output path already exists."""
    existing = [path for path in paths if path.exists()]
    if existing:
        print(
            "Error: --no-overwrite was requested and output file(s) already exist:",
            file=sys.stderr,
        )
        for path in existing:
            print(f"  {path}", file=sys.stderr)
        raise SystemExit(1)


def write_json(enriched: dict[str, Any], out_path: Path) -> None:
    ensure_parent_directory(out_path)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(enriched, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"  JSON saved to : {out_path}")


def iter_linkage_rows(enriched: dict[str, Any]):
    """Yield normalized linkage rows."""
    for name, data in enriched.get("lipids", {}).items():
        ids = data.get("external_ids", {}) if isinstance(data, dict) else {}
        yield {
            "LiMBS_name": name,
            "LIPID_MAPS_ID": ids.get("lmid")
            or ids.get("lipidmaps_class")
            or "",
            "ChEBI_ID": ids.get("chebi_id") or "",
            "InChIKey": ids.get("inchikey") or "",
            "PubChem_CID": ids.get("pubchem_cid") or "",
            "Confidence": ids.get("confidence", "unknown"),
            "URL": ids.get("lipidmaps_url") or "",
            "Note": ids.get("note") or "",
        }


def write_tsv(enriched: dict[str, Any], out_path: Path) -> None:
    ensure_parent_directory(out_path)
    fieldnames = [
        "LiMBS_name",
        "LIPID_MAPS_ID",
        "ChEBI_ID",
        "InChIKey",
        "PubChem_CID",
        "Confidence",
        "URL",
        "Note",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(iter_linkage_rows(enriched))
    print(f"  TSV table saved to: {out_path}")


def write_csv(enriched: dict[str, Any], out_path: Path) -> None:
    ensure_parent_directory(out_path)
    fieldnames = [
        "LiMBS_name",
        "LIPID_MAPS_ID",
        "ChEBI_ID",
        "InChIKey",
        "PubChem_CID",
        "Confidence",
        "URL",
        "Note",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(iter_linkage_rows(enriched))
    print(f"  CSV table saved to: {out_path}")


# Command-line interface

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Enrich a validated LiMBS file with curated external lipid-database "
            "identifiers."
        )
    )
    parser.add_argument(
        "input",
        help="LiMBS .txt file",
    )
    parser.add_argument(
        "--parser",
        default=None,
        help="Path to LiMBS parser (auto-detected if omitted)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output JSON path (default: <input>_linked.json)",
    )
    parser.add_argument(
        "--format",
        choices=["json", "tsv", "csv", "all"],
        default="json",
        help="Output format: json, tsv, csv, or all (default: json)",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print summary only; skip per-lipid details",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Exit with error if any requested output file already exists",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    input_path = Path(args.input).expanduser()
    if not input_path.is_file():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        raise SystemExit(1)

    try:
        parser_path = find_parser(str(input_path), args.parser)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)

    output_paths = default_output_paths(str(input_path))
    if args.out:
        output_paths["json"] = Path(args.out).expanduser()

    requested_paths: list[Path] = []
    if args.format in ("json", "all"):
        requested_paths.append(output_paths["json"])
    if args.format in ("tsv", "all"):
        requested_paths.append(output_paths["tsv"])
    if args.format in ("csv", "all"):
        requested_paths.append(output_paths["csv"])

    if args.no_overwrite:
        check_no_overwrite(requested_paths)

    print(f"Parser : {parser_path}")
    print(f"Input  : {input_path}")

    try:
        limbs_json = parse_limbs_to_json(str(input_path), parser_path)
        enriched = enrich_with_external_ids(limbs_json)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)

    print_linkage_report(enriched, verbose=not args.summary)

    print("\nOutput files:")
    if args.format in ("json", "all"):
        write_json(enriched, output_paths["json"])
    if args.format in ("tsv", "all"):
        write_tsv(enriched, output_paths["tsv"])
    if args.format in ("csv", "all"):
        write_csv(enriched, output_paths["csv"])

    print()


if __name__ == "__main__":
    main()
