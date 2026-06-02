#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class LeafletComposition:
    """Lipid ratios for vesicle outer and inner leaflets."""
    upper: Dict[str, float] = field(default_factory=dict)  # Domain 0 / outer
    lower: Dict[str, float] = field(default_factory=dict)  # Domain 1 / inner


@dataclass
class LiMBSSystem:
    """Parsed LiMBS vesicle system."""
    leaflets: LeafletComposition
    geometry: str
    box: List[float]
    solvent: List[str]
    salt_conc: float
    salt_ion: str
    temperature: Optional[float]
    ph: Optional[float]
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


def read_limbs_file(path: str) -> str:
    """Read LiMBS notation from file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    lines: List[str] = []
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line and not line.startswith("#"):
                lines.append(line)

    if not lines:
        raise ValueError(f"No valid LiMBS content found in: {path}")

    joined = " ".join(lines)
    return re.sub(r"[ \t]+", " ", joined)


class LiMBSParser:
    """Parser for the system block of a LiMBS vesicle notation."""

    def parse(self, raw: str) -> LiMBSSystem:
        geometry = self._parse_geometry(raw)
        leaflets = self._parse_leaflets(raw)
        box = self._parse_box(raw)
        solvent = self._parse_solvent(raw)
        salt_conc, salt_ion = self._parse_salt(raw)
        temperature = self._parse_float_key(raw, r"temp:([0-9.]+)")
        ph = self._parse_float_key(raw, r"pH:([0-9.]+)")

        system = LiMBSSystem(
            leaflets=leaflets,
            geometry=geometry,
            box=box,
            solvent=solvent,
            salt_conc=salt_conc,
            salt_ion=salt_ion,
            temperature=temperature,
            ph=ph,
        )
        self._validate(system)
        return system

    def _parse_geometry(self, block: str) -> str:
        match = re.search(r"type\s*:\s*(\w+)", block, re.IGNORECASE)
        geometry = match.group(1).lower() if match else "vesicle"

        if geometry != "vesicle":
            raise ValueError(
                f"This script is TS2CG-vesicle only. Use type:vesicle, not type:{geometry}"
            )
        return geometry

    def _parse_leaflets(self, block: str) -> LeafletComposition:
        match = re.search(
            r"leaflets\s*-u\{([^}]+)\}\s*,?\s*-l\{([^}]+)\}",
            block,
            re.IGNORECASE,
        )
        if not match:
            raise ValueError(
                "Leaflet block not found. Expected format: "
                "leaflets -u{DPPC:0.64,CHOL:0.36} -l{DPPC:0.64,CHOL:0.36}"
            )

        upper = self._parse_leaflet_side(match.group(1))
        lower = self._parse_leaflet_side(match.group(2))
        return LeafletComposition(upper=upper, lower=lower)

    def _parse_leaflet_side(self, text: str) -> Dict[str, float]:
        result: Dict[str, float] = {}

        for item in text.split(","):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                raise ValueError(f"Invalid leaflet item '{item}'. Expected LIPID:value")

            name, value = item.split(":", 1)
            name = name.strip()
            value = value.strip()

            if not name:
                raise ValueError(f"Missing lipid name in leaflet item '{item}'")

            try:
                result[name] = float(value)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid value for lipid '{name}': {value}. Use a number."
                ) from exc

        if not result:
            raise ValueError("Empty leaflet composition found.")

        return result

    def _parse_box(self, block: str) -> List[float]:
        match = re.search(r"box:\[([^\]]+)\]", block, re.IGNORECASE)
        if not match:
            return [20.0, 20.0, 20.0]

        values = [float(x.strip()) for x in match.group(1).split(",")]
        if len(values) != 3:
            raise ValueError("Box must contain three values, for example box:[20,20,20]")
        return values

    def _parse_solvent(self, block: str) -> List[str]:
        match = re.search(r"sol:([^\s|,\|]+)", block, re.IGNORECASE)
        return match.group(1).strip().split("-") if match else ["W"]

    def _parse_salt(self, block: str) -> Tuple[float, str]:
        match = re.search(r"salt:([0-9.]+)\s*([A-Za-z0-9+\-]+)?", block, re.IGNORECASE)
        if not match:
            return 0.15, "NaCl"
        return float(match.group(1)), (match.group(2) or "NaCl").strip()

    def _parse_float_key(self, block: str, pattern: str) -> Optional[float]:
        match = re.search(pattern, block, re.IGNORECASE)
        return float(match.group(1)) if match else None

    def _validate(self, system: LiMBSSystem) -> None:
        if any(x <= 0 for x in system.box):
            system.errors.append(f"Box dimensions must be positive: {system.box}")

        for leaflet_name, leaflet in (
            ("upper/domain 0", system.leaflets.upper),
            ("lower/domain 1", system.leaflets.lower),
        ):
            total = sum(leaflet.values())
            if total <= 0:
                system.errors.append(f"{leaflet_name} composition total must be positive.")
            if total > 10:
                system.warnings.append(
                    f"{leaflet_name} total is {total}. "
                    "For TS2CG vesicles, ratios are usually preferred over absolute counts."
                )

        if set(system.leaflets.upper) != set(system.leaflets.lower):
            system.warnings.append(
                "Upper and lower leaflets contain different lipid names. "
            )


class TS2CGBuilder:
    """Write input.str and run TS2CG PLM + PCG."""

    TS2CG_COMMAND = "TS2CG"

    def build(
        self,
        system: LiMBSSystem,
        dry_run: bool = False,
        tsi_file: str = "sphere.tsi",
        martini_lib: str = "Martini3.LIB",
        str_file: str = "input.str",
        out_name: str = "vesicle",
        bilayer_thickness: float = 3.0,
        rescale: Tuple[float, float, float] = (4.8, 4.8, 4.8),
        mashno: int = 4,
        bond_length: float = 0.2,
    ) -> bool:
        if system.errors:
            print("\nCannot build as errors were found:")
            for err in system.errors:
                print(f"  - {err}")
            return False

        self.write_input_str(system, str_file)

        if not os.path.isfile(tsi_file):
            print(f"\nError: TS2CG geometry file not found: {tsi_file}")
            return False

        if not os.path.isfile(martini_lib):
            print(f"\nError: Martini lipid library not found: {martini_lib}")
            return False

        plm_cmd = [
            self.TS2CG_COMMAND,
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
            self.TS2CG_COMMAND,
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

        print("\n[TS2CG input.str]")
        print(f"  Wrote: {str_file}")

        print("\n[Step 1] TS2CG PLM")
        if not self._run(plm_cmd, dry_run):
            return False

        print("\n[Step 2] TS2CG PCG")
        if not self._run(pcg_cmd, dry_run):
            return False

        print(f"\nFinished. Expected output prefix: {out_name}")
        return True

    def write_input_str(self, system: LiMBSSystem, path: str) -> None:
        """
        Write TS2CG input.str for a single-domain vesicle.

        Both monolayers are written under Domain 0 because the default
        sphere.tsi/point file usually contains only Domain 0 points.

        Format:
            [Lipids List]
            Domain 0
            DPPC 1 1 0.64   # upper/outer monolayer
            DPPC 2 1 0.64   # lower/inner monolayer
            End
            ~
        """
        with open(path, "w") as fh:
            fh.write("[Lipids List]\n")
            fh.write("Domain 0\n")

            # Upper/outer monolayer: leaflet id 1
            for lipid, value in system.leaflets.upper.items():
                fh.write(f"{lipid} 1 1 {value:g}\n")

            # Lower/inner monolayer: leaflet id 2
            for lipid, value in system.leaflets.lower.items():
                fh.write(f"{lipid} 2 1 {value:g}\n")

            fh.write("End\n")
           # fh.write("~\n")

    def _run(self, command: List[str], dry_run: bool) -> bool:
        print("  " + " ".join(command))

        if dry_run:
            print("  dry-run: command not executed")
            return True

        try:
            completed = subprocess.run(command, check=False)
        except FileNotFoundError:
            print(f"\nError: command not found: {command[0]}")
            return False

        if completed.returncode != 0:
            print(f"\nCommand failed  {completed.returncode}")
            return False

        return True


def print_summary(system: LiMBSSystem) -> None:
    print("\nParsed LiMBS vesicle system")
    print(f"  Geometry      : {system.geometry}")
    print(f"  Box           : {system.box} nm")
    print(f"  Solvent       : {system.solvent}")
    print(f"  Salt          : {system.salt_conc} M {system.salt_ion}")
    print(f"  Temperature   : {system.temperature}")
    print(f"  pH            : {system.ph}")
    print(f"  Domain 0 / upper / outer: {system.leaflets.upper}")
    print(f"  Domain 1 / lower / inner: {system.leaflets.lower}")

    if system.warnings:
        print("\nWarnings:")
        for warning in system.warnings:
            print(f"  - {warning}")

    if system.errors:
        print("\nErrors:")
        for error in system.errors:
            print(f"  - {error}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a TS2CG vesicle input from LiMBS notation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("input", help="LiMBS notation file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running TS2CG")
    parser.add_argument("--validate-only", action="store_true", help="Parse and write nothing else")
    parser.add_argument("--write-str-only", action="store_true", help="Only write input.str, do not run TS2CG")

    parser.add_argument("--tsi", default="sphere.tsi", help="TS2CG .tsi geometry file")
    parser.add_argument("--lib", default="Martini3.LIB", help="TS2CG Martini lipid library file")
    parser.add_argument("--str-file", default="input.str", help="Output TS2CG composition file")
    parser.add_argument("--out-name", default="vesicle", help="TS2CG output prefix")

    parser.add_argument("--bilayer-thickness", type=float, default=3.0, help="TS2CG bilayer thickness")
    parser.add_argument("--rescale", type=float, nargs=3, default=(4.8, 4.8, 4.8), help="TS2CG rescale factors")
    parser.add_argument("--mashno", type=int, default=4, help="TS2CG Mashno value")
    parser.add_argument("--bond-length", type=float, default=0.2, help="TS2CG PCG bond length")

    args = parser.parse_args()

    try:
        raw = read_limbs_file(args.input)
        system = LiMBSParser().parse(raw)
    except (FileNotFoundError, ValueError) as exc:
        print(f"\nError: {exc}")
        sys.exit(1)

    print_summary(system)

    if system.errors:
        sys.exit(1)

    if args.validate_only:
        print("\nValidation passed.")
        sys.exit(0)

    builder = TS2CGBuilder()

    if args.write_str_only:
        builder.write_input_str(system, args.str_file)
        print(f"\nWrote {args.str_file}")
        sys.exit(0)

    success = builder.build(
        system=system,
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
