#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Optional



# Data classes

@dataclass
class LipidSpec:
    """Lipid identity and physicochemical properties from the header block."""
    name: str
    tails: str = ""
    headgroup: str = ""
    charge: float = 0.0
    apl: Optional[float] = None


@dataclass
class LeafletComposition:
    """Lipid counts for a planar bilayer."""
    upper: dict = field(default_factory=dict)
    lower: dict = field(default_factory=dict)


@dataclass
class CGLayer:
    """One level of a CG mapping."""
    index: int
    content: str


@dataclass
class LipidCGBlock:
    """CG mapping block for one lipid species."""
    name: str
    layers: list = field(default_factory=list)
    fragments: dict = field(default_factory=dict)


@dataclass
class LiMBSSystem:
    # System block
    lipid_specs: list
    leaflets: LeafletComposition
    geometry: str
    box: list
    solvent: list
    salt_conc: float
    salt_ion: str
    rand: float
    solr: float
    protein_name: Optional[str] = None
    protein_file: Optional[str] = None
    temperature: Optional[float] = None
    ph: Optional[float] = None

    # Chemical block
    lipid_cg_blocks: list = field(default_factory=list)
    is_cg: bool = False

    # Validation
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)



# File reader

def read_limbs_file(path: str) -> str:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: '{path}'")

    lines = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            stripped = raw_line.strip()
            if stripped and not stripped.startswith("#"):
                lines.append(stripped)

    if not lines:
        raise ValueError(f"No valid content found in: '{path}'")

    joined = " ".join(lines)
    joined = re.sub(r"[ \t]+", " ", joined)
    return joined


# Parser

class LiMBSParser:
    SUPPORTED_GEOMETRIES = {"planar"}

    def parse(self, raw: str) -> LiMBSSystem:
        tokens = self._split_top_level(raw)
        system_tokens = []
        chemical_tokens = []

        for tok in tokens:
            tok = tok.strip()
            if not tok:
                continue
            if re.match(r"[A-Za-z0-9_]+\s*=\s*CG\s*\|", tok):
                chemical_tokens.append(tok)
            else:
                system_tokens.append(tok)

        system_str = " || ".join(system_tokens)
        lipid_specs = self._parse_lipid_specs(system_str)
        leaflets = self._parse_leaflets(system_str)
        geometry = self._parse_geometry(system_str)
        box = self._parse_box(system_str)
        solvent = self._parse_solvent(system_str)
        salt_conc, salt_ion = self._parse_salt(system_str)
        rand = self._parse_float_key(system_str, r"rand\s*:\s*([0-9.]+)", default=0.1)
        solr = self._parse_float_key(system_str, r"solr\s*:\s*([0-9.]+)", default=0.5)
        temperature = self._parse_float_key(system_str, r"(?:temp|temperature)\s*:\s*([0-9.]+)", default=None)
        ph = self._parse_float_key(system_str, r"pH\s*:\s*([0-9.]+)", default=None)
        protein_name, protein_file = self._parse_protein(system_str)

        lipid_cg_blocks = [self._parse_cg_block(t) for t in chemical_tokens]
        is_cg = len(lipid_cg_blocks) > 0

        system = LiMBSSystem(
            lipid_specs=lipid_specs,
            leaflets=leaflets,
            geometry=geometry,
            box=box,
            solvent=solvent,
            salt_conc=salt_conc,
            salt_ion=salt_ion,
            rand=rand,
            solr=solr,
            protein_name=protein_name,
            protein_file=protein_file,
            temperature=temperature,
            ph=ph,
            lipid_cg_blocks=lipid_cg_blocks,
            is_cg=is_cg,
        )
        self._validate(system)
        return system

    def _split_top_level(self, s: str) -> list:
        """Split on '||' not inside [] or {}."""
        tokens, buf, depth, i = [], [], 0, 0
        while i < len(s):
            if s[i] in "[{":
                depth += 1
                buf.append(s[i])
                i += 1
            elif s[i] in "]}":
                depth -= 1
                buf.append(s[i])
                i += 1
            elif s[i:i + 2] == "||" and depth == 0:
                tokens.append("".join(buf).strip())
                buf = []
                i += 2
            else:
                buf.append(s[i])
                i += 1
        if buf:
            tokens.append("".join(buf).strip())
        return [t for t in tokens if t]

    def _parse_lipid_specs(self, block: str) -> list:
        header = block.split("leaflets")[0] if "leaflets" in block else block
        specs, seen = [], set()

        pattern = re.compile(
            r"\b([A-Z][A-Z0-9]+)"                        # lipid name
            r":([\w:\-.]+)"                              # tail spec
            r"(?:[,\s]+([\w]+)-head)?"                   # optional headgroup
            r"(?:[,\s]+charge\s*([+\-−]?\s*[0-9.]+))?"  # optional charge
        )

        for m in pattern.finditer(header):
            name = m.group(1)
            if name in seen:
                continue
            seen.add(name)
            tails = m.group(2) or ""
            headgroup = m.group(3) or ""
            charge_raw = (m.group(4) or "0.0").replace("−", "-").replace(" ", "")
            try:
                charge = float(charge_raw)
            except ValueError:
                charge = 0.0
            specs.append(LipidSpec(name=name, tails=tails, headgroup=headgroup, charge=charge))
        return specs

    def _parse_leaflets(self, block: str) -> LeafletComposition:
        m = re.search(r"leaflets\s*-u\s*\{([^}]+)\}\s*,?\s*-l\s*\{([^}]+)\}", block, re.IGNORECASE)
        if not m:
            raise ValueError(
                "Leaflet specification '-u{...}, -l{...}' not found. "
                "Check that both leaflets use a leading dash: -u{...}, -l{...}"
            )
        return LeafletComposition(
            upper=self._parse_leaflet_side(m.group(1)),
            lower=self._parse_leaflet_side(m.group(2)),
        )

    def _parse_leaflet_side(self, s: str) -> dict:
        result = {}
        for item in s.split(","):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                raise ValueError(f"Bad leaflet item: '{item}'. Expected format LIPID:count")
            name, val = item.split(":", 1)
            name, val = name.strip(), val.strip()
            result[name] = int(val) if val.isdigit() else float(val)
        return result

    def _parse_geometry(self, block: str) -> str:
        m = re.search(r"type\s*:\s*(\w+)", block, re.IGNORECASE)
        geom = m.group(1).lower() if m else "planar"
        if geom not in self.SUPPORTED_GEOMETRIES:
            raise ValueError(
                f"Unsupported geometry '{geom}'."
            )
        return geom

    def _parse_box(self, block: str) -> list:
        m = re.search(r"box\s*:\s*\[([^\]]+)\]", block, re.IGNORECASE)
        if not m:
            return [10.0, 10.0, 8.0]
        box = [float(x.strip()) for x in m.group(1).split(",")]
        if len(box) != 3:
            raise ValueError("box must contain exactly three values, for example box:[10,10,8]")
        return box

    def _parse_solvent(self, block: str) -> list:
        m = re.search(r"sol\s*:\s*([^\s|,]+)", block, re.IGNORECASE)
        return m.group(1).strip().split("-") if m else ["W"]

    def _parse_salt(self, block: str):
        m = re.search(r"salt\s*:\s*([0-9.]+)\s*([A-Za-z]+)?", block, re.IGNORECASE)
        return (float(m.group(1)), (m.group(2) or "NaCl").strip()) if m else (0.15, "NaCl")

    def _parse_float_key(self, block: str, pattern: str, default=None):
        m = re.search(pattern, block, re.IGNORECASE)
        return float(m.group(1)) if m else default

    def _parse_protein(self, block: str):
        m = re.search(r"Protein\s*:\s*([^\s|;,]+)", block, re.IGNORECASE)
        if not m:
            return None, None
        name = m.group(1).strip()
        if os.path.isfile(name):
            return os.path.splitext(os.path.basename(name))[0], name
        for ext in (".gro", ".pdb"):
            if os.path.isfile(name + ext):
                return name, name + ext
        return name, None

    def _parse_cg_block(self, token: str) -> LipidCGBlock:
        """Parse one LIPIDNAME=CG|[layer1|layer2] token."""
        name_m = re.match(r"([A-Za-z0-9_]+)\s*=\s*CG\s*\|", token)
        name = name_m.group(1) if name_m else "UNKNOWN"

        inner_m = re.search(r"\[\s*(.*)\s*\]", token, re.DOTALL)
        if not inner_m:
            return LipidCGBlock(name=name)

        layer_strings = self._split_layers(inner_m.group(1).strip())
        layers = [CGLayer(index=i + 1, content=s.strip()) for i, s in enumerate(layer_strings) if s.strip()]
        fragments = self._extract_fragments(token)
        return LipidCGBlock(name=name, layers=layers, fragments=fragments)

    def _split_layers(self, s: str) -> list:
        """Split on '|' not inside [] or {}."""
        parts, buf, depth = [], [], 0
        for c in s:
            if c in "[{":
                depth += 1
                buf.append(c)
            elif c in "]}":
                depth -= 1
                buf.append(c)
            elif c == "|" and depth == 0:
                parts.append("".join(buf))
                buf = []
            else:
                buf.append(c)
        if buf:
            parts.append("".join(buf))
        return parts

    def _extract_fragments(self, block: str) -> dict:
        return {m.group(1): m.group(2).strip() for m in re.finditer(r"\{#(\w+)\s*=\s*([^}]+)\}", block)}

    def _validate(self, system: LiMBSSystem) -> None:
        if not system.leaflets.upper:
            system.errors.append("Upper leaflet is empty.")
        if not system.leaflets.lower:
            system.errors.append("Lower leaflet is empty.")
        if any(v <= 0 for v in system.box):
            system.errors.append("All box dimensions must be positive.")
        if system.protein_name and not system.protein_file:
            system.warnings.append(
                f"Protein '{system.protein_name}' was requested, but no matching .gro or .pdb file was found."
            )
        leaflet_lipids = set(system.leaflets.upper) | set(system.leaflets.lower)
        spec_lipids = {ls.name for ls in system.lipid_specs}
        missing_specs = sorted(leaflet_lipids - spec_lipids)
        if missing_specs:
            system.warnings.append(
                "These leaflet lipids do not have header lipid specs: " + ", ".join(missing_specs)
            )


# Builder

class LiMBSBuilder:
    INSANE_SCRIPT = "insane_peptoid_py3.py"

    def build(self, s: LiMBSSystem, dry_run=False, out_gro="membrane.gro", out_top="membrane.top") -> bool:
        if s.errors:
            print("\nCannot build. Resolve validation errors first:")
            for err in s.errors:
                print(f"  - {err}")
            return False

        if s.geometry != "planar":
            print("\nCannot build.")
            return False

        return self._build_insane(s, dry_run, out_gro, out_top)

    def _build_insane(self, s: LiMBSSystem, dry_run: bool, out_gro: str, out_top: str) -> bool:
        args = []
        if s.protein_file:
            args += ["-f", s.protein_file]
        for lip, count in s.leaflets.upper.items():
            args += ["-u", f"{lip}:{count}"]
        for lip, count in s.leaflets.lower.items():
            args += ["-l", f"{lip}:{count}"]
        for sv in s.solvent:
            args += ["-sol", sv]
        args += ["-x", str(s.box[0]), "-y", str(s.box[1]), "-z", str(s.box[2])]
        args += ["-salt", str(s.salt_conc), "-rand", str(s.rand), "-solr", str(s.solr)]
        args += ["-o", out_gro, "-p", out_top]
        return self._run(["python3", self.INSANE_SCRIPT] + args, dry_run)

    def _run(self, cmd: list, dry_run: bool) -> bool:
        print("\nCommand:")
        print("  " + " ".join(cmd))
        if dry_run:
            print("  [dry-run mode: command not executed]")
            return True
        retcode = subprocess.call(cmd)
        ok = retcode == 0
        print("\nBuild completed." if ok else f"\nBuild failed with return code {retcode}.")
        return ok


# Summary

def print_summary(result: LiMBSSystem) -> None:
    print("\nLiMBS Parsed System")
    print(f"  Geometry     : {result.geometry}")
    print(f"  Box (nm)     : {result.box}")
    print(f"  Solvent      : {result.solvent}")
    print(f"  Salt         : {result.salt_conc} M {result.salt_ion}")
    print(f"  Temp / pH    : {result.temperature} / {result.ph}")
    print(f"  rand / solr  : {result.rand} / {result.solr}")
    print(f"  Resolution   : {'Coarse-grained' if result.is_cg else 'Atomistic'}")
    print(f"  Lipid specs  : {[ls.name for ls in result.lipid_specs]}")
    print(f"  Upper leaflet: {result.leaflets.upper}")
    print(f"  Lower leaflet: {result.leaflets.lower}")

    if result.protein_name:
        print(f"  Protein      : {result.protein_name} ({result.protein_file})")

    if result.is_cg:
        print(f"  CG blocks    : {[b.name for b in result.lipid_cg_blocks]}")
        for b in result.lipid_cg_blocks:
            print(f"    {b.name}: {len(b.layers)} layer(s), {len(b.fragments)} fragment(s)")

    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print(f"  - {warning}")

    if result.errors:
        print("\nErrors:")
        for error in result.errors:
            print(f"  - {error}")
    print()


def system_to_json(result: LiMBSSystem) -> dict:
    return {
        "geometry": result.geometry,
        "box_nm": result.box,
        "solvent": result.solvent,
        "salt": {"conc": result.salt_conc, "ion": result.salt_ion},
        "temperature": result.temperature,
        "ph": result.ph,
        "rand": result.rand,
        "solr": result.solr,
        "is_cg": result.is_cg,
        "upper_leaflet": result.leaflets.upper,
        "lower_leaflet": result.leaflets.lower,
        "lipid_specs": [
            {"name": ls.name, "tails": ls.tails, "headgroup": ls.headgroup, "charge": ls.charge}
            for ls in result.lipid_specs
        ],
        "cg_blocks": [
            {"name": b.name, "layers": len(b.layers), "fragments": list(b.fragments.keys())}
            for b in result.lipid_cg_blocks
        ],
        "protein": result.protein_name,
        "protein_file": result.protein_file,
        "warnings": result.warnings,
        "errors": result.errors,
    }



def main() -> None:
    ap = argparse.ArgumentParser(
        description="LiMBS: read, parse and build planar lipid membranes from LiMBS notation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples

  python3 LiMBS.py bilayer.txt
        """,
    )
    ap.add_argument("input", help="LiMBS notation file. Lines starting with # are ignored.")
    ap.add_argument("--dry-run", action="store_true", help="Print the builder command without executing it")
    ap.add_argument("--validate-only", action="store_true", help="Parse and validate only; do not build")
    ap.add_argument("--out-gro", default="membrane.gro", help="Output coordinate file name")
    ap.add_argument("--out-top", default="membrane.top", help="Output topology file name")
    ap.add_argument("--json", action="store_true", help="Print the parsed system as JSON")
    args = ap.parse_args()

    print(f"\n[1/3] Reading LiMBS file: {args.input}")
    try:
        raw = read_limbs_file(args.input)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}")
        sys.exit(1)

    print(f"      {len(raw)} characters read across {raw.count('||') + 1} top-level tokens.")

    print("\n[2/3] Parsing LiMBS notation...")
    parser = LiMBSParser()
    try:
        result = parser.parse(raw)
    except ValueError as exc:
        print(f"Parse error: {exc}")
        sys.exit(1)

    print_summary(result)

    if args.json:
        print(json.dumps(system_to_json(result), indent=2))

    if args.validate_only:
        if result.errors:
            sys.exit(1)
        print("Validation completed.")
        return

    print("\n[3/3] Building planar membrane...")
    builder = LiMBSBuilder()
    success = builder.build(result, dry_run=args.dry_run, out_gro=args.out_gro, out_top=args.out_top)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

