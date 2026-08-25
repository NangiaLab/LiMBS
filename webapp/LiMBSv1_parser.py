#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LiMBS1.0 parser, validator, JSON exporter, and canonicalizer.

This version supports the current LiMBS line notation:

System block:
    LiMBS1.0 ||
    rsln: AA | CG | SCG ||
    lsc: 1 ||
    uls: 1 ||
    lls: 1 ||
    {DPPC:tails=16:0-16:0,head=PC,charge=0.0} ||
    leaflet:count{-u{DPPC:56}-l{DPPC:56}} ||
    type:planar ||
    box:[10.0,10.0,8.0]nm ||
    sol:W ||
    salt:0.15M XY

Separator:
    |++|   separates the system block from the species-definition block

Species block:
    DPPC=AA|[...]       one AA species
    |+|
    DPPC=CG|[...]       one CG species
    |+|
    DPPC=SCG|[...]      one SCG species

Resolution logic:
    AA  : usually 1 layer, heavy-atom notation allowed
    CG  : usually 2 layers, atomistic/fragment representation -> CG bead mapping
    SCG : usually 3 layers, atomistic/fragment representation -> CG mapping -> SCG semantic grouping
"""

import argparse
import json
import re
import sys
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


# Data models

class LiMBSError(Exception):
    pass


@dataclass
class Message:
    code: str
    message: str


@dataclass
class SpeciesBlock:
    name: str
    resolution: str
    body: str
    layer_count: int = 0
    layers: List[Dict[str, Any]] = field(default_factory=list)
    fragments_defined: List[str] = field(default_factory=list)
    fragments_used: List[str] = field(default_factory=list)
    beads: List[str] = field(default_factory=list)
    semantic_groups: List[str] = field(default_factory=list)
    semantic_group_definitions: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class LiMBSSystem:
    version: str = "1.0"

    # New compact header fields
    resolution: str = ""
   # forcefield: str = ""
    lipid_species_count: Optional[int] = None
    upper_leaflet_species_count: Optional[int] = None
    lower_leaflet_species_count: Optional[int] = None

    # System fields
    lipids: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    leaflets: Dict[str, Any] = field(default_factory=dict)
    membrane_type: str = ""
    box: Dict[str, Any] = field(default_factory=dict)
    solvent: str = ""
    salt: Dict[str, Any] = field(default_factory=dict)

    # Species definitions
    species_blocks: Dict[str, SpeciesBlock] = field(default_factory=dict)

    errors: List[Message] = field(default_factory=list)
    warnings: List[Message] = field(default_factory=list)


# Parser

class LiMBSParser:
    VERSION_TOKEN = "LiMBS1.0"

    SUPPORTED_RESOLUTIONS = {"AA", "CG", "SCG"}
    SUPPORTED_TYPES = {"planar", "vesicle"}
    SUPPORTED_LEAFLET_MODES = {"count", "ratio"}

    def parse(self, text: str) -> LiMBSSystem:
        text = self._clean_text(text)

        if not text:
            raise LiMBSError("LMB-E000: Empty LiMBS input.")

        if text.count("|++|") != 1:
            raise LiMBSError(
                "LMB-E002: LiMBS string must contain exactly one |++| separator "
                "between the system block and the species-definition block."
            )

        system_text, species_text = text.split("|++|", 1)

        system = self._parse_system_block(system_text)
        system.species_blocks = self._parse_species_blocks(species_text)

        self.validate(system)
        return system

    def _clean_text(self, text: str) -> str:
        """
        Removes comments and blank lines, normalizes whitespace, and normalizes minus signs.
        """
        lines = []

        for raw in text.splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            lines.append(line)

        cleaned = " ".join(lines)
        cleaned = re.sub(r"\s+", " ", cleaned)
        cleaned = cleaned.replace("−", "-")
        return cleaned.strip()

    def _parse_system_block(self, system_text: str) -> LiMBSSystem:
        parts = self._split_top_level_double_pipe(system_text)

        if not parts:
            raise LiMBSError("LMB-E000: Empty system block.")

        if parts[0].strip() != self.VERSION_TOKEN:
            raise LiMBSError("LMB-E001: LiMBS string must begin with LiMBS1.0.")

        system = LiMBSSystem(version="1.0")
        seen_fields = set()

        for part in parts[1:]:
            part = part.strip()

            if not part:
                continue

            if part.startswith("rsln:"):
                self._check_duplicate(seen_fields, "rsln")
                system.resolution = self._parse_resolution(part)


            elif part.startswith("lsc:"):
                self._check_duplicate(seen_fields, "lsc")
                system.lipid_species_count = self._parse_positive_int_field(part, "lsc")

            elif part.startswith("uls:"):
                self._check_duplicate(seen_fields, "uls")
                system.upper_leaflet_species_count = self._parse_positive_int_field(part, "uls")

            elif part.startswith("lls:"):
                self._check_duplicate(seen_fields, "lls")
                system.lower_leaflet_species_count = self._parse_positive_int_field(part, "lls")

            elif self._is_lipid_metadata_block(part):
                self._check_duplicate(seen_fields, "lipid_metadata")
                system.lipids = self._parse_lipids(part)

            elif part.startswith("leaflet:") or part.startswith("leaflets:"):
                self._check_duplicate(seen_fields, "leaflet")
                system.leaflets = self._parse_leaflets(part)

            elif part.startswith("type:"):
                self._check_duplicate(seen_fields, "type")
                system.membrane_type = self._parse_type(part)

            elif part.startswith("box:"):
                self._check_duplicate(seen_fields, "box")
                system.box = self._parse_box(part)

            elif part.startswith("sol:"):
                self._check_duplicate(seen_fields, "sol")
                system.solvent = self._parse_solvent(part)

            elif part.startswith("salt:"):
                self._check_duplicate(seen_fields, "salt")
                system.salt = self._parse_salt(part)

            else:
                system.warnings.append(
                    Message("LMB-W001", f"Unknown system field ignored: {part}")
                )

        return system

    def _check_duplicate(self, seen_fields: set, field: str) -> None:
        if field in seen_fields:
            raise LiMBSError(f"LMB-E003: Duplicate system field: {field}")
        seen_fields.add(field)

    def _parse_resolution(self, part: str) -> str:
        m = re.fullmatch(r"rsln:\s*(AA|CG|SCG)", part, flags=re.IGNORECASE)

        if not m:
            raise LiMBSError("LMB-E004: Invalid resolution. Expected rsln: AA, rsln: CG, or rsln: SCG.")

        return m.group(1).upper()

    def _parse_positive_int_field(self, part: str, key: str) -> int:
        m = re.fullmatch(rf"{re.escape(key)}:\s*(\d+)", part)

        if not m:
            raise LiMBSError(f"LMB-E005: Invalid {key} field. Expected {key}: 1")

        value = int(m.group(1))

        if value <= 0:
            raise LiMBSError(f"LMB-E006: {key} must be positive.")

        return value

    def _is_lipid_metadata_block(self, part: str) -> bool:
        return (
            part.startswith("{")
            and part.endswith("}")
            and "tails=" in part
            and "head=" in part
            and "charge=" in part
        )

    def _parse_lipids(self, part: str) -> Dict[str, Dict[str, Any]]:
        m = re.fullmatch(r"\{(.+)\}", part)

        if not m:
            raise LiMBSError("LMB-E010: Invalid lipid metadata block.")

        content = m.group(1).strip()
        lipids = {}

        pattern = re.compile(
            r"([A-Za-z][A-Za-z0-9_]*):"
            r"tails=([^,{};]+),"
            r"head=([A-Za-z0-9_]+),"
            r"charge=([+-]?\d+(?:\.\d+)?)"
        )

        matches = list(pattern.finditer(content))

        if not matches:
            raise LiMBSError("LMB-E011: No valid lipid entries found in lipid metadata block.")

        leftover = content
        for match in matches:
            leftover = leftover.replace(match.group(0), "", 1)

        leftover = leftover.replace(",", "").replace(";", "").strip()

        if leftover:
            raise LiMBSError(f"LMB-E011: Unparsed lipid text in lipid metadata block: {leftover}")

        for match in matches:
            name, tails, head, charge = match.groups()

            if name in lipids:
                raise LiMBSError(f"LMB-E012: Duplicate lipid entry: {name}")

            lipids[name] = {
                "tails": tails,
                "head": head,
                "charge": float(charge),
            }

        return lipids

    def _parse_leaflets(self, part: str) -> Dict[str, Any]:
        m = re.fullmatch(
            r"leaflets?:(count|ratio)\{-u\{([^}]*)\}-l\{([^}]*)\}\}",
            part,
            flags=re.IGNORECASE,
        )

        if not m:
            raise LiMBSError(
                "LMB-E020: Invalid leaflet block. Expected: "
                "leaflet:count{-u{DPPC:56}-l{DPPC:56}} "
                "or leaflet:ratio{-u{DPPC:1.0}-l{DPPC:1.0}}"
            )

        mode, upper_text, lower_text = m.groups()
        mode = mode.lower()

        return {
            "mode": mode,
            "upper": self._parse_composition(upper_text),
            "lower": self._parse_composition(lower_text),
        }

    def _parse_composition(self, text: str) -> Dict[str, Any]:
        comp = {}
        text = text.strip()

        if not text:
            return comp

        for item in text.split(","):
            item = item.strip()

            if ":" not in item:
                raise LiMBSError(f"LMB-E021: Invalid leaflet item: {item}")

            name, value = item.split(":", 1)
            name = name.strip()
            value = value.strip()

            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
                raise LiMBSError(f"LMB-E022: Invalid lipid name in leaflet: {name}")

            try:
                parsed_value = float(value) if "." in value else int(value)
            except ValueError:
                raise LiMBSError(f"LMB-E023: Invalid leaflet value for {name}: {value}")

            comp[name] = parsed_value

        return comp

    def _parse_type(self, part: str) -> str:
        m = re.fullmatch(r"type:\s*([A-Za-z]+)", part)

        if not m:
            raise LiMBSError("LMB-E030: Invalid type block. Use type:planar or type:vesicle.")

        return m.group(1).lower()

    def _parse_box(self, part: str) -> Dict[str, Any]:
        m = re.fullmatch(
            r"box:\["
            r"([+-]?\d+(?:\.\d+)?),"
            r"([+-]?\d+(?:\.\d+)?),"
            r"([+-]?\d+(?:\.\d+)?)"
            r"\]([A-Za-z0-9]+)",
            part,
        )

        if not m:
            raise LiMBSError("LMB-E040: Invalid box block. Expected: box:[10.0,10.0,8.0]nm")

        x, y, z, unit = m.groups()

        return {
            "values": [float(x), float(y), float(z)],
            "unit": unit,
        }

    def _parse_solvent(self, part: str) -> str:
        m = re.fullmatch(r"sol:\s*([A-Za-z][A-Za-z0-9_+-]*)", part)

        if not m:
            raise LiMBSError("LMB-E050: Invalid solvent block. Example: sol:W, sol:WW, or sol:TIP3P")

        return m.group(1)

    def _parse_salt(self, part: str) -> Dict[str, Any]:
        """
        Supports both generic salt species and chemical salt formulas.

        Examples:
            salt:0.15M XY
            salt:0.15M NaCl
            salt:0.15M KCl
            salt:0.15M CaCl2
            salt:0.15M MgCl2
            salt:0.15M Na Cl

        The parser does not assume a fixed salt list. It parses chemical
        formulas into element/count pairs when possible.
        """
        m = re.fullmatch(
            r"salt:\s*(\d+(?:\.\d+)?)M\s+(.+)",
            part,
            flags=re.IGNORECASE,
        )

        if not m:
            raise LiMBSError("LMB-E060: Invalid salt block. Example: salt:0.15M NaCl")

        concentration, species_text = m.groups()
        species_text = species_text.strip()

        if not species_text:
            raise LiMBSError("LMB-E061: Salt species cannot be empty.")

        species_tokens = species_text.split()

        # If salt is written as one formula, parse it as a formula.
        # NaCl -> {Na:1, Cl:1}; CaCl2 -> {Ca:1, Cl:2}.
        formula = None
        if len(species_tokens) == 1:
            formula = self._parse_chemical_formula(species_tokens[0])

        # If salt is written as separate species, treat each token as one species.
        # Example: Na Cl -> {Na:1, Cl:1}
        separate_ions = None
        if len(species_tokens) > 1:
            separate_ions = {token: 1 for token in species_tokens}

        return {
            "concentration": float(concentration),
            "unit": "M",
            "species_text": species_text,
            "species": species_tokens,
            "formula": formula,
            "ions": formula if formula is not None else separate_ions,
        }

    def _parse_chemical_formula(self, formula: str) -> Optional[Dict[str, int]]:
        """
        Parse a simple inorganic formula such as NaCl, KCl, CaCl2, or MgCl2.

        Returns None for generic labels such as XY or NACL if they do not look
        like a real element formula.
        """
        # Common element symbols are enough for salt validation and prevent
        # interpreting arbitrary labels like XY as real chemistry.
        valid_elements = {
            "H", "Li", "Na", "K", "Rb", "Cs",
            "Mg", "Ca", "Sr", "Ba",
            "F", "Cl", "Br", "I",
            "Zn", "Fe", "Cu", "Mn", "Co", "Ni",
        }

        tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
        if not tokens:
            return None

        reconstructed = "".join(element + count for element, count in tokens)
        if reconstructed != formula:
            return None

        counts: Dict[str, int] = {}
        for element, count_text in tokens:
            if element not in valid_elements:
                return None

            count = int(count_text) if count_text else 1
            counts[element] = counts.get(element, 0) + count

        return counts

    def _parse_species_blocks(self, text: str) -> Dict[str, SpeciesBlock]:
        text = text.strip()

        if not text:
            raise LiMBSError("LMB-E070: Species-definition block is empty.")

        raw_blocks = self._split_top_level_species_separator(text)
        blocks: Dict[str, SpeciesBlock] = {}

        for raw in raw_blocks:
            block = raw.strip()

            if not block:
                continue

            m = re.fullmatch(
                r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(AA|CG|SCG)\s*\|\s*(.+)",
                block,
                flags=re.IGNORECASE | re.DOTALL,
            )

            if not m:
                raise LiMBSError(f"LMB-E071: Invalid species definition: {block[:120]}")

            name, resolution, body = m.groups()
            resolution = resolution.upper()
            body = body.strip()

            if name in blocks:
                raise LiMBSError(f"LMB-E072: Duplicate species definition: {name}")

            blocks[name] = self._analyze_species_body(name, resolution, body)

        if not blocks:
            raise LiMBSError("LMB-E073: No valid species definitions found.")

        return blocks

    def _analyze_species_body(self, name: str, resolution: str, body: str) -> SpeciesBlock:
        layers = self._split_layers(body)

        block = SpeciesBlock(
            name=name,
            resolution=resolution,
            body=body,
            layer_count=len(layers),
            layers=[{"index": i + 1, "body": layer.strip()} for i, layer in enumerate(layers)],
            fragments_defined=sorted(set(re.findall(r"\{#([A-Za-z0-9_]+)\s*=", body))),
            fragments_used=sorted(set(re.findall(r"\[#([A-Za-z0-9_]+)\]", body))),
            beads=sorted(set(re.findall(r"(?<!\^)\^([A-Za-z0-9_]+)", body))),
            semantic_groups=sorted(set(re.findall(r"\^\^([A-Za-z0-9_]+)", body))),
            semantic_group_definitions=self._parse_semantic_groups(body),
        )

        return block

    def _parse_semantic_groups(self, body: str) -> Dict[str, List[str]]:
        groups = {}

        for group_name, rhs in re.findall(r"\{\^\^([A-Za-z0-9_]+)\s*=\s*([^}]+)\}", body):
            members = []

            # Members such as [^NC3]
            members.extend([f"^{x}" for x in re.findall(r"(?<!\^)\^([A-Za-z0-9_]+)", rhs)])

            # Members such as [#TailA]
            members.extend([f"#{x}" for x in re.findall(r"\[#([A-Za-z0-9_]+)\]", rhs)])

            groups[group_name] = members

        return groups

    def _split_layers(self, body: str) -> List[str]:
        """
        Splits body layers by top-level single |.

        For:
            DPPC=SCG|[
            ...
            |
            ...
            |
            ...
            ]

        returns three layers.
        """
        inner = body.strip()

        if inner.startswith("[") and inner.endswith("]"):
            inner = inner[1:-1].strip()

        parts = []
        buffer = []
        depth = 0
        i = 0

        while i < len(inner):
            ch = inner[i]

            if ch in "[{(":
                depth += 1
                buffer.append(ch)
                i += 1

            elif ch in "]})":
                depth -= 1
                if depth < 0:
                    raise LiMBSError("LMB-E080: Unbalanced brackets in species definition.")
                buffer.append(ch)
                i += 1

            elif ch == "|" and depth == 0:
                parts.append("".join(buffer).strip())
                buffer = []
                i += 1

            else:
                buffer.append(ch)
                i += 1

        if depth != 0:
            raise LiMBSError("LMB-E081: Unbalanced brackets in species definition.")

        if buffer:
            parts.append("".join(buffer).strip())

        return [p for p in parts if p]

    def _split_top_level_double_pipe(self, text: str) -> List[str]:
        return self._split_top_level_token(text, "||")

    def _split_top_level_species_separator(self, text: str) -> List[str]:
        return self._split_top_level_token(text, "|+|")

    def _split_top_level_token(self, text: str, token: str) -> List[str]:
        parts = []
        buffer = []
        depth = 0
        i = 0

        while i < len(text):
            ch = text[i]

            if ch in "[{(":
                depth += 1
                buffer.append(ch)
                i += 1

            elif ch in "]})":
                depth -= 1
                buffer.append(ch)
                i += 1

            elif text[i:i + len(token)] == token and depth == 0:
                parts.append("".join(buffer).strip())
                buffer = []
                i += len(token)

            else:
                buffer.append(ch)
                i += 1

        if buffer:
            parts.append("".join(buffer).strip())

        return [p for p in parts if p]

    # Validation

    def validate(self, system: LiMBSSystem) -> None:
        self._validate_required_system_fields(system)
        self._validate_leaflets_and_lipids(system)
        self._validate_species_definitions(system)
        self._validate_header_counts(system)

    def _validate_required_system_fields(self, system: LiMBSSystem) -> None:
        if not system.resolution:
            system.errors.append(
                Message("LMB-E100", "Missing resolution field: rsln: AA, CG, or SCG.")
            )

        elif system.resolution not in self.SUPPORTED_RESOLUTIONS:
            system.errors.append(
                Message("LMB-E101", f"Unsupported resolution '{system.resolution}'.")
            )

        if not system.lipids:
            system.errors.append(Message("LMB-E102", "Missing lipid metadata block."))

        if not system.leaflets:
            system.errors.append(Message("LMB-E103", "Missing leaflet block."))

        if not system.membrane_type:
            system.errors.append(Message("LMB-E104", "Missing type block."))

        elif system.membrane_type not in self.SUPPORTED_TYPES:
            system.errors.append(
                Message(
                    "LMB-E105",
                    f"Unsupported membrane type '{system.membrane_type}'. "
                    f"Supported types: {sorted(self.SUPPORTED_TYPES)}",
                )
            )

        if not system.box:
            system.errors.append(Message("LMB-E106", "Missing box block."))

        else:
            if system.box.get("unit") != "nm":
                system.errors.append(
                    Message("LMB-E107", "Box unit must be nm in LiMBS1.0.")
                )

            if any(v <= 0 for v in system.box.get("values", [])):
                system.errors.append(
                    Message("LMB-E108", "All box dimensions must be positive.")
                )

        if not system.solvent:
            system.errors.append(Message("LMB-E109", "Missing solvent block."))

        else:
            if system.resolution == "CG" and system.solvent != "W":
                system.errors.append(
                    Message(
                        "LMB-E180",
                        f"CG systems require solvent 'W'. Found '{system.solvent}'."
                    )
                )

            elif system.resolution == "SCG" and system.solvent != "WW":
                system.errors.append(
                    Message(
                        "LMB-E181",
                        f"SCG systems require solvent 'WW'. Found '{system.solvent}'."
                    )
                )

            elif system.resolution == "AA":
                allowed_aa_solvents = {"TIP3P", "TIP4P", "SPC", "SPCE"}

                if system.solvent not in allowed_aa_solvents:
                    system.warnings.append(
                        Message(
                            "LMB-W180",
                            f"Unrecognized AA solvent '{system.solvent}'."
                        )
                    )

    def _validate_leaflets_and_lipids(self, system: LiMBSSystem) -> None:
        if not system.leaflets:
            return

        mode = system.leaflets.get("mode")
        upper = system.leaflets.get("upper", {})
        lower = system.leaflets.get("lower", {})

        if mode not in self.SUPPORTED_LEAFLET_MODES:
            system.errors.append(Message("LMB-E110", "Leaflet mode must be count or ratio."))

        if not upper:
            system.errors.append(Message("LMB-E111", "Upper leaflet is empty."))

        if not lower:
            system.errors.append(Message("LMB-E112", "Lower leaflet is empty."))

        if mode == "count":
            for leaflet_name, comp in [("upper", upper), ("lower", lower)]:
                for lipid, value in comp.items():
                    if not isinstance(value, int):
                        system.errors.append(
                            Message(
                                "LMB-E113",
                                f"Count mode requires integer values. "
                                f"{leaflet_name} leaflet has {lipid}:{value}",
                            )
                        )
                    elif value <= 0:
                        system.errors.append(
                            Message(
                                "LMB-E114",
                                f"Leaflet counts must be positive. "
                                f"{leaflet_name} leaflet has {lipid}:{value}",
                            )
                        )

        if mode == "ratio":
            for leaflet_name, comp in [("upper", upper), ("lower", lower)]:
                for lipid, value in comp.items():
                    if float(value) <= 0:
                        system.errors.append(
                            Message(
                                "LMB-E115",
                                f"Leaflet ratios must be positive. "
                                f"{leaflet_name} leaflet has {lipid}:{value}",
                            )
                        )

        leaflet_lipids = set(upper) | set(lower)
        lipid_specs = set(system.lipids)

        missing_specs = leaflet_lipids - lipid_specs
        if missing_specs:
            system.errors.append(
                Message(
                    "LMB-E116",
                    "Leaflet lipids missing from lipid metadata block: "
                    + ", ".join(sorted(missing_specs)),
                )
            )

        unused_specs = lipid_specs - leaflet_lipids
        if unused_specs:
            system.warnings.append(
                Message(
                    "LMB-W116",
                    "Lipids defined in metadata block but absent from leaflets: "
                    + ", ".join(sorted(unused_specs)),
                )
            )

    def _validate_species_definitions(self, system: LiMBSSystem) -> None:
        chemical_defs = set(system.species_blocks)
        lipid_defs = set(system.lipids)

        leaflet_lipids = set()
        if system.leaflets:
            leaflet_lipids = set(system.leaflets.get("upper", {})) | set(system.leaflets.get("lower", {}))

        # Lipids appearing in leaflets must have species definitions.
        missing_lipid_chem = leaflet_lipids - chemical_defs
        if missing_lipid_chem:
            system.errors.append(
                Message(
                    "LMB-E120",
                    "Leaflet lipids missing species definitions: "
                    + ", ".join(sorted(missing_lipid_chem)),
                )
            )

        # Solvent should have a species definition.
        if system.solvent and system.solvent not in chemical_defs:
            system.warnings.append(
                Message(
                    "LMB-W121",
                    f"Solvent '{system.solvent}' is declared but has no species definition.",
                )
            )

        # Salt validation:
        # 1. If salt is a named species such as XY and XY is defined, accept it.
        # 2. If salt is a formula such as NaCl, CaCl2, or MgCl2, check that the
        #    required ions appear somewhere in the non-lipid species definitions.
        # 3. If salt is written as separate ions, e.g. Na Cl, each token may be
        #    a species block or may appear inside a generic salt block.
        if system.salt:
            self._validate_salt_species(system, chemical_defs, lipid_defs)

        # Extra definitions are fine for solvent/ions. Warn only if they are unrelated.
        for name, block in system.species_blocks.items():
            if block.resolution not in self.SUPPORTED_RESOLUTIONS:
                system.errors.append(
                    Message("LMB-E130", f"{name} uses unsupported resolution {block.resolution}.")
                )

            # If this is a lipid metadata species, its species resolution should usually match system resolution.
            if name in lipid_defs and system.resolution and block.resolution != system.resolution:
                system.warnings.append(
                    Message(
                        "LMB-W130",
                        f"{name} is defined as {block.resolution}, while system resolution is {system.resolution}.",
                    )
                )

            self._validate_fragments(system, name, block)
            self._validate_layers(system, name, block)
            self._validate_semantic_groups(system, name, block)


    def _validate_salt_species(self, system: LiMBSSystem, chemical_defs: set, lipid_defs: set) -> None:
        salt_text = system.salt.get("species_text", "")
        salt_tokens = system.salt.get("species", [])
        salt_ions = system.salt.get("ions")

        # Direct species definition: salt:0.15M XY and XY=CG|[...]
        if len(salt_tokens) == 1 and salt_tokens[0] in chemical_defs:
            return

        # Gather ion labels that appear in non-lipid species blocks.
        available_ion_labels = set()
        for species_name, block in system.species_blocks.items():
            if species_name in lipid_defs:
                continue

            # Bead labels such as ^Na, ^Cl, ^Ca, ^Mg.
            available_ion_labels.update(block.beads)

            # Also detect raw ion text inside the chemical layer, such as [Na+:^Na], [Cl-:^Cl], [Ca2+:^Ca].
            for raw_label in re.findall(r"\[([A-Z][a-z]?)(?:[0-9]*[+-])?:\^[A-Za-z0-9_]+\]", block.body):
                available_ion_labels.add(raw_label)

        # Formula mode: NaCl -> Na and Cl required; CaCl2 -> Ca and Cl required.
        if salt_ions:
            missing_ions = [ion for ion in salt_ions if ion not in available_ion_labels and ion not in chemical_defs]
            if missing_ions:
                system.warnings.append(
                    Message(
                        "LMB-W122",
                        f"Salt '{salt_text}' requires ion labels not found in non-lipid species blocks: "
                        + ", ".join(missing_ions),
                    )
                )
            return

        # Generic label that is not a formula and not a defined species.
        missing_salt = [sp for sp in salt_tokens if sp not in chemical_defs]
        if missing_salt:
            system.warnings.append(
                Message(
                    "LMB-W122",
                    "Salt species declared but not defined as species blocks: "
                    + ", ".join(missing_salt),
                )
            )

    def _validate_fragments(self, system: LiMBSSystem, name: str, block: SpeciesBlock) -> None:
        defined_fragments = set(block.fragments_defined)
        used_fragments = set(block.fragments_used)
        undefined_fragments = used_fragments - defined_fragments

        if undefined_fragments:
            system.errors.append(
                Message(
                    "LMB-E140",
                    f"{name} uses undefined fragments: "
                    + ", ".join(sorted(undefined_fragments)),
                )
            )

    def _validate_layers(self, system: LiMBSSystem, name: str, block: SpeciesBlock) -> None:
        """
        These are warnings rather than hard errors so the parser remains flexible.
        """
        if block.resolution == "AA":
            if block.layer_count < 1:
                system.errors.append(Message("LMB-E150", f"{name}=AA has no layer."))

            if block.layer_count > 1:
                system.warnings.append(
                    Message(
                        "LMB-W150",
                        f"{name}=AA has {block.layer_count} layers. This is allowed, but AA usually has one layer.",
                    )
                )

        elif block.resolution == "CG":
            if block.layer_count < 2:
                system.warnings.append(
                    Message(
                        "LMB-W151",
                        f"{name}=CG has {block.layer_count} layer(s). CG usually has two layers.",
                    )
                )

            if block.layer_count > 2:
                system.warnings.append(
                    Message(
                        "LMB-W152",
                        f"{name}=CG has {block.layer_count} layers. If the third layer is semantic grouping, consider SCG.",
                    )
                )

        elif block.resolution == "SCG":
            if block.layer_count < 3:
                system.warnings.append(
                    Message(
                        "LMB-W153",
                        f"{name}=SCG has {block.layer_count} layer(s). SCG usually has three layers.",
                    )
                )

    def _validate_semantic_groups(self, system: LiMBSSystem, name: str, block: SpeciesBlock) -> None:
        if not block.semantic_groups:
            return

        groups_used = set(block.semantic_groups)
        groups_defined = set(block.semantic_group_definitions)

        undefined_groups = groups_used - groups_defined
        if undefined_groups:
            system.errors.append(
                Message(
                    "LMB-E160",
                    f"{name} uses undefined semantic groups: "
                    + ", ".join(sorted(undefined_groups)),
                )
            )

        defined_fragments = {f"#{x}" for x in block.fragments_defined}
        defined_beads = {f"^{x}" for x in block.beads}
        allowed_members = defined_fragments | defined_beads

        for group, members in block.semantic_group_definitions.items():
            for member in members:
                if member not in allowed_members:
                    system.errors.append(
                        Message(
                            "LMB-E161",
                            f"{name} semantic group ^^{group} contains undefined member {member}.",
                        )
                    )

    def _validate_header_counts(self, system: LiMBSSystem) -> None:
        if system.lipid_species_count is not None and system.lipids:
            actual = len(system.lipids)

            if system.lipid_species_count != actual:
                system.warnings.append(
                    Message(
                        "LMB-W170",
                        f"lsc is {system.lipid_species_count}, but lipid metadata defines {actual} lipid species.",
                    )
                )

        if system.leaflets:
            upper_species = len(system.leaflets.get("upper", {}))
            lower_species = len(system.leaflets.get("lower", {}))

            if (
                system.upper_leaflet_species_count is not None
                and system.upper_leaflet_species_count != upper_species
            ):
                system.warnings.append(
                    Message(
                        "LMB-W171",
                        f"uls is {system.upper_leaflet_species_count}, but upper leaflet contains {upper_species} lipid species.",
                    )
                )

            if (
                system.lower_leaflet_species_count is not None
                and system.lower_leaflet_species_count != lower_species
            ):
                system.warnings.append(
                    Message(
                        "LMB-W172",
                        f"lls is {system.lower_leaflet_species_count}, but lower leaflet contains {lower_species} lipid species.",
                    )
                )


# Export / canonicalization

def system_to_json(system: LiMBSSystem) -> Dict[str, Any]:
    data = asdict(system)
    data["valid"] = len(system.errors) == 0
    data["errors"] = [asdict(e) for e in system.errors]
    data["warnings"] = [asdict(w) for w in system.warnings]
    return data
def canonical_number(value, max_decimals=6) -> str:
    text = f"{float(value):.{max_decimals}f}"
    text = text.rstrip("0").rstrip(".")
    if "." not in text:
        text += ".0"
    return text


def canonicalize(system: LiMBSSystem) -> str:
    if system.errors:
        raise LiMBSError("Cannot canonicalize invalid LiMBS system.")

    lipid_entries = []
    for lipid in sorted(system.lipids):
        data = system.lipids[lipid]
        lipid_entries.append(
            f"{lipid}:tails={data['tails']},head={data['head']},charge={data['charge']:.1f}"
        )

    def comp_string(comp: Dict[str, Any]) -> str:
        return ",".join(f"{lipid}:{comp[lipid]}" for lipid in sorted(comp))

    leaf = system.leaflets

    system_parts = [
        "LiMBS1.0",
        f"rsln:{system.resolution}",
    ]

    if system.lipid_species_count is not None:
        system_parts.append(f"lsc:{system.lipid_species_count}")

    if system.upper_leaflet_species_count is not None:
        system_parts.append(f"uls:{system.upper_leaflet_species_count}")

    if system.lower_leaflet_species_count is not None:
        system_parts.append(f"lls:{system.lower_leaflet_species_count}")

    system_parts.extend(
        [
            f"{{{';'.join(lipid_entries)}}}",
            f"leaflet:{leaf['mode']}{{-u{{{comp_string(leaf['upper'])}}}-l{{{comp_string(leaf['lower'])}}}}}",
            f"type:{system.membrane_type}",
            (
                f"box:[{system.box['values'][0]},"
                f"{system.box['values'][1]},"
                f"{system.box['values'][2]}]"
                f"{system.box['unit']}"
            ),
            f"sol:{system.solvent}",
        ]
    )

    if system.salt:
        system_parts.append(
            f"salt:{system.salt['concentration']:.2f}{system.salt['unit']} {system.salt['species_text']}"
        )

    species_parts = []

    def species_sort_key(name: str):
        if name in system.lipids:
            return (0, name)
        return (1, name)

    for name in sorted(system.species_blocks, key=species_sort_key):
        block = system.species_blocks[name]
        species_parts.append(f"{name}={block.resolution}|{block.body}")

    return " ||\n".join(system_parts) + "\n\n|++|\n\n" + "\n\n|+|\n\n".join(species_parts)

def print_summary(system: LiMBSSystem) -> None:
    print("\nLiMBS1.0 parsed system")
    print(f"  Version          : {system.version}")
    print(f"  Resolution       : {system.resolution}")
   # print(f"  Force field      : {system.forcefield}")
    print(f"  Lipid sp. count  : {system.lipid_species_count}")
    print(f"  U leaflet species: {system.upper_leaflet_species_count}")
    print(f"  L leaflet species: {system.lower_leaflet_species_count}")
    print(f"  Type             : {system.membrane_type}")
    print(f"  Box              : {system.box.get('values')} {system.box.get('unit')}")
    print(f"  Solvent          : {system.solvent}")

    if system.salt:
        print(
            f"  Salt             : {system.salt['concentration']} "
            f"{system.salt['unit']} {system.salt['species_text']}"
        )

    if system.leaflets:
        print(f"  Leaflet mode     : {system.leaflets['mode']}")
        print(f"  Upper leaflet    : {system.leaflets['upper']}")
        print(f"  Lower leaflet    : {system.leaflets['lower']}")

    print(f"  Lipid metadata   : {list(system.lipids.keys())}")
    print(f"  Species defs     : {list(system.species_blocks.keys())}")

    for name, block in system.species_blocks.items():
        print(
            f"  {name:<12}: {block.resolution}, "
            f"{block.layer_count} layer(s)"
        )

        if block.fragments_defined:
            print(f"      Fragments defined: {block.fragments_defined}")

        if block.fragments_used:
            print(f"      Fragments used   : {block.fragments_used}")

        if block.semantic_groups:
            print(f"      Semantic groups  : {block.semantic_groups}")

    if system.warnings:
        print("\nWarnings:")
        for warning in system.warnings:
            print(f"  {warning.code}: {warning.message}")

    if system.errors:
        print("\nErrors:")
        for error in system.errors:
            print(f"  {error.code}: {error.message}")

    print()


# CLI

def main() -> None:
    ap = argparse.ArgumentParser(
        description="LiMBS1.0 parser, validator, JSON exporter, and canonicalizer for AA/CG/SCG notation."
    )

    ap.add_argument("input", help="Input LiMBS text file")
    ap.add_argument("--validate-only", action="store_true", help="Only validate the LiMBS string")
    ap.add_argument("--json", action="store_true", help="Print parsed JSON")
    ap.add_argument("--canonical", action="store_true", help="Print canonical LiMBS string")
    ap.add_argument("--quiet", action="store_true", help="Suppress summary output")

    args = ap.parse_args()

    try:
        with open(args.input, "r", encoding="utf-8") as fh:
            text = fh.read()
    except FileNotFoundError:
        print(f"LMB-E900: File not found: {args.input}")
        sys.exit(1)

    parser = LiMBSParser()

    try:
        system = parser.parse(text)
    except LiMBSError as exc:
        print(str(exc))
        sys.exit(1)


    if system.errors:
        print("Validation failed.")
        sys.exit(1)

    if args.json:
        print(json.dumps(system_to_json(system), indent=2))
        return

    if args.canonical:
        try:
            print(canonicalize(system))
        except LiMBSError as exc:
            print(str(exc))
            sys.exit(1)
        return

    if not args.quiet:
        print_summary(system)

    print("Validation passed.")

    if args.validate_only:
        return


if __name__ == "__main__":
    main()
