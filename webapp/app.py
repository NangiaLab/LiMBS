#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LiMBS1.0 web application.

Current scope
Generate LiMBS1.0 notation from the web form.
Generate resolution-specific lipid/environment species blocks.
Validate generated notation with LiMBSv1_parser.py.
Display canonical LiMBS and structured JSON.
Provide manual or assisted box-dimension guidance.
Use the curated LiMBS APL reference library when the user explicitly chooses it.
For mixed planar membranes, provide a clearly labelled composition-weighted starting estimate only when every component has a suitable UI reference.
Recommend a downstream construction interface only for supported CG systems:
     CG planar  -> INSANE
     CG vesicle -> TS2CG
Provide downloadable LiMBS builder-interface scripts.

Reference APL values are construction aids, not universal lipid constants or equilibrium membrane properties.  Mixed-membrane estimates are explicitly labelled approximations.  User-provided APL values always remain available.
The web application generates previews and files; it does not execute INSANE or TS2CG.
"""

from __future__ import annotations

import json
import math
import re
import shlex
from html import escape
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from flask import Flask, abort, render_template, request, send_file, url_for

from lipid_registry import LIPIDS
from cg_templates import generate_cg_block


# Optional curated APL reference library

try:
    from apl_reference_combined import get_reference_for_ui
except ImportError:
    get_reference_for_ui = None


# Optional AA / SCG template libraries

try:
    from aa_templates import generate_aa_block, build_aa_environment_blocks
except ImportError:
    generate_aa_block = None
    build_aa_environment_blocks = None

try:
    from scg_templates import generate_scg_block
except ImportError:
    generate_scg_block = None


# LiMBS parser

try:
    from LiMBSv1_parser import (
        LiMBSParser,
        LiMBSError,
        canonicalize,
        system_to_json,
    )
except ImportError as exc:
    raise RuntimeError(
        "Could not import LiMBSv1_parser.py. "
        "Place LiMBSv1_parser.py in the same directory as app.py."
    ) from exc


app = Flask(__name__)
PROJECT_DIR = Path(__file__).resolve().parent


# General helpers


def error_page(title: str, message: str, status_code: int = 400):
    """Return a simple, readable application error page."""
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{escape(title)}</title>
        <style>
            body {{
                font-family: Arial, Helvetica, sans-serif;
                margin: 40px;
                line-height: 1.5;
                color: #222;
            }}
            .error-box {{
                max-width: 900px;
                border: 1px solid #d0d0d0;
                border-left: 5px solid #b42318;
                border-radius: 6px;
                padding: 18px;
                background: #fff8f7;
            }}
            pre {{
                white-space: pre-wrap;
                word-break: break-word;
            }}
            a {{
                color: #1a5fb4;
                text-decoration: none;
            }}
        </style>
    </head>
    <body>
        <div class="error-box">
            <h2>{escape(title)}</h2>
            <pre>{escape(message)}</pre>
            <p><a href="/">Return to generator</a></p>
        </div>
    </body>
    </html>
    """
    return html, status_code


def format_decimal(value: float, max_decimals: int = 6) -> str:
    """
    Format a numeric value without unnecessary trailing zeroes.

    Examples:
        0.800000 -> 0.8
        0.200000 -> 0.2
        8.062000 -> 8.062
    """
    text = f"{float(value):.{max_decimals}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def format_box_value(value: float) -> str:
    """
    Format box dimensions for generated LiMBS.

    Whole-number dimensions retain one decimal place:
        10 -> 10.0
        26 -> 26.0

    Non-integer values retain useful precision:
        8.062 -> 8.062
    """
    value = float(value)

    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1e-9):
        return f"{value:.1f}"

    return f"{value:.6f}".rstrip("0").rstrip(".")


def parse_positive_float(raw: str, field_name: str, allow_zero: bool = False) -> float:
    """Parse a positive numeric form value."""
    try:
        value = float(str(raw).strip())
    except (TypeError, ValueError):
        raise ValueError(f"{field_name} must be a number.")

    if not math.isfinite(value):
        raise ValueError(f"{field_name} must be finite.")

    if allow_zero:
        if value < 0:
            raise ValueError(f"{field_name} cannot be negative.")
    elif value <= 0:
        raise ValueError(f"{field_name} must be greater than zero.")

    return value


# APL and composition helpers


def active_composition(composition: Dict[str, Any]) -> Dict[str, float]:
    """Return positive lipid entries as normalized uppercase keys."""
    clean: Dict[str, float] = {}

    for lipid, raw_value in composition.items():
        value = float(raw_value)
        if value > 0:
            clean[str(lipid).strip().upper()] = value

    return clean


def default_apl_model(resolution: str) -> Optional[str]:
    """Return the curated reference model used for assisted APL lookup."""
    resolution = str(resolution).strip().upper()

    if resolution == "CG":
        return "Martini3"
    if resolution == "AA":
        return "CHARMM36"
    return None


def build_apl_availability_map() -> Dict[str, Dict[str, bool]]:
    """
    Return per-lipid APL-reference availability for the web interface.

    This reports only whether a suitable UI-recommended reference exists.
    It does not describe structural/template support.  In particular, a lipid
    such as CHOL can be structurally supported while intentionally having no
    standalone APL reference.
    """
    availability: Dict[str, Dict[str, bool]] = {}

    for lipid in LIPIDS:
        availability[lipid] = {
            "CG": False,
            "AA": False,
            "SCG": False,
        }

        if get_reference_for_ui is None:
            continue

        for resolution, model in (("CG", "Martini3"), ("AA", "CHARMM36")):
            try:
                ref = get_reference_for_ui(
                    lipid=lipid,
                    resolution=resolution,
                    model=model,
                )
            except (KeyError, TypeError, ValueError):
                ref = None

            availability[lipid][resolution] = (
                ref is not None and ref.get("value_nm2") is not None
            )

    return availability


def read_optional_temperature() -> Optional[float]:
    """Read an optional target temperature in K for reference selection."""
    raw = request.form.get("temperature_k", "").strip()
    if not raw:
        return None

    return parse_positive_float(raw, "Temperature")


def get_leaflet_reference_apl(
    composition: Dict[str, Any],
    resolution: str,
    model: str,
    target_temperature_K: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Return a pure-lipid reference or a composition-weighted starting estimate.

    The mixed estimate is calculated only when every lipid has an appropriate
    UI reference.  Missing/caution-only references are never silently replaced.
    """
    comp = active_composition(composition)

    if not comp:
        return {
            "available": False,
            "reason": "Leaflet composition is empty.",
            "missing_lipids": [],
        }

    if get_reference_for_ui is None:
        return {
            "available": False,
            "reason": (
                "The curated APL reference library is not installed. "
                "Choose a user-specified APL."
            ),
            "missing_lipids": list(comp),
        }

    total = sum(comp.values())
    components: List[Dict[str, Any]] = []
    missing: List[str] = []

    for lipid, amount in sorted(comp.items()):
        ref = get_reference_for_ui(
            lipid=lipid,
            resolution=resolution,
            model=model,
            target_temperature_K=target_temperature_K,
        )

        if ref is None or ref.get("value_nm2") is None:
            missing.append(lipid)
            continue

        fraction = amount / total
        components.append(
            {
                "lipid": lipid,
                "amount": amount,
                "fraction": fraction,
                "value_nm2": float(ref["value_nm2"]),
                "temperature_K": ref.get("temperature_K"),
                "model": ref.get("model", model),
                "source": ref.get("source", {}),
                "display_label": ref.get("display_label"),
            }
        )

    if missing:
        return {
            "available": False,
            "reason": (
                "No suitable LiMBS UI reference APL is available for: "
                + ", ".join(missing)
                + ". Choose a user-specified APL for this membrane."
            ),
            "missing_lipids": missing,
            "components": components,
        }

    if len(components) == 1:
        value_nm2 = components[0]["value_nm2"]
        method = "pure-lipid-reference"
        label = "LiMBS reference APL"
    else:
        value_nm2 = sum(
            item["fraction"] * item["value_nm2"]
            for item in components
        )
        method = "composition-weighted-estimate"
        label = "LiMBS composition-weighted starting estimate"

    return {
        "available": True,
        "value_nm2": value_nm2,
        "method": method,
        "label": label,
        "components": components,
        "missing_lipids": [],
        "warning": (
            "Reference APL values are construction aids. A composition-weighted "
            "value is an initial approximation derived from pure-lipid references "
            "and is not an equilibrium APL for the mixed membrane."
        ),
    }


def read_user_apl(leaflet_name: str) -> float:
    """
    Read a user-supplied APL.

    New forms may provide leaflet-specific fields.  The shared custom_apl field
    and the legacy reference_apl field are retained as fallbacks.
    """
    candidates = [
        request.form.get(f"custom_apl_{leaflet_name}", "").strip(),
        request.form.get("custom_apl", "").strip(),
        request.form.get("reference_apl", "").strip(),
    ]

    raw = next((value for value in candidates if value), "")
    if not raw:
        raise ValueError(
            "Enter a user-specified area per lipid in nm²/lipid."
        )

    label = f"{leaflet_name.capitalize()} leaflet area per lipid"
    return parse_positive_float(raw, label)


def lateral_dimensions_from_area(
    area_nm2: float,
    shape: str,
    aspect_ratio: float,
) -> Tuple[float, float]:
    """Convert a target XY area into square, rectangular, or cubic lateral dimensions."""
    if shape in {"square", "cubic"}:
        side = math.sqrt(area_nm2)
        return side, side

    if shape == "rectangular":
        x = math.sqrt(area_nm2 * aspect_ratio)
        y = math.sqrt(area_nm2 / aspect_ratio)
        return x, y

    raise ValueError(
        "Planar box shape must be square, rectangular, or cubic."
    )


# Lipid metadata helpers


def get_tail_descriptor(lipid: str, data: Dict[str, Any]) -> str:
    """
    Return the human-facing tail descriptor used in LiMBS metadata.

    The registry header is used first so the displayed biological lipid
    naming remains independent of the internal Martini tail_a/tail_b ordering.
    """
    header = str(data.get("header", ""))

    match = re.match(r"^[^:]+:([^,]+)", header)
    if match:
        return match.group(1).strip()

    tail_a = data.get("tail_a")
    tail_b = data.get("tail_b")

    if tail_a is None and tail_b is None:
        return "sterol"

    if tail_a is not None and tail_b is not None:
        return (
            f"{tail_a[0]}:{tail_a[1]}-"
            f"{tail_b[0]}:{tail_b[1]}"
        )

    return "unknown"


def get_head_descriptor(lipid: str, data: Dict[str, Any]) -> str:
    """Return the compact LiMBS headgroup label."""
    lipid_class = str(data.get("class", "")).upper()

    if lipid.upper() == "CHOL" or lipid_class == "STEROL":
        return "CHOL"

    return lipid_class


def build_lipid_metadata(selected_lipids: List[str]) -> str:
    """Build the {...} LiMBS lipid metadata block."""
    entries = []

    for lipid in sorted(selected_lipids):
        data = LIPIDS[lipid]
        tails = get_tail_descriptor(lipid, data)
        head = get_head_descriptor(lipid, data)
        charge = float(data.get("charge", 0.0))

        entries.append(
            f"{lipid}:tails={tails},head={head},charge={charge:.1f}"
        )

    return "{" + ";".join(entries) + "}"


# Composition


def read_lipid_composition(
    leaflet_mode: str,
) -> Tuple[List[str], Dict[str, Any], Dict[str, Any]]:
    """
    Read upper/outer and lower/inner lipid composition from the submitted form.

    Zero means "lipid absent from this leaflet" and is omitted from LiMBS.
    """
    if leaflet_mode not in {"count", "ratio"}:
        raise ValueError("Leaflet composition mode must be count or ratio.")

    selected_lipids: List[str] = []
    upper: Dict[str, Any] = {}
    lower: Dict[str, Any] = {}

    for lipid in LIPIDS:
        upper_raw = request.form.get(f"{lipid}_upper", "0").strip()
        lower_raw = request.form.get(f"{lipid}_lower", "0").strip()

        try:
            upper_number = float(upper_raw or "0")
            lower_number = float(lower_raw or "0")
        except ValueError:
            raise ValueError(
                f"{lipid}: leaflet values must be numeric."
            )

        if not math.isfinite(upper_number) or not math.isfinite(lower_number):
            raise ValueError(f"{lipid}: leaflet values must be finite.")

        if upper_number < 0 or lower_number < 0:
            raise ValueError(f"{lipid}: leaflet values cannot be negative.")

        if leaflet_mode == "count":
            if (
                upper_number > 0
                and not float(upper_number).is_integer()
            ):
                raise ValueError(
                    f"{lipid}: count mode requires integer molecule counts."
                )

            if (
                lower_number > 0
                and not float(lower_number).is_integer()
            ):
                raise ValueError(
                    f"{lipid}: count mode requires integer molecule counts."
                )

            upper_value = int(upper_number)
            lower_value = int(lower_number)

        else:
            upper_value = float(upper_number)
            lower_value = float(lower_number)

        if upper_number > 0 or lower_number > 0:
            selected_lipids.append(lipid)

        if upper_number > 0:
            upper[lipid] = upper_value

        if lower_number > 0:
            lower[lipid] = lower_value

    if not selected_lipids:
        raise ValueError(
            "No lipids were selected. Enter at least one non-zero lipid value."
        )

    if not upper:
        raise ValueError(
            "The upper/outer leaflet is empty. "
            "Enter at least one non-zero lipid value."
        )

    if not lower:
        raise ValueError(
            "The lower/inner leaflet is empty. "
            "Enter at least one non-zero lipid value."
        )

    return selected_lipids, upper, lower


def build_leaflet_block(
    leaflet_mode: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
) -> str:
    """Build the LiMBS leaflet block."""

    def format_composition(comp: Dict[str, Any]) -> str:
        parts = []

        for lipid in sorted(comp):
            value = comp[lipid]

            if leaflet_mode == "count":
                value_text = str(int(value))
            else:
                value_text = format_decimal(float(value))

            parts.append(f"{lipid}:{value_text}")

        return ",".join(parts)

    return (
        f"leaflet:{leaflet_mode}"
        f"{{-u{{{format_composition(upper)}}}"
        f"-l{{{format_composition(lower)}}}}}"
    )


# Species-definition generation


CG_ENVIRONMENT_BLOCKS = {
    "W": """W=CG|[
[OOOO:^W]
|
[^W]
]""",

    "NaCl": """NaCl=CG|[
[Na+:^Na][Cl-:^Cl]
|
[^Na][^Cl]
]""",

    "KCl": """KCl=CG|[
[K+:^K][Cl-:^Cl]
|
[^K][^Cl]
]""",

    "CaCl2": """CaCl2=CG|[
[Ca2+:^Ca][Cl-:^Cl][Cl-:^Cl]
|
[^Ca][^Cl][^Cl]
]""",

    "MgCl2": """MgCl2=CG|[
[Mg2+:^Mg][Cl-:^Cl][Cl-:^Cl]
|
[^Mg][^Cl][^Cl]
]""",
}


def build_cg_environment_blocks(solvent: str, salt_species: str) -> List[str]:
    """Build CG solvent and salt species definitions."""
    blocks: List[str] = []

    if salt_species:
        if salt_species not in CG_ENVIRONMENT_BLOCKS:
            raise ValueError(
                f"No CG LiMBS environment template is currently defined "
                f"for salt '{salt_species}'."
            )
        blocks.append(CG_ENVIRONMENT_BLOCKS[salt_species])

    if solvent:
        if solvent not in CG_ENVIRONMENT_BLOCKS:
            raise ValueError(
                f"No CG LiMBS environment template is currently defined "
                f"for solvent '{solvent}'."
            )
        blocks.append(CG_ENVIRONMENT_BLOCKS[solvent])

    return blocks


def build_species_blocks(
    resolution: str,
    selected_lipids: List[str],
    solvent: str,
    salt_species: str,
) -> str:
    """
    Build all species definitions after |++|.

    AA and SCG are deliberately not fabricated from CG templates.
    """
    blocks: List[str] = []

    if resolution == "CG":
        for lipid in selected_lipids:
            data = LIPIDS[lipid]

            blocks.append(
                generate_cg_block(
                    lipid,
                    data["class"],
                    data.get("tail_a"),
                    data.get("tail_b"),
                )
            )

        blocks.extend(
            build_cg_environment_blocks(
                solvent=solvent,
                salt_species=salt_species,
            )
        )

    elif resolution == "AA":
        if generate_aa_block is None or build_aa_environment_blocks is None:
            raise ValueError(
                "AA was selected, but aa_templates.py is not installed. "
                "The application will not invent atomistic species definitions."
            )

        for lipid in selected_lipids:
            blocks.append(
                generate_aa_block(lipid, LIPIDS[lipid])
            )
        blocks.extend(
            build_aa_environment_blocks(
                solvent=solvent,
                salt_species=salt_species,
            )
        )

    elif resolution == "SCG":
        if generate_scg_block is None:
            raise ValueError(
                "SCG was selected, but scg_templates.py is not installed. "
                "The application will not invent semantic-group definitions."
            )

        for lipid in selected_lipids:
            blocks.append(generate_scg_block(lipid, LIPIDS[lipid]))

    else:
        raise ValueError(f"Unsupported resolution: {resolution}")

    return "\n\n|+|\n\n".join(blocks)


# Solvent / salt


def read_solvent(resolution: str) -> str:
    """Read resolution-aware solvent selection."""
    solvent = request.form.get("sol", "").strip()

    if resolution == "AA" and solvent == "OTHER":
        solvent = request.form.get("custom_solvent", "").strip()

        if not solvent:
            raise ValueError(
                "A custom solvent identifier is required when 'Other' is selected."
            )

    if not solvent:
        raise ValueError("Solvent cannot be empty.")

    return solvent


def read_salt() -> Tuple[float, str]:
    """Read salt concentration and salt name."""
    salt_conc = parse_positive_float(
        request.form.get("salt_conc", "0.15"),
        "Salt concentration",
        allow_zero=True,
    )

    salt_species = request.form.get("salt_ion", "NaCl").strip()

    if not salt_species:
        raise ValueError("Salt species cannot be empty.")

    return salt_conc, salt_species


# Box dimensions


def determine_box(
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
) -> Tuple[List[float], str, Dict[str, Any]]:
    """
    Determine simulation-box dimensions.

    Assisted planar dimensions may use either:
      1. an explicit user APL, or
      2. a LiMBS reference/composition-weighted starting estimate.

    Assisted planar geometry may be square, rectangular, or cubic.
    All assisted values are construction starting estimates, not equilibrium
    membrane properties.
    """
    box_mode = request.form.get("box_mode", "manual").strip().lower()

    if box_mode not in {"manual", "estimate"}:
        raise ValueError("Dimension method must be manual or estimate.")

    box_info: Dict[str, Any] = {
        "mode": box_mode,
        "membrane_type": membrane_type,
        "resolution": resolution,
    }

    # Manual dimensions
    if box_mode == "manual":
        x = parse_positive_float(request.form.get("box_x", "10.0"), "X dimension")
        y = parse_positive_float(request.form.get("box_y", "10.0"), "Y dimension")
        z = parse_positive_float(request.form.get("box_z", "8.0"), "Z dimension")

        box = [x, y, z]

        if membrane_type == "planar":
            area = x * y
            box_info["projected_area"] = area

            if leaflet_mode == "count":
                upper_count = sum(int(v) for v in upper.values())
                lower_count = sum(int(v) for v in lower.values())
                upper_apl = area / upper_count
                lower_apl = area / lower_count

                box_info["upper_implied_apl"] = upper_apl
                box_info["lower_implied_apl"] = lower_apl

                guidance = (
                    "Manual dimensions were used exactly as entered. "
                    f"The projected XY area is {area:.3f} nm². "
                    f"The implied APL is {upper_apl:.3f} nm²/lipid for the "
                    f"upper leaflet and {lower_apl:.3f} nm²/lipid for the lower "
                    "leaflet. These are diagnostics only."
                )
            else:
                guidance = (
                    "Manual dimensions were used exactly as entered. "
                    f"The projected XY area is {area:.3f} nm². "
                    "Ratio mode does not define a final lipid population, so an "
                    "implied APL cannot be obtained from the box alone."
                )
        else:
            guidance = (
                "Manual vesicle box dimensions were used exactly as entered. "
                "Choose dimensions large enough for the complete vesicle and "
                "adequate solvent separation from periodic images."
            )

        return box, guidance, box_info

    # Assisted planar estimate
    if membrane_type == "planar":
        shape = request.form.get(
            "planar_shape",
            request.form.get("box_shape", "square"),
        ).strip().lower()

        if shape not in {"square", "rectangular", "cubic"}:
            raise ValueError(
                "Planar box shape must be square, rectangular, or cubic."
            )

        if shape == "rectangular":
            aspect_ratio = parse_positive_float(
                request.form.get(
                    "xy_aspect_ratio",
                    request.form.get("aspect_ratio", "1.0"),
                ),
                "X:Y aspect ratio",
            )
        else:
            aspect_ratio = 1.0

        # For square/rectangular planar boxes, Z is independently supplied.
        # For a cubic box, Z is linked to the lateral side and therefore does
        # not need a separate user input.
        if shape == "cubic":
            z = None
        else:
            z = parse_positive_float(
                request.form.get("estimated_planar_z", "8.0"),
                "Planar Z dimension",
            )

        raw_apl_source = request.form.get("apl_source")
        if raw_apl_source is None:
            apl_source = "custom" if request.form.get("reference_apl") else "limbs"
        else:
            apl_source = raw_apl_source.strip().lower()

        if apl_source in {"reference", "recommended", "limbs_reference"}:
            apl_source = "limbs"
        if apl_source in {"own", "manual", "user"}:
            apl_source = "custom"

        if apl_source not in {"limbs", "custom"}:
            raise ValueError(
                "APL source must be LiMBS reference/estimate or user-specified."
            )

        target_temperature = read_optional_temperature()
        model = request.form.get("apl_model", "").strip() or default_apl_model(resolution)

        if apl_source == "limbs":
            if model is None:
                raise ValueError(
                    f"LiMBS does not currently provide a curated assisted APL "
                    f"reference model for {resolution}. Choose a user-specified APL."
                )

            upper_apl_info = get_leaflet_reference_apl(
                upper,
                resolution=resolution,
                model=model,
                target_temperature_K=target_temperature,
            )
            lower_apl_info = get_leaflet_reference_apl(
                lower,
                resolution=resolution,
                model=model,
                target_temperature_K=target_temperature,
            )

            unavailable = []
            if not upper_apl_info.get("available"):
                unavailable.append("Upper leaflet: " + upper_apl_info["reason"])
            if not lower_apl_info.get("available"):
                unavailable.append("Lower leaflet: " + lower_apl_info["reason"])

            if unavailable:
                raise ValueError(
                    "LiMBS-assisted APL is unavailable for this membrane. "
                    + " ".join(unavailable)
                )

            upper_apl = float(upper_apl_info["value_nm2"])
            lower_apl = float(lower_apl_info["value_nm2"])

        else:
            upper_apl = read_user_apl("upper")
            lower_apl = read_user_apl("lower")
            upper_apl_info = {
                "available": True,
                "value_nm2": upper_apl,
                "method": "user-specified",
                "label": "User-specified APL",
                "components": [],
            }
            lower_apl_info = {
                "available": True,
                "value_nm2": lower_apl,
                "method": "user-specified",
                "label": "User-specified APL",
                "components": [],
            }

        if leaflet_mode == "count":
            upper_count = sum(int(v) for v in upper.values())
            lower_count = sum(int(v) for v in lower.values())
        else:
            target_count = int(
                parse_positive_float(
                    request.form.get("target_lipids_per_leaflet", "100"),
                    "Target lipids per leaflet",
                )
            )
            upper_count = target_count
            lower_count = target_count
            box_info["target_lipids_per_leaflet"] = target_count

        upper_area = upper_count * upper_apl
        lower_area = lower_count * lower_apl
        target_area = max(upper_area, lower_area)

        mean_area = (upper_area + lower_area) / 2.0
        area_mismatch_fraction = (
            abs(upper_area - lower_area) / mean_area
            if mean_area > 0
            else 0.0
        )

        x, y = lateral_dimensions_from_area(
            target_area,
            shape=shape,
            aspect_ratio=aspect_ratio,
        )

        if shape == "cubic":
            z = x

        source_phrase = (
            "LiMBS reference/composition-weighted starting APL"
            if apl_source == "limbs"
            else "user-specified APL"
        )

        if shape == "cubic":
            dimension_phrase = (
                f"X = Y = Z = {x:.3f} nm. "
                "Because cubic mode links the membrane-normal dimension to the "
                "lateral side, verify that the resulting Z dimension provides "
                "sufficient solvent separation for the intended system. "
            )
        elif shape == "square":
            dimension_phrase = (
                f"X = Y = {x:.3f} nm and Z = {z:.3f} nm. "
            )
        else:
            dimension_phrase = (
                f"X = {x:.3f} nm, Y = {y:.3f} nm, and Z = {z:.3f} nm. "
            )

        guidance = (
            f"Assisted planar dimensions were calculated using the {source_phrase}. "
            f"Upper starting APL: {upper_apl:.5f} nm²/lipid; "
            f"lower starting APL: {lower_apl:.5f} nm²/lipid. "
            f"The corresponding preferred leaflet areas are {upper_area:.3f} and "
            f"{lower_area:.3f} nm². LiMBS used the larger area ({target_area:.3f} "
            f"nm²) to avoid undersizing the starting lateral area, giving "
            + dimension_phrase
            + "These are construction values and should not be interpreted as "
            "equilibrium membrane properties."
        )

        if area_mismatch_fraction > 0.10:
            guidance += (
                f" The two leaflets imply areas that differ by "
                f"{100.0 * area_mismatch_fraction:.1f}%; consider adjusting "
                "leaflet counts/composition or supplying leaflet-specific APL values."
            )

        box_info.update(
            {
                "apl_source": apl_source,
                "apl_model": model,
                "temperature_K": target_temperature,
                "planar_shape": shape,
                "xy_aspect_ratio": aspect_ratio,
                "z_linked_to_lateral": (shape == "cubic"),
                "upper_apl": upper_apl,
                "lower_apl": lower_apl,
                "upper_apl_info": upper_apl_info,
                "lower_apl_info": lower_apl_info,
                "upper_preferred_area": upper_area,
                "lower_preferred_area": lower_area,
                "area_mismatch_fraction": area_mismatch_fraction,
                "projected_area": target_area,
            }
        )

        return [x, y, z], guidance, box_info

    # Assisted vesicle estimate
    diameter = parse_positive_float(
        request.form.get("vesicle_diameter", "20.0"),
        "Vesicle diameter",
    )
    padding = parse_positive_float(
        request.form.get("vesicle_padding", "3.0"),
        "Vesicle padding",
        allow_zero=True,
    )

    side = diameter + (2.0 * padding)
    box_info.update(
        {
            "vesicle_diameter": diameter,
            "vesicle_padding": padding,
        }
    )

    guidance = (
        "The assisted vesicle box uses vesicle diameter + 2 × solvent padding. "
        f"For a diameter of {format_decimal(diameter)} nm and "
        f"{format_decimal(padding)} nm padding per side, the starting box is "
        f"{format_decimal(side)} × {format_decimal(side)} × "
        f"{format_decimal(side)} nm. This is a construction estimate; verify "
        "solvent separation for the final built vesicle."
    )

    return [side, side, side], guidance, box_info


# LiMBS construction


def build_system_block(
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    selected_lipids: List[str],
    upper: Dict[str, Any],
    lower: Dict[str, Any],
    box: List[float],
    solvent: str,
    salt_conc: float,
    salt_species: str,
) -> str:
    """Build the LiMBS system-level block."""
    metadata = build_lipid_metadata(selected_lipids)
    leaflet_block = build_leaflet_block(leaflet_mode, upper, lower)

    box_text = ",".join(format_box_value(v) for v in box)

    parts = [
        "LiMBS1.0",
        f"rsln:{resolution}",
        f"lsc:{len(selected_lipids)}",
        f"uls:{len(upper)}",
        f"lls:{len(lower)}",
        metadata,
        leaflet_block,
        f"type:{membrane_type}",
        f"box:[{box_text}]nm",
        f"sol:{solvent}",
        f"salt:{salt_conc:.2f}M {salt_species}",
    ]

    return " ||\n".join(parts)


def build_limbs_notation(
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    selected_lipids: List[str],
    upper: Dict[str, Any],
    lower: Dict[str, Any],
    box: List[float],
    solvent: str,
    salt_conc: float,
    salt_species: str,
) -> str:
    """Build complete LiMBS notation."""
    system_block = build_system_block(
        resolution=resolution,
        membrane_type=membrane_type,
        leaflet_mode=leaflet_mode,
        selected_lipids=selected_lipids,
        upper=upper,
        lower=lower,
        box=box,
        solvent=solvent,
        salt_conc=salt_conc,
        salt_species=salt_species,
    )

    species_block = build_species_blocks(
        resolution=resolution,
        selected_lipids=selected_lipids,
        solvent=solvent,
        salt_species=salt_species,
    )

    return f"{system_block}\n\n|++|\n\n{species_block}"


# Parser validation


def validate_generated_notation(notation: str) -> Dict[str, Any]:
    """Parse, validate, canonicalize, and export JSON."""
    parser = LiMBSParser()

    try:
        system = parser.parse(notation)
    except LiMBSError as exc:
        return {
            "passed": False,
            "message": f"LiMBS parser error: {exc}",
            "canonical": None,
            "json": None,
            "warnings": [],
        }
    except Exception as exc:
        return {
            "passed": False,
            "message": f"Unexpected parser error: {exc}",
            "canonical": None,
            "json": None,
            "warnings": [],
        }

    data = system_to_json(system)

    errors = data.get("errors", [])
    warnings = data.get("warnings", [])

    if errors:
        error_text = "\n".join(
            f"[{item.get('code', 'ERROR')}] {item.get('message', item)}"
            for item in errors
        )

        return {
            "passed": False,
            "message": error_text,
            "canonical": None,
            "json": data,
            "warnings": warnings,
        }

    try:
        canonical = canonicalize(system)
    except Exception as exc:
        return {
            "passed": False,
            "message": f"Canonicalization failed: {exc}",
            "canonical": None,
            "json": data,
            "warnings": warnings,
        }

    return {
        "passed": True,
        "message": (
            "Generated LiMBS1.0 specification passed parser validation."
        ),
        "canonical": canonical,
        "json": data,
        "warnings": warnings,
    }


# INSANE workflow preview


def insane_compatibility(
    resolution: str,
    membrane_type: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
) -> Tuple[bool, List[str]]:
    """Check whether the current system matches the website's INSANE interface."""
    reasons = []

    if resolution != "CG":
        reasons.append("The current INSANE interface is limited to CG systems.")

    if membrane_type != "planar":
        reasons.append(
            "The current INSANE interface is intended for planar membranes."
        )

    if not upper or not lower:
        reasons.append("Both leaflets must contain at least one lipid.")

    return (len(reasons) == 0), reasons


def build_insane_command_preview(
    leaflet_mode: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
    box: List[float],
    solvent: str,
    salt_conc: float,
) -> Tuple[str, str]:
    """
    Build a conservative INSANE preview.

    INSANE's -u/-l syntax is expressed as relative abundance in current
    documentation. Therefore count-mode LiMBS values are passed as relative
    weights in this preview rather than claimed to be exact molecule counts.
    """
    args = ["python3", "insane.py"]

    for lipid in sorted(upper):
        args.extend(
            [
                "-u",
                f"{lipid}:{format_decimal(float(upper[lipid]))}",
            ]
        )

    for lipid in sorted(lower):
        args.extend(
            [
                "-l",
                f"{lipid}:{format_decimal(float(lower[lipid]))}",
            ]
        )

    args.extend(
        [
            "-x", format_box_value(box[0]),
            "-y", format_box_value(box[1]),
            "-z", format_box_value(box[2]),
            "-sol", solvent,
            "-salt", format_decimal(salt_conc),
            "-o", "membrane.gro",
            "-p", "membrane.top",
        ]
    )

    command = " ".join(shlex.quote(item) for item in args)

    if leaflet_mode == "count":
        note = (
            "LiMBS count mode records absolute molecule counts. "
            "The displayed INSANE command uses those numbers as composition "
            "weights in the builder preview; it should not be interpreted as "
            "a guarantee of exact final lipid counts."
        )
    else:
        note = (
            "LiMBS ratio values are passed to INSANE as relative leaflet "
            "abundances. Final construction remains controlled by INSANE."
        )

    return command, note


# TS2CG workflow preview


def ts2cg_compatibility(
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
) -> Tuple[bool, List[str]]:
    """Check whether the current system matches the website's TS2CG interface."""
    reasons = []

    if resolution != "CG":
        reasons.append("The current TS2CG interface is limited to CG systems.")

    if membrane_type != "vesicle":
        reasons.append(
            "The current TS2CG interface is intended for vesicle systems."
        )

    if leaflet_mode != "ratio":
        reasons.append(
            "The current TS2CG composition template uses ratio-mode leaflets."
        )

    if not upper or not lower:
        reasons.append("Both vesicle leaflets must contain lipid composition.")

    return (len(reasons) == 0), reasons


def build_ts2cg_preview(
    upper: Dict[str, Any],
    lower: Dict[str, Any],
) -> Tuple[str, str, str]:
    """
    Build a scientifically honest TS2CG template.

    Official TS2CG input.str lipid lines use:
        LipidName  RatioUp  RatioDown  Area/Lipid

    LiMBS does not currently store lipid-specific TS2CG APL values, so the
    preview uses explicit placeholders instead of inventing numbers.
    """
    lipid_names = sorted(set(upper) | set(lower))

    lines = [
        "[Lipids List]",
        ";LipidName  RatioUp  RatioDown  Area/Lipid",
        "Domain 0",
    ]

    for lipid in lipid_names:
        ratio_up = format_decimal(float(upper.get(lipid, 0.0)))
        ratio_down = format_decimal(float(lower.get(lipid, 0.0)))

        lines.append(
            f"{lipid}  {ratio_up}  {ratio_down}  <APL_{lipid}>"
        )

    lines.append("End")

    composition_template = "\n".join(lines)

    command_template = (
        "TS2CG PLM "
        "-TSfile Sphere.tsi "
        "-bilayerThickness <builder_value> "
        "-rescalefactor <rx> <ry> <rz> "
        "-Mashno <builder_value>\n\n"
        "TS2CG PCG "
        "-dts point "
        "-str input.str "
        "-Bondlength <builder_value> "
        "-LLIB Martini3.LIB "
        "-defout vesicle"
    )

    note = (
        "TS2CG requires an area-per-lipid value for each lipid in input.str. "
        "Those values depend on the builder/model setup and are not part of "
        "the current LiMBS core specification, so this website shows explicit "
        "APL placeholders rather than inventing values. Likewise, TS2CG "
        "surface rescaling is left as a builder-specific parameter and is not "
        "derived automatically from the LiMBS simulation box."
    )

    return composition_template, command_template, note


# Workflow selection


def evaluate_builder(
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    upper: Dict[str, Any],
    lower: Dict[str, Any],
    box: List[float],
    solvent: str,
    salt_conc: float,
) -> Dict[str, Any]:
    """
    Return one simple automatic construction recommendation.

    AA/SCG results deliberately show no INSANE/TS2CG recommendation.
    CG planar systems recommend the LiMBS-INSANE interface.
    CG vesicles recommend the LiMBS-TS2CG interface.
    """
    if resolution != "CG":
        return {
            "show": False,
            "name": None,
            "status": None,
            "summary": None,
            "preview_title": None,
            "preview": None,
            "note": None,
            "script_key": None,
            "script_name": None,
        }

    if membrane_type == "planar":
        compatible, reasons = insane_compatibility(
            resolution=resolution,
            membrane_type=membrane_type,
            upper=upper,
            lower=lower,
        )

        return {
            "show": True,
            "name": "INSANE",
            "status": "Recommended" if compatible else "Check input",
            "summary": "Recommended construction interface for the selected CG planar membrane.",
            "preview_title": "LiMBS-INSANE usage",
            "preview": "python LiMBS_INSANE.py system_CG.txt --dry-run",
            "note": (
                "Download the LiMBS-INSANE interface script and place it with "
                "LiMBSv1_parser.py. The wrapper validates the LiMBS file and "
                "generates the corresponding INSANE command."
                if compatible
                else " ".join(reasons)
            ),
            "script_key": "insane",
            "script_name": "LiMBS_INSANE.py",
        }

    if membrane_type == "vesicle":
        compatible, reasons = ts2cg_compatibility(
            resolution=resolution,
            membrane_type=membrane_type,
            leaflet_mode=leaflet_mode,
            upper=upper,
            lower=lower,
        )

        lipid_names = sorted(set(upper) | set(lower))
        apl_args = " ".join(f"--apl {lipid}=<value>" for lipid in lipid_names)
        preview = (
            f"python LiMBS_TS2CG.py vesicle_CG.txt {apl_args} --dry-run"
            if apl_args
            else "python LiMBS_TS2CG.py vesicle_CG.txt --dry-run"
        )

        return {
            "show": True,
            "name": "TS2CG",
            "status": "Recommended" if compatible else "Requires ratio mode",
            "summary": "Recommended construction interface for the selected CG vesicle.",
            "preview_title": "LiMBS-TS2CG usage",
            "preview": preview,
            "note": (
                "Download the LiMBS-TS2CG interface script and place it with "
                "LiMBSv1_parser.py. The current adapter uses leaflet:ratio and "
                "requires an explicit APL value for each lipid when generating "
                "TS2CG input.str."
                if compatible
                else " ".join(reasons)
            ),
            "script_key": "ts2cg",
            "script_name": "LiMBS_TS2CG.py",
        }

    return {
        "show": False,
        "name": None,
        "status": None,
        "summary": None,
        "preview_title": None,
        "preview": None,
        "note": None,
        "script_key": None,
        "script_name": None,
    }


# Result-page rendering


def render_result_page(
    *,
    resolution: str,
    membrane_type: str,
    leaflet_mode: str,
    selected_lipids: List[str],
    box: List[float],
    solvent: str,
    salt_conc: float,
    salt_species: str,
    box_guidance: str,
    box_info: Dict[str, Any],
    notation: str,
    validation: Dict[str, Any],
    workflow: Dict[str, Any],
) -> str:
    """Render the generated LiMBS result page."""
    box_summary = (
        f"{format_box_value(box[0])} × "
        f"{format_box_value(box[1])} × "
        f"{format_box_value(box[2])} nm"
    )

    validation_status = "PASSED" if validation["passed"] else "FAILED"
    warnings = validation.get("warnings", [])

    if warnings:
        warning_text = "\n".join(
            f"[{item.get('code', 'WARNING')}] {item.get('message', item)}"
            for item in warnings
        )
    else:
        warning_text = "None"

    canonical_html = (
        f"<pre>{escape(validation['canonical'])}</pre>"
        if validation.get("canonical")
        else "<p>Canonical LiMBS unavailable because validation failed.</p>"
    )

    json_data = validation.get("json")
    json_html = (
        f"<pre>{escape(json.dumps(json_data, indent=2, ensure_ascii=False))}</pre>"
        if json_data is not None
        else "<p>Structured JSON unavailable because parsing failed.</p>"
    )

    apl_html = ""
    if box_info.get("mode") == "estimate" and membrane_type == "planar":
        upper_apl = box_info.get("upper_apl")
        lower_apl = box_info.get("lower_apl")
        apl_source = box_info.get("apl_source")

        if upper_apl is not None and lower_apl is not None:
            source_label = (
                "LiMBS reference / composition-weighted estimate"
                if apl_source == "limbs"
                else "User-specified APL"
            )
            apl_html = f"""
            <div class="section">
                <h2>Assisted APL Summary</h2>
                <table>
                    <tr><th>APL source</th><td>{escape(source_label)}</td></tr>
                    <tr><th>Upper leaflet</th><td>{upper_apl:.5f} nm²/lipid</td></tr>
                    <tr><th>Lower leaflet</th><td>{lower_apl:.5f} nm²/lipid</td></tr>
                </table>
                <div class="guidance-box">
                    Reference and composition-weighted APL values are starting
                    construction estimates, not equilibrium membrane properties.
                </div>
            </div>
            """

    builder_html = ""
    if workflow.get("show"):
        preview_html = ""
        if workflow.get("preview_title") and workflow.get("preview"):
            preview_html = (
                f"<h3>{escape(workflow['preview_title'])}</h3>"
                f"<pre>{escape(workflow['preview'])}</pre>"
            )

        download_html = ""
        if workflow.get("script_key") and workflow.get("script_name"):
            download_url = url_for(
                "download_builder_script",
                builder_name=workflow["script_key"],
            )
            download_html = f"""
                <p>
                    <a class="download-button" href="{escape(download_url)}">
                        Download {escape(workflow['script_name'])}
                    </a>
                </p>
            """

        builder_html = f"""
        <div class="section">
            <h2>Construction Tool</h2>
            <p><strong>{escape(workflow['name'])}</strong></p>
            <p><strong>Status:</strong> {escape(workflow['status'])}</p>
            <p>{escape(workflow['summary'])}</p>
            {download_html}
            {preview_html}
            <div class="note-box">
                {escape(workflow.get('note') or '')}
            </div>
        </div>
        """

    return f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Generated LiMBS1.0 Specification</title>
        <style>
            * {{ box-sizing: border-box; }}
            body {{
                font-family: Arial, Helvetica, sans-serif;
                margin: 0;
                color: #222;
                background: #fff;
                line-height: 1.5;
            }}
            .navbar {{
                background: #f5f5f5;
                border-bottom: 1px solid #ddd;
                padding: 14px 32px;
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            .brand {{ font-size: 20px; font-weight: bold; }}
            .nav-links a {{ margin-left: 20px; color: #1a5fb4; text-decoration: none; }}
            .container {{ max-width: 1250px; margin: 35px auto; padding: 0 25px 50px; }}
            .section {{ border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin-bottom: 24px; }}
            table {{ border-collapse: collapse; width: 100%; max-width: 760px; }}
            th, td {{ border: 1px solid #ddd; padding: 9px 12px; text-align: left; }}
            th {{ width: 220px; background: #f5f5f5; }}
            pre {{
                white-space: pre-wrap;
                word-break: break-word;
                overflow-x: auto;
                background: #f7f7f7;
                border: 1px solid #ddd;
                border-radius: 6px;
                padding: 15px;
                font-size: 13px;
            }}
            .pass {{ color: #176b2c; font-weight: bold; }}
            .fail {{ color: #b42318; font-weight: bold; }}
            .note-box {{ background: #fffde7; border-left: 4px solid #e0a800; padding: 13px 15px; margin-top: 15px; }}
            .guidance-box {{ background: #f7fbff; border-left: 4px solid #4a90e2; padding: 13px 15px; margin-top: 15px; }}
            a {{ color: #1a5fb4; }}
            .download-button {{
                display: inline-block;
                padding: 9px 14px;
                border: 1px solid #1a5fb4;
                border-radius: 6px;
                text-decoration: none;
                font-weight: 600;
            }}
        </style>
    </head>
    <body>
        <div class="navbar">
            <div class="brand">LiMBS1.0</div>
            <div class="nav-links">
                <a href="/">Build</a>
                <a href="/lipids">Lipid Library</a>
            </div>
        </div>

        <div class="container">
            <h1>Generated LiMBS1.0 Specification</h1>

            <div class="section">
                <h2>System Summary</h2>
                <table>
                    <tr><th>Resolution</th><td>{escape(resolution)}</td></tr>
                    <tr><th>Geometry</th><td>{escape(membrane_type)}</td></tr>
                    <tr><th>Leaflet mode</th><td>{escape(leaflet_mode)}</td></tr>
                    <tr><th>Lipid species</th><td>{len(selected_lipids)}</td></tr>
                    <tr><th>Box</th><td>{escape(box_summary)}</td></tr>
                    <tr><th>Solvent</th><td>{escape(solvent)}</td></tr>
                    <tr><th>Salt</th><td>{salt_conc:.2f} M {escape(salt_species)}</td></tr>
                </table>
            </div>

            <div class="section">
                <h2>Box-Dimension Guidance</h2>
                <div class="guidance-box">{escape(box_guidance)}</div>
            </div>

            {apl_html}

            <div class="section">
                <h2>Generated LiMBS</h2>
                <pre>{escape(notation)}</pre>
            </div>

            <div class="section">
                <h2>LiMBS Parser Validation</h2>
                <p>
                    <strong>Validation:</strong>
                    <span class="{'pass' if validation['passed'] else 'fail'}">{validation_status}</span>
                </p>
                <p>{escape(validation['message'])}</p>
                <p><strong>Parser warnings:</strong></p>
                <pre>{escape(warning_text)}</pre>
            </div>

            <div class="section">
                <h2>Canonical LiMBS</h2>
                {canonical_html}
            </div>

            <div class="section">
                <h2>Structured JSON</h2>
                {json_html}
            </div>

            {builder_html}

            <p><a href="/">Generate another LiMBS1.0 specification</a></p>
        </div>
    </body>
    </html>
    """


# Flask routes


BUILDER_SCRIPT_FILES = {
    "insane": {
        "download_name": "LiMBS_INSANE.py",
        "candidates": ("LiMBS_INSANE.py", "limbs_insane.py"),
    },
    "ts2cg": {
        "download_name": "LiMBS_TS2CG.py",
        "candidates": ("LiMBS_TS2CG.py", "limbs_ts2cg.py", "LiMBS_vesicle.py"),
    },
}


@app.route("/")
def home():
    return render_template(
        "form.html",
        lipids=LIPIDS,
        apl_reference_available=(get_reference_for_ui is not None),
        apl_availability=build_apl_availability_map(),
    )


@app.route("/lipids")
def lipid_library():
    return render_template("lipids.html", lipids=LIPIDS)


@app.route("/download/builder/<builder_name>")
def download_builder_script(builder_name: str):
    """Download the LiMBS builder-interface script from the application folder."""
    key = str(builder_name).strip().lower()
    config = BUILDER_SCRIPT_FILES.get(key)
    if config is None:
        abort(404)

    script_path = None
    for candidate in config["candidates"]:
        path = PROJECT_DIR / candidate
        if path.is_file():
            script_path = path
            break

    if script_path is None:
        return error_page(
            "Builder script unavailable",
            (
                f"Could not find {config['download_name']} in the LiMBS application "
                "directory. Place the interface script beside app.py."
            ),
            status_code=404,
        )

    return send_file(
        script_path,
        as_attachment=True,
        download_name=config["download_name"],
        mimetype="text/x-python",
    )


@app.route("/generate", methods=["POST"])
def generate():
    try:
        resolution = request.form.get("resolution", "CG").strip().upper()
        membrane_type = request.form.get("type", "planar").strip().lower()
        leaflet_mode = request.form.get("leaflet_mode", "count").strip().lower()

        if resolution not in {"AA", "CG", "SCG"}:
            raise ValueError("Resolution must be AA, CG, or SCG.")

        if membrane_type not in {"planar", "vesicle"}:
            raise ValueError("Geometry must be planar or vesicle.")

        if leaflet_mode not in {"count", "ratio"}:
            raise ValueError("Leaflet mode must be count or ratio.")

        selected_lipids, upper, lower = read_lipid_composition(
            leaflet_mode=leaflet_mode,
        )

        box, box_guidance, box_info = determine_box(
            resolution=resolution,
            membrane_type=membrane_type,
            leaflet_mode=leaflet_mode,
            upper=upper,
            lower=lower,
        )

        solvent = read_solvent(resolution)
        salt_conc, salt_species = read_salt()

        notation = build_limbs_notation(
            resolution=resolution,
            membrane_type=membrane_type,
            leaflet_mode=leaflet_mode,
            selected_lipids=selected_lipids,
            upper=upper,
            lower=lower,
            box=box,
            solvent=solvent,
            salt_conc=salt_conc,
            salt_species=salt_species,
        )

        validation = validate_generated_notation(notation)

        workflow = evaluate_builder(
            resolution=resolution,
            membrane_type=membrane_type,
            leaflet_mode=leaflet_mode,
            upper=upper,
            lower=lower,
            box=box,
            solvent=solvent,
            salt_conc=salt_conc,
        )

        return render_result_page(
            resolution=resolution,
            membrane_type=membrane_type,
            leaflet_mode=leaflet_mode,
            selected_lipids=selected_lipids,
            box=box,
            solvent=solvent,
            salt_conc=salt_conc,
            salt_species=salt_species,
            box_guidance=box_guidance,
            box_info=box_info,
            notation=notation,
            validation=validation,
            workflow=workflow,
        )

    except ValueError as exc:
        return error_page("Input error", str(exc), status_code=400)

    except Exception as exc:
        return error_page("Application error", str(exc), status_code=500)


# Development server


if __name__ == "__main__":
    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True,
        use_reloader=False,
    )
