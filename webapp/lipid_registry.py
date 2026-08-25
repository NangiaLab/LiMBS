"""
lipid_registry.py

Reusable lipid metadata for LiMBS 1.0.

Important:
- tail_a and tail_b preserve the ordering used by the corresponding
  supported Martini 3 topology.
- Do not reorder tails based only on the lipid name or header.
"""

from typing import Dict, Any


LIPIDS: Dict[str, Dict[str, Any]] = {'DOPA': {'header': 'DOPA:18:1-18:1, PA-head, charge -1.0',
          'class': 'PA',
          'tail_a': (18, 1),
          'tail_b': (18, 1),
          'charge': -1.0},
 'POPA': {'header': 'POPA:16:0-18:1, PA-head, charge -1.0',
          'class': 'PA',
          'tail_a': (16, 0),
          'tail_b': (18, 1),
          'charge': -1.0},
 'PIPA': {'header': 'PIPA:16:0-20:2, PA-head, charge -1.0',
          'class': 'PA',
          'tail_a': (16, 0),
          'tail_b': (20, 2),
          'charge': -1.0},
 'PAPA': {'header': 'PAPA:16:0-20:4, PA-head, charge -1.0',
          'class': 'PA',
          'tail_a': (16, 0),
          'tail_b': (20, 4),
          'charge': -1.0},
 'POPC': {'header': 'POPC:16:0-18:1, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (18, 1),
          'tail_b': (16, 0),
          'charge': 0.0},
 'PIPC': {'header': 'PIPC:16:0-20:2, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (20, 2),
          'tail_b': (16, 0),
          'charge': 0.0},
 'PEPC': {'header': 'PEPC:16:0-22:1, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (22, 1),
          'tail_b': (16, 0),
          'charge': 0.0},
 'PAPC': {'header': 'PAPC:16:0-20:4, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (20, 4),
          'tail_b': (16, 0),
          'charge': 0.0},
 'PUPC': {'header': 'PUPC:16:0-22:6, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (22, 6),
          'tail_b': (16, 0),
          'charge': 0.0},
 'DPPC': {'header': 'DPPC:16:0-16:0, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (16, 0),
          'tail_b': (16, 0),
          'charge': 0.0},
 'DOPC': {'header': 'DOPC:18:1-18:1, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (18, 1),
          'tail_b': (18, 1),
          'charge': 0.0},
 'DIPC': {'header': 'DIPC:18:2-18:2, PC-head, charge 0.0',
          'class': 'PC',
          'tail_a': (18, 2),
          'tail_b': (18, 2),
          'charge': 0.0},
 'POPE': {'header': 'POPE:16:0-18:1, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (18, 1),
          'tail_b': (16, 0),
          'charge': 0.0},
 'PIPE': {'header': 'PIPE:16:0-20:2, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (16, 0),
          'tail_b': (20, 2),
          'charge': 0.0},
 'PQPE': {'header': 'PQPE:16:0-20:3, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (16, 0),
          'tail_b': (20, 3),
          'charge': 0.0},
 'PAPE': {'header': 'PAPE:16:0-20:4, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (16, 0),
          'tail_b': (20, 4),
          'charge': 0.0},
 'DAPE': {'header': 'DAPE:20:4-20:4, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (20, 4),
          'tail_b': (20, 4),
          'charge': 0.0},
 'PUPE': {'header': 'PUPE:16:0-22:6, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (16, 0),
          'tail_b': (22, 6),
          'charge': 0.0},
 'DUPE': {'header': 'DUPE:12:0-12:0, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (12, 0),
          'tail_b': (12, 0),
          'charge': 0.0},
 'DOPE': {'header': 'DOPE:18:1-18:1, PE-head, charge 0.0',
          'class': 'PE',
          'tail_a': (18, 1),
          'tail_b': (18, 1),
          'charge': 0.0},
 'POPS': {'header': 'POPS:16:0-18:1, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (16, 0),
          'tail_b': (18, 1),
          'charge': -1.0},
 'PIPS': {'header': 'PIPS:16:0-20:2, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (16, 0),
          'tail_b': (20, 2),
          'charge': -1.0},
 'PQPS': {'header': 'PQPS:16:0-20:3, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (16, 0),
          'tail_b': (20, 3),
          'charge': -1.0},
 'PAPS': {'header': 'PAPS:16:0-20:4, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (16, 0),
          'tail_b': (20, 4),
          'charge': -1.0},
 'DAPS': {'header': 'DAPS:20:4-20:4, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (20, 4),
          'tail_b': (20, 4),
          'charge': -1.0},
 'DUPS': {'header': 'DUPS:12:0-12:0, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (12, 0),
          'tail_b': (12, 0),
          'charge': -1.0},
 'DOPS': {'header': 'DOPS:18:1-18:1, PS-head, charge -1.0',
          'class': 'PS',
          'tail_a': (18, 1),
          'tail_b': (18, 1),
          'charge': -1.0},
 'POPI': {'header': 'POPI:16:0-18:1, PI-head, charge -1.0',
          'class': 'PI',
          'tail_a': (16, 0),
          'tail_b': (18, 1),
          'charge': -1.0},
 'PIPI': {'header': 'PIPI:16:0-18:2, PI-head, charge -1.0',
          'class': 'PI',
          'tail_a': (16, 0),
          'tail_b': (18, 2),
          'charge': -1.0},
 'PAPI': {'header': 'PAPI:16:0-20:4, PI-head, charge -1.0',
          'class': 'PI',
          'tail_a': (16, 0),
          'tail_b': (20, 4),
          'charge': -1.0},
 'PUPI': {'header': 'PUPI:16:0-22:6, PI-head, charge -1.0',
          'class': 'PI',
          'tail_a': (16, 0),
          'tail_b': (22, 6),
          'charge': -1.0},
 'DOPG': {'header': 'DOPG:18:1-18:1, PG-head, charge -1.0',
          'class': 'PG',
          'tail_a': (18, 1),
          'tail_b': (18, 1),
          'charge': -1.0},
 'POPG': {'header': 'POPG:16:0-18:1, PG-head, charge -1.0',
          'class': 'PG',
          'tail_a': (16, 0),
          'tail_b': (18, 1),
          'charge': -1.0},
 'PODG': {'header': 'PODG:16:0-18:1, DG-head, charge 0.0',
          'class': 'DG',
          'tail_a': (16, 0),
          'tail_b': (18, 1),
          'charge': 0.0},
 'PIDG': {'header': 'PIDG:16:0-18:2, DG-head, charge 0.0',
          'class': 'DG',
          'tail_a': (16, 0),
          'tail_b': (18, 2),
          'charge': 0.0},
 'PADG': {'header': 'PADG:16:0-20:4, DG-head, charge 0.0',
          'class': 'DG',
          'tail_a': (16, 0),
          'tail_b': (20, 4),
          'charge': 0.0},
 'PUDG': {'header': 'PUDG:16:0-22:6, DG-head, charge 0.0',
          'class': 'DG',
          'tail_a': (16, 0),
          'tail_b': (22, 6),
          'charge': 0.0},
 'DPSM': {'header': 'DPSM:18:1-18:0, SM-head, charge 0.0',
          'class': 'SM',
          'tail_a': (18, 1),
          'tail_b': (18, 0),
          'charge': 0.0},
 'DXSM': {'header': 'DXSM:24:1-24:0, SM-head, charge 0.0',
          'class': 'SM',
          'tail_a': (24, 1),
          'tail_b': (24, 0),
          'charge': 0.0},
 'PSM': {'header': 'PSM:18:1-16:0, SM-head, charge 0.0',
         'class': 'SM',
         'tail_a': (18, 1),
         'tail_b': (16,0),
         'charge': 0.0},
 'NSM': {'header': 'NSM:18:1-24:1, SM-head, charge 0.0',
         'class': 'SM',
         'tail_a': (18, 1),
         'tail_b': (24, 1),
         'charge': 0.0},
 'XSM': {'header': 'XSM:18:1-24:1, SM-head, charge 0.0',
         'class': 'SM',
         'tail_a': (18, 1),
         'tail_b': (24, 1),
         'charge': 0.0},
 'OSM': {'header': 'OSM:18:1-16:1, SM-head, charge 0.0',
         'class': 'SM',
         'tail_a': (18, 1),
         'tail_b': (16, 1),
         'charge': 0.0},
 'USM': {'header': 'USM:18:1-18:0, SM-head, charge 0.0',
         'class': 'SM',
         'tail_a': (18, 1),
         'tail_b': (18, 0),
         'charge': 0.0},
 'CHOL': {'header': 'Cholesterol:sterol, OH-head, charge 0.0',
          'class': 'STEROL',
          'tail_a': None,
          'tail_b': None,
          'charge': 0.0}}


REQUIRED_KEYS = [
    "header",
    "class",
    "tail_a",
    "tail_b",
    "charge",
]


CLASS_DEFAULT_CHARGE = {
    "PC": 0.0,
    "PE": 0.0,
    "PA": -1.0,
    "PS": -1.0,
    "PI": -1.0,
    "PG": -1.0,
    "DG": 0.0,
    "SM": 0.0,
    "STEROL": 0.0,
    "CE": 0.0,
}


def get_lipid(name: str) -> Dict[str, Any]:
    """Return the registry entry for a lipid name (case-insensitive)."""
    key = name.upper()
    if key not in LIPIDS:
        raise KeyError(f"Unknown lipid: {name}")
    return LIPIDS[key]


def lipid_exists(name: str) -> bool:
    """Return True when a lipid is present in the predefined registry."""
    return name.upper() in LIPIDS


def get_charge(name: str) -> float:
    """Return the formal charge stored for a predefined lipid."""
    lipid = get_lipid(name)
    if "charge" in lipid and lipid["charge"] is not None:
        return float(lipid["charge"])

    lipid_class = lipid.get("class")
    if lipid_class in CLASS_DEFAULT_CHARGE:
        return CLASS_DEFAULT_CHARGE[lipid_class]

    raise ValueError(f"Charge unknown for lipid: {name}")


def normalize_lipid_entry(name: str, lipid: Dict[str, Any]) -> Dict[str, Any]:
    """Add convenience metadata without changing authoritative lipid metadata."""
    lipid.setdefault("name", name)
    lipid.setdefault("headgroup", lipid.get("class"))
    lipid.setdefault("charge", CLASS_DEFAULT_CHARGE.get(lipid.get("class")))
    lipid.setdefault("tail_order", "martini3/internal")
    return lipid


def validate_lipid_registry() -> None:
    """Validate the structure and required metadata of the predefined registry."""
    errors = []

    for name, lipid in LIPIDS.items():
        normalize_lipid_entry(name, lipid)

        for key in REQUIRED_KEYS:
            if key not in lipid:
                errors.append(f"{name}: missing required key '{key}'")

        if lipid.get("charge") is None:
            errors.append(f"{name}: charge is unknown")
        elif not isinstance(lipid.get("charge"), (float, int)):
            errors.append(f"{name}: charge must be numeric")

        lipid_class = lipid.get("class")
        if lipid_class not in CLASS_DEFAULT_CHARGE:
            errors.append(f"{name}: unknown lipid class '{lipid_class}'")

        for tail_key in ("tail_a", "tail_b"):
            tail = lipid.get(tail_key)
            if tail is not None:
                if (
                    not isinstance(tail, tuple)
                    or len(tail) != 2
                    or not all(isinstance(value, int) for value in tail)
                ):
                    errors.append(
                        f"{name}: {tail_key} must be None or a "
                        "(carbon_count, unsaturation_count) tuple"
                    )

    if errors:
        raise ValueError(
            "lipid_registry.py validation failed:\n" + "\n".join(errors)
        )


validate_lipid_registry()


if __name__ == "__main__":
    print(f"Validated {len(LIPIDS)} predefined lipids.")
