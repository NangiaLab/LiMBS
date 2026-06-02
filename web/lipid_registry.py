LIPIDS = {
    # PA
    "DOPA": {"header": "DOPA:16:1-18:1, PA-head, charge -1.0", "class": "PA", "tail_a": (16, 1), "tail_b": (18, 1), "apl": 0.64},
    "POPA": {"header": "POPA:16:0-18:1, PA-head, charge -1.0", "class": "PA", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.64},
    "PIPA": {"header": "PIPA:16:0-18:2, PA-head, charge -1.0", "class": "PA", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.64},
    "PAPA": {"header": "PAPA:16:0-20:4, PA-head, charge -1.0", "class": "PA", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.64},

    # PC
    "POPC": {"header": "POPC:16:0-18:1, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.64},
    "PIPC": {"header": "PIPC:16:0-18:2, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.64},
    "PEPC": {"header": "PEPC:16:0-20:2, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (20, 2), "apl": 0.64},
    "PAPC": {"header": "PAPC:16:0-20:4, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.64},
    "PUPC": {"header": "PUPC:16:0-22:6, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (22, 6), "apl": 0.64},
    "UPC":  {"header": "UPC:20:5-22:6, PC-head, charge 0.0", "class": "PC", "tail_a": (20, 5), "tail_b": (22, 6), "apl": 0.64},
    "APC":  {"header": "APC:20:4-22:5, PC-head, charge 0.0", "class": "PC", "tail_a": (20, 4), "tail_b": (22, 5), "apl": 0.64},
    "IPC":  {"header": "IPC:16:2-18:0, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 2), "tail_b": (18, 0), "apl": 0.64},
    "OPC":  {"header": "OPC:16:1-18:1, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 1), "tail_b": (18, 1), "apl": 0.64},
    "PPC":  {"header": "PPC:16:0-18:0, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (18, 0), "apl": 0.64},
    "DPPC": {"header": "DPPC:16:0-sat, PC-head, charge 0.0", "class": "PC", "tail_a": (16, 0), "tail_b": (16, 0), "apl": 0.65},

    # PE
    "POPE": {"header": "POPE:16:0-18:1, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.62},
    "PIPE": {"header": "PIPE:16:0-18:2, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.62},
    "PQPE": {"header": "PQPE:16:0-20:3, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 0), "tail_b": (20, 3), "apl": 0.62},
    "PAPE": {"header": "PAPE:16:0-20:4, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.62},
    "DAPE": {"header": "DAPE:20:4-22:5, PE-head, charge 0.0", "class": "PE", "tail_a": (20, 4), "tail_b": (22, 5), "apl": 0.62},
    "PUPE": {"header": "PUPE:16:0-22:6, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 0), "tail_b": (22, 6), "apl": 0.62},
    "DUPE": {"header": "DUPE:20:5-22:6, PE-head, charge 0.0", "class": "PE", "tail_a": (20, 5), "tail_b": (22, 6), "apl": 0.62},
    "DOPE": {"header": "DOPE:16:1-18:1, PE-head, charge 0.0", "class": "PE", "tail_a": (16, 1), "tail_b": (18, 1), "apl": 0.62},

    # PS
    "POPS": {"header": "POPS:16:0-18:1, PS-head, charge -1.0", "class": "PS", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.65},
    "PIPS": {"header": "PIPS:16:0-18:2, PS-head, charge -1.0", "class": "PS", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.65},
    "PQPS": {"header": "PQPS:16:0-20:3, PS-head, charge -1.0", "class": "PS", "tail_a": (16, 0), "tail_b": (20, 3), "apl": 0.65},
    "PAPS": {"header": "PAPS:16:0-20:4, PS-head, charge -1.0", "class": "PS", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.65},
    "DAPS": {"header": "DAPS:20:4-22:5, PS-head, charge -1.0", "class": "PS", "tail_a": (20, 4), "tail_b": (22, 5), "apl": 0.65},
    "DUPS": {"header": "DUPS:20:5-22:6, PS-head, charge -1.0", "class": "PS", "tail_a": (20, 5), "tail_b": (22, 6), "apl": 0.65},
    "DOPS": {"header": "DOPS:18:1-18:1, PS-head, charge -1.0", "class": "PS", "tail_a": (18, 1), "tail_b": (18, 1), "apl": 0.65},

    # PI
    "POPI": {"header": "POPI:16:0-18:1, PI-head, charge -1.0", "class": "PI", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.70},
    "PIPI": {"header": "PIPI:16:0-18:2, PI-head, charge -1.0", "class": "PI", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.70},
    "PAPI": {"header": "PAPI:16:0-20:4, PI-head, charge -1.0", "class": "PI", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.70},
    "PUPI": {"header": "PUPI:16:0-22:6, PI-head, charge -1.0", "class": "PI", "tail_a": (16, 0), "tail_b": (22, 6), "apl": 0.70},

    # PG
    "DOPG": {"header": "DOPG:16:1-18:1, PG-head, charge -1.0", "class": "PG", "tail_a": (16, 1), "tail_b": (18, 1), "apl": 0.65},
    "POPG": {"header": "POPG:16:0-18:1, PG-head, charge -1.0", "class": "PG", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.65},

    # DG
    "PODG": {"header": "PODG:16:0-18:1, DG-head, charge 0.0", "class": "DG", "tail_a": (16, 0), "tail_b": (18, 1), "apl": 0.65},
    "PIDG": {"header": "PIDG:16:0-18:2, DG-head, charge 0.0", "class": "DG", "tail_a": (16, 0), "tail_b": (18, 2), "apl": 0.65},
    "PADG": {"header": "PADG:16:0-20:4, DG-head, charge 0.0", "class": "DG", "tail_a": (16, 0), "tail_b": (20, 4), "apl": 0.65},
    "PUDG": {"header": "PUDG:16:0-22:6, DG-head, charge 0.0", "class": "DG", "tail_a": (16, 0), "tail_b": (22, 6), "apl": 0.65},

    # SM and sterol
    "DPSM": {"header": "DPSM:18:1-18:0, SM-head, charge 0.0", "class": "SM", "tail_a": (18, 1), "tail_b": (18, 0), "apl": 0.65},
    "DXSM": {"header": "DXSM:24:1-24:0, SM-head, charge 0.0", "class": "SM", "tail_a": (24, 1), "tail_b": (24, 0), "apl": 0.65},
    "CHOL": {"header": "Cholesterol:sterol, OH-head, charge 0.0", "class": "STEROL", "tail_a": None, "tail_b": None, "apl": 0.40},
}
