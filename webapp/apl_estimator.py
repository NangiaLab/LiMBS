# apl_estimator.py

from apl_reference_combined import get_reference_for_ui


def get_mixture_apl_estimate(
    composition,
    resolution,
    model,
    target_temperature_K=None,
):
    """
    Calculate a composition-weighted starting APL estimate.

    Parameters
    ----------
    composition : dict
        Example:
        {
            "POPC": 50,
            "POPE": 30,
            "POPS": 20,
        }

        Values may be counts or fractions. They are normalized internally.

    resolution : str
        "CG" or "AA"

    model : str
        e.g. "Martini3", "CHARMM36", or "experimental"

    Returns
    -------
    dict
        Information required by the web interface.
    """

    clean = {
        str(lipid).strip().upper(): float(amount)
        for lipid, amount in composition.items()
        if float(amount) > 0
    }

    if not clean:
        return {
            "available": False,
            "reason": "No lipids were supplied."
        }

    total = sum(clean.values())

    components = []
    missing = []

    for lipid, amount in clean.items():

        ref = get_reference_for_ui(
            lipid=lipid,
            resolution=resolution,
            model=model,
            target_temperature_K=target_temperature_K,
        )

        if ref is None:
            missing.append(lipid)
            continue

        fraction = amount / total

        components.append({
            "lipid": lipid,
            "amount": amount,
            "fraction": fraction,
            "apl_nm2": ref["value_nm2"],
            "temperature_K": ref.get("temperature_K"),
            "source": ref.get("source"),
        })

    # Do not silently estimate when any component lacks
    # a suitable reference.
    if missing:
        return {
            "available": False,
            "missing_lipids": missing,
            "reason": (
                "A composition-weighted estimate cannot be calculated "
                "because one or more components do not have a suitable "
                "LiMBS reference APL."
            ),
        }

    apl_estimate = sum(
        item["fraction"] * item["apl_nm2"]
        for item in components
    )

    return {
        "available": True,
        "value_nm2": apl_estimate,
        "components": components,
        "method": "composition-weighted pure-lipid reference estimate",
        "warning": (
            "This value is a composition-weighted construction estimate "
            "derived from pure-lipid reference systems. It is not an "
            "equilibrium APL for the mixed membrane."
        ),
    }
