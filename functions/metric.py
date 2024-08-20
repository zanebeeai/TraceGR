import sympy as sp


def getMetric(
    lineElement, coordSystem="Cartesian", subs=None, overrideConst=False
):  # the override lets the code run faster if you know for sure your line element will work out
    if coordSystem not in [
        "Cartesian",
        "PlanePolar",
        "SphericalPolar",
        "CylindricalPolar",
    ]:
        raise ValueError("Unknown coordinate system")

    lineElement = sp.expand(lineElement)
    coords = (t, x, y, z)

    dim = len(coords)
    g = sp.zeros(dim)

    for mu in range(dim):
        for nu in range(dim):
            coeff = lineElement.coeff(sp.diff(coords[mu]) * sp.diff(coords[nu]))
            if mu != nu and coeff != 0:
                g[mu, nu] = coeff.subs(subs) / 2
            else:
                g[mu, nu] = coeff.subs(subs)

    # Check for unexpected terms in the line element
    if not overrideConst:
        reconstructed_line_element = sum(
            g[i, j] * sp.diff(coords[i]) * sp.diff(coords[j])
            for i in range(dim)
            for j in range(dim)
        )
        if sp.simplify(lineElement.subs(subs) - reconstructed_line_element) != 0:
            raise ValueError(
                "Line element contains terms that are not pure differentials of the coordinates used"
            )
    return g
