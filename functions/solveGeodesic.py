import numpy as np
import sympy as sp
import typing
import time


def intialPosition(
    t: sp.Symbol,
    x: sp.Symbol,
    y: sp.Symbol,
    z: sp.Symbol,
    px: sp.Symbol,
    py: sp.Symbol,
    pz: sp.Symbol,
    G: sp.Symbol,
    M: sp.Symbol,
    c: sp.Symbol,
    G_val: float,
    M_val: float,
    c_val: float,
    metric: float,
    pos: float,
    angles_len: int,
    pixSize: int,
    focal_length: int,
    resolution: typing.Tuple[int, int],
) -> np.ndarray:
    count = 0
    angles_len = 200
    angles = np.linspace(0, 2 * np.pi, angles_len)

    values = np.zeros(
        (resolution.prod(), 8)
    )  # 8 for t, x, y, z, px, py, pz, and an additional parameter if needed

    # Set initial positions (all rays start from the camera position)
    values[:, 0] = 0
    values[:, 1:4] = pos

    # Compute directions for each ray
    half_width = resolution[0] / 2
    half_height = resolution[1] / 2

    for i in range(resolution[0]):
        for j in range(resolution[1]):
            # Pixel center coordinates
            pixel_x = (i - half_width + 0.5) * pixSize
            pixel_y = (j - half_height + 0.5) * pixSize

            # Calculate unnormalized direction vector
            direction = np.array([pixel_x, pixel_y, focal_length])

            # Normalize the direction vector to have unit norm
            norm = np.linalg.norm(direction)
            print(norm)
            normalized_direction = direction / norm
            print(np.linalg.norm(normalized_direction))

            # Assign normalized direction to values
            values[i * resolution[1] + j, 5:8] = normalized_direction  # x', y', z'

    Bquad = 2 * (metric[0, 1] * px + metric[0, 2] * py + metric[0, 3] * pz)
    Cquad = (
        metric[1, 1] * px**2
        + metric[2, 2] * py**2
        + metric[3, 3] * pz**2
        + 2 * (metric[1, 2] * px * py + metric[1, 3] * px * pz + metric[2, 3] * py * pz)
    )

    pt1 = (-Bquad + sp.sqrt(Bquad**2 - 4 * metric[0, 0] * Cquad)) / (2 * metric[0, 0])
    pt2 = (-Bquad - sp.sqrt(Bquad**2 - 4 * metric[0, 0] * Cquad)) / (2 * metric[0, 0])

    pt1.subs([(t, -t)])

    pt1_lambda = sp.lambdify((t, x, y, z, G, M, c, px, py, pz), pt1, "numpy")
    pt2_lambda = sp.lambdify((t, x, y, z, G, M, c, px, py, pz), pt2, "numpy")
    args = (
        values[:, 0],
        values[:, 1],
        values[:, 2],
        values[:, 3],
        G_val,
        M_val,
        c_val,
        values[:, 5],
        values[:, 6],
        values[:, 7],
    )

    values[:, 4] = pt1_lambda(*args)

    return values


def solveGeodesic(
    values: np.ndarray,
    G: sp.Symbol,
    M: sp.Symbol,
    c: sp.Symbol,
    G_val: float,
    M_val: float,
    c_val: float,
    dptdl: sp.Expr,
    dpxdl: sp.Expr,
    dpydl: sp.Expr,
    dpzdl: sp.Expr,
    resolution: typing.Tuple[int, int],
    RENDER_RAD: float,
):
    lam = sp.symbols("lambda")
    # define initial conditions
    t = sp.Function("t")(lam)
    x = sp.Function("x")(lam)
    y = sp.Function("y")(lam)
    z = sp.Function("z")(lam)

    dt = sp.diff(t)
    dx = sp.diff(x)
    dy = sp.diff(y)
    dz = sp.diff(z)
    geods = []
    dL = 0.02
    meters = 800
    iterations = int(meters // dL)
    collision_info = []
    # append a numpy copy of the array to geod
    geods.append(np.copy(values))
    active_rays = np.ones(values.shape[0], dtype=bool)

    dptdl_lambda = sp.lambdify((t, x, y, z, dt, dx, dy, dz, G, M, c), dptdl, "numpy")
    dpxdl_lambda = sp.lambdify((t, x, y, z, dt, dx, dy, dz, G, M, c), dpxdl, "numpy")
    dpydl_lambda = sp.lambdify((t, x, y, z, dt, dx, dy, dz, G, M, c), dpydl, "numpy")
    dpzdl_lambda = sp.lambdify((t, x, y, z, dt, dx, dy, dz, G, M, c), dpzdl, "numpy")

    print("starting iterations")
    start = time.time()
    L = 0
    active_indices = np.where(active_rays)[
        0
    ]  # Initialize active indices before the loop

    progTotal = resolution.prod()

    for i in range(iterations):
        L += dL

        if active_indices.size == 0:
            break  # Early exit if no rays are active

        # Perform all operations only on active rays
        args = tuple(values[active_indices, k] for k in range(values.shape[1])) + (
            G_val,
            M_val,
            c_val,
        )

        dp = [
            dptdl_lambda(*args),
            dpxdl_lambda(*args),
            dpydl_lambda(*args),
            dpzdl_lambda(*args),
        ]
        dp_array = np.array(dp).T

        # Update velocities and positions
        values[active_indices, 4:] += dL * dp_array
        #     values[active_indices, 0:4] += dL * values[active_indices, 4:]
        values[active_indices, :4] += (
            dL
            * values[active_indices, 4:]
            / np.linalg.norm(values[active_indices, 4:], axis=1, keepdims=True)
        )

        # Check for collisions using a vectorized operation
        distances = np.linalg.norm(values[active_indices, 1:4], axis=1)
        collided = distances > RENDER_RAD

        # Deactivate collided rays in a vectorized way
        collision_indices = active_indices[collided]
        active_rays[collision_indices] = False
        for idx in collision_indices:
            pixel_idx = (idx // resolution[1], idx % resolution[1])
            collision_info.append((pixel_idx, values[idx, 1:4].tolist()))

        # Update active indices after processing collisions
        active_indices = np.where(active_rays)[0]

        # Progress update
        print("Percentage done: {:3.5f} %".format((i / iterations) * 100), end="\r")

    geods.append(np.copy(values))
    geods = np.array(geods)

    return geods
