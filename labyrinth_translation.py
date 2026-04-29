import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import math

NUMBER_OF_TILES_IN_HEIGHT    = 15
NUMBER_OF_LABYRINTHS_TO_FIND = 1
SIZE_OF_LINE_IN_LABY         = 15

IMAGE_FILE = "lincoln.png"

target_image = np.array(Image.open(IMAGE_FILE).convert('L'))  # grayscale
height_target_orig, width_target_orig = target_image.shape
print(f"Image loaded: {width_target_orig}x{height_target_orig}")

nH = NUMBER_OF_TILES_IN_HEIGHT
nW = math.floor(width_target_orig / height_target_orig * (nH + 1))

height_target = (height_target_orig // nH) * nH
width_target  = (width_target_orig  // nW) * nW

target_image = np.array(
    Image.fromarray(target_image).resize((width_target, height_target), Image.LANCZOS)
)

pixels_in_height = height_target // nH
pixels_in_width  = width_target  // nW
print(f"Grid: {nH} rows x {nW} cols  ({pixels_in_height}x{pixels_in_width} px per tile)")

# ---------------------------------------------------------------------------
# Pixel averages over sub-blocks
# ---------------------------------------------------------------------------
pixel_average = np.zeros((nH, nW))
for i in range(nH):
    for j in range(nW):
        r0 = i * pixels_in_height
        c0 = j * pixels_in_width
        pixel_average[i, j] = target_image[r0:r0+pixels_in_height,
                                            c0:c0+pixels_in_width].mean()

# ---------------------------------------------------------------------------
# Tile type definitions (0-based)
#   0 = DR : open bottom + right
#   1 = DL : open bottom + left
#   2 = UR : open top    + right
#   3 = UL : open top    + left
#   4 = LR : open left   + right  (horizontal straight)
#   5 = UD : open top    + bottom (vertical straight)
#   6 = blank
# ---------------------------------------------------------------------------
NUMBER_OF_TOTAL_TILES = 7
n_vars = nH * nW * NUMBER_OF_TOTAL_TILES

TOP_OPEN    = [2, 3, 5]
BOTTOM_OPEN = [0, 1, 5]
LEFT_OPEN   = [1, 3, 4]
RIGHT_OPEN  = [0, 2, 4]

def var_idx(i, j, k):
    return (i * nW + j) * NUMBER_OF_TOTAL_TILES + k

# Reverse map: variable index -> (i, j, k)
index_to_grid = np.zeros((n_vars, 3), dtype=int)
for i in range(nH):
    for j in range(nW):
        for k in range(NUMBER_OF_TOTAL_TILES):
            index_to_grid[var_idx(i, j, k)] = [i, j, k]

# ---------------------------------------------------------------------------
# Objective: path tiles prefer bright pixels; blank tiles prefer dark pixels - CURRENTLY INVERSED
# ---------------------------------------------------------------------------
obj = np.zeros(n_vars)
for i in range(nH):
    for j in range(nW):
        pa = pixel_average[i, j]
        for k in range(NUMBER_OF_TOTAL_TILES - 1):
            obj[var_idx(i, j, k)] = (255 - pa) ** 2
        obj[var_idx(i, j, 6)] = pa ** 2 

# ---------------------------------------------------------------------------
# Equality constraints
# ---------------------------------------------------------------------------
eq_rows_data = []
eq_rhs       = []

def add_eq(row_dict, rhs):
    eq_rows_data.append(row_dict)
    eq_rhs.append(rhs)

# One tile per grid position
for i in range(nH):
    for j in range(nW):
        add_eq({var_idx(i, j, k): 1 for k in range(NUMBER_OF_TOTAL_TILES)}, 1)

# Horizontal matching: right-open count of (i,j) == left-open count of (i,j+1)
for i in range(nH):
    for j in range(nW - 1):
        row = {}
        for k in RIGHT_OPEN:
            row[var_idx(i, j,   k)] = row.get(var_idx(i, j,   k), 0) + 1
        for k in LEFT_OPEN:
            row[var_idx(i, j+1, k)] = row.get(var_idx(i, j+1, k), 0) - 1
        add_eq(row, 0)

# Vertical matching: bottom-open count of (i,j) == top-open count of (i+1,j)
for i in range(nH - 1):
    for j in range(nW):
        row = {}
        for k in BOTTOM_OPEN:
            row[var_idx(i,   j, k)] = row.get(var_idx(i,   j, k), 0) + 1
        for k in TOP_OPEN:
            row[var_idx(i+1, j, k)] = row.get(var_idx(i+1, j, k), 0) - 1
        add_eq(row, 0)

# ---------------------------------------------------------------------------
# Build Aeq matrix
# ---------------------------------------------------------------------------
n_eq = len(eq_rows_data)
Aeq  = np.zeros((n_eq, n_vars))
beq  = np.array(eq_rhs, dtype=float)
for ri, rd in enumerate(eq_rows_data):
    for col, val in rd.items():
        Aeq[ri, col] = val

# ---------------------------------------------------------------------------
# Variable bounds: border tiles that open outward are forbidden (ub = 0)
# ---------------------------------------------------------------------------
var_lb = np.zeros(n_vars)
var_ub = np.ones(n_vars)

for j in range(nW):
    for k in TOP_OPEN:    var_ub[var_idx(0,    j, k)] = 0
for j in range(nW):
    for k in BOTTOM_OPEN: var_ub[var_idx(nH-1, j, k)] = 0
for i in range(nH):
    for k in LEFT_OPEN:   var_ub[var_idx(i, 0,    k)] = 0
for i in range(nH):
    for k in RIGHT_OPEN:  var_ub[var_idx(i, nW-1, k)] = 0

# ---------------------------------------------------------------------------
# Tile traversal map for sub-tour detection
# TILE_EXITS[k][(entry_dr, entry_dc)] = (exit_dr, exit_dc)
# ---------------------------------------------------------------------------
TILE_EXITS = {
    0: {(-1,  0): ( 0,  1),  ( 0, -1): ( 1,  0)},  # DR
    1: {(-1,  0): ( 0, -1),  ( 0,  1): ( 1,  0)},  # DL
    2: {( 1,  0): ( 0,  1),  ( 0, -1): (-1,  0)},  # UR
    3: {( 1,  0): ( 0, -1),  ( 0,  1): (-1,  0)},  # UL
    4: {( 0,  1): ( 0, -1),  ( 0, -1): ( 0,  1)},  # LR
    5: {(-1,  0): ( 1,  0),  ( 1,  0): (-1,  0)},  # UD
}

def trace_loop(start_var, x_sol):
    """Trace the closed loop containing start_var. Returns list of var indices.

    Convention for TILE_EXITS:
        key   = the step direction taken to ARRIVE at the tile (exit_dir of prev tile)
        value = the step direction to EXIT the tile (toward next tile)
    """
    i0, j0, k0 = index_to_grid[start_var]
    for first_entry in TILE_EXITS[k0]:
        loop = [start_var]
        # 'entry' is the step we took to arrive at (ci,cj,ck)
        entry = first_entry
        ci, cj, ck = i0, j0, k0
        success = False
        for _ in range(n_vars):
            exit_dir = TILE_EXITS[ck][entry]          # step to take next
            ni, nj   = ci + exit_dir[0], cj + exit_dir[1]
            if not (0 <= ni < nH and 0 <= nj < nW):
                break
            # When we arrive at (ni,nj) we stepped exit_dir to get there,
            # so the entry key for the new tile is exit_dir itself.
            next_entry = exit_dir
            found_next = False
            for k in range(NUMBER_OF_TOTAL_TILES - 1):
                v = var_idx(ni, nj, k)
                if x_sol[v] > 0.5 and next_entry in TILE_EXITS[k]:
                    if v == start_var:
                        success = True
                        break
                    loop.append(v)
                    entry = next_entry
                    ci, cj, ck = ni, nj, k
                    found_next = True
                    break
            if success or not found_next:
                break
        if success:
            return loop
    return [start_var]

# ---------------------------------------------------------------------------
# ILP solver (Gurobi)
# ---------------------------------------------------------------------------
def solve_ilp(obj, Aeq, beq, var_lb, var_ub, A_ub_extra, b_ub_extra):
    m = gp.Model()
    m.setParam('OutputFlag', 0)  # suppress Gurobi console output

    # Binary decision variables with per-variable bounds
    x = m.addMVar(n_vars, lb=var_lb, ub=var_ub, vtype=GRB.BINARY, name="x")

    # Objective
    m.setObjective(obj @ x, GRB.MINIMIZE)

    # Equality constraints
    m.addMConstr(Aeq, x, '=', beq)

    # Inequality constraints (sub-tour elimination, grown dynamically)
    if A_ub_extra:
        A_ub = np.array(A_ub_extra)
        b_ub = np.array(b_ub_extra)
        m.addMConstr(A_ub, x, '<', b_ub)

    m.optimize()

    if m.Status != GRB.OPTIMAL:
        raise RuntimeError(f"Gurobi failed with status {m.Status}")

    return x.X

# ---------------------------------------------------------------------------
# Draw a labyrinth solution
# ---------------------------------------------------------------------------
def draw_labyrinth(x_sol, title="Labyrinth"):
    fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')
    ax.set_facecolor('black')
    ax.set_aspect('equal')
    ax.axis('off')
    fig.suptitle(title, color='white')

    for v in range(n_vars):
        if x_sol[v] < 0.5:
            continue
        i, j, k = index_to_grid[v]
        if k == 6:
            continue

        if k == 0:
            xp = [j+1,        (2*j+1)/2, (2*j+1)/2]
            yp = [-(2*i+1)/2, -(2*i+1)/2, -(i+1)]
        elif k == 1:
            xp = [j,          (2*j+1)/2, (2*j+1)/2]
            yp = [-(2*i+1)/2, -(2*i+1)/2, -(i+1)]
        elif k == 2:
            xp = [(2*j+1)/2, (2*j+1)/2, j+1]
            yp = [-i,        -(2*i+1)/2, -(2*i+1)/2]
        elif k == 3:
            xp = [(2*j+1)/2, (2*j+1)/2, j]
            yp = [-i,        -(2*i+1)/2, -(2*i+1)/2]
        elif k == 4:
            xp = [j,   j+1]
            yp = [-(2*i+1)/2, -(2*i+1)/2]
        elif k == 5:
            xp = [(2*j+1)/2, (2*j+1)/2]
            yp = [-i, -(i+1)]

        ax.plot(xp, yp, '-w', linewidth=SIZE_OF_LINE_IN_LABY)

    plt.tight_layout(pad=0)
    plt.show()
    print(f"Displayed: {title}")

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
A_ub_rows = []
b_ub_list = []

# Require at least 4 path tiles (smallest possible closed loop)
# This prevents the all-blank trivial solution
min_path_row = np.zeros(n_vars)
for i in range(nH):
    for j in range(nW):
        for k in range(NUMBER_OF_TOTAL_TILES - 1):  # exclude blank (k=6)
            min_path_row[var_idx(i, j, k)] = 1
A_ub_rows.append((-min_path_row).tolist())
b_ub_list.append(-4)

for laby_num in range(1, NUMBER_OF_LABYRINTHS_TO_FIND + 1):
    print(f"\n--- Searching for labyrinth {laby_num} ---")
    found_labyrinth = False

    while not found_labyrinth:
        x_sol = solve_ilp(obj, Aeq, beq, var_lb, var_ub, A_ub_rows, b_ub_list)

        active_vars = np.where(x_sol > 0.5)[0]
        blank_set   = set(active_vars[index_to_grid[active_vars, 2] == 6].tolist())
        path_vars   = [v for v in active_vars if v not in blank_set]
        n_path      = len(path_vars)
        print(f"  ILP solved: {n_path} path tiles, {len(blank_set)} blank tiles")

        all_found   = set(blank_set)
        any_subtour = False  # did we find at least one sub-tour this round?

        while True:
            remaining = [v for v in active_vars if v not in all_found]
            if not remaining:
                break
            loop = trace_loop(remaining[0], x_sol)

            # Always mark these tiles so we don't re-visit them in this inner loop
            all_found.update(loop)

            if len(loop) == n_path:
                print(f"  Single loop found ({n_path} tiles)!")
                found_labyrinth = True
                break

            print(f"  Sub-tour of length {len(loop)}, adding elimination constraint...")
            row = np.zeros(n_vars)
            for v in loop:
                row[v] = 1
            A_ub_rows.append(row.tolist())
            b_ub_list.append(len(loop) - 2)
            any_subtour = True

        # If we found sub-tours but no single loop, re-solve the ILP with new constraints
        if any_subtour and not found_labyrinth:
            continue

    draw_labyrinth(x_sol, title=f"Labyrinth {laby_num}")

    # Exclude this solution so the next labyrinth is different
    row = np.zeros(n_vars)
    for v in path_vars:
        row[v] = 1
    A_ub_rows.append(row.tolist())
    b_ub_list.append(len(path_vars) - 2)

print("\nDone.")