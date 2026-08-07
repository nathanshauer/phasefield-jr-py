# Simulation System - PhaseFieldJr

This document describes how to use the simulation configuration system for the PhaseFieldJr project.

## File Structure

- **`config_simulations.py`**: defines all data structures (`SimulationConfig`, `MaterialParameters`, etc.) and the loader functions. It does **not** contain any simulation data itself.
- **`simulations/<name>.py`**: each file defines exactly one simulation through a module-level `CONFIG` variable. Adding a new simulation means adding a new file here — `config_simulations.py` never needs to be edited.
- **`PhaseFieldJr.py`**: the phase-field solver. Its `main(config_name)` function loads the requested configuration, builds the mesh/material/BCs, and runs the staggered elasticity/phase-field scheme.

## How the System Works

When you run:

```bash
python PhaseFieldJr.py example_1_simplebar
```

`PhaseFieldJr.py` calls `get_simulation_config('example_1_simplebar')`, which uses `importlib` to load **only** `simulations/example_1_simplebar.py` and return its `CONFIG` variable. No other simulation file is ever imported or executed.

`list_available_simulations()` works by scanning filenames in the `simulations/` folder.

If no argument is given, `config_name` defaults to `'default'`, which does not exist, and the run fails with a clear error message listing what is available. Always pass an explicit simulation name.

To run any example, open a terminal, navigate to the project directory, and execute the following command:
```bash
python PhaseFieldJr.py <example_file>
```

## Important Implementation Notes

- **Only `mesh_type='gmsh'` is currently supported.** Any other value raises a `ValueError`. Make sure you have the Meshio library. If you do not, go to the terminal and run:
```bash
pip install meshio
```
- **BC values scale with pseudo-time.** `applyBoundaryConditions()` multiplies `bc.xval`/`bc.yval` by `pseudotime` on every step, so the value in the config is effectively a *rate*. With `totaltime=1.5` and `xval=0.08`, the imposed displacement keeps growing past `0.08` once `pseudotime > 1.0`.
- **`readGmshMesh()` returns `(nodes, elements, physical_groups)`**, where `physical_groups` maps each Gmsh physical tag (`int`) to a list of node IDs. `create_bc_nodes_from_config()` uses tags when `node_filter` is an `int`/`str`, or falls back to geometric filtering when it is a callable.
- **Thread count is system-aware.** `get_default_nthreads(reserve=2)` reads `os.cpu_count()` and reserves 2 cores for the OS by default. The value is computed once at import time as `DEFAULT_NTHREADS`. Inside `main()`, `mythread` is set from `config.simulation_params.get('nthreads', DEFAULT_NTHREADS)` — so you can override the thread count per simulation by adding `'nthreads': N` to `simulation_params`.


## Boundary Condition Types

| `bc_type` | Description | Typical use |
|---|---|---|
| 0 | Dirichlet in x and y | full clamp |
| 1 | Dirichlet in x | horizontal displacement |
| 2 | Dirichlet in y | vertical displacement |
| 3 | Neumann | concentrated load |

## Pre-configured Simulations

### `example_1_simplebar`
Bar under tension.
- Mesh: `simplebar.msh` · plane stress
- Material: E=30, ν=0.2, Gc=1.2e-4, l0=10.0, length=200.0
- BCs: tag `100` clamped (x,y); tag `200` Dirichlet-x `xval=0.08`
- `dt=0.01`, `totaltime=1.5`, `maxiter=500`, `stagtol=1e-4`
- Graph: `stress_vs_time` at `element_id=450` → `outputs/ex1_stress_vs_timeSimpleBar.png`

### `example_2_notchedplate`
Notched plate under tension.
- Mesh: `notchedplate.msh` · plane strain
- Material: E=210, ν=0.3, Gc=2.7e-3, l0=0.003, length=1.0
- BCs: tag `20` Dirichlet-y `yval=0.01` (top); tag `10` Dirichlet-y (base); tag `30` Dirichlet-x (center)
- `dt=0.01`, `totaltime=0.9`, `maxiter=500`, `stagtol=1e-4`
- Graph: `force_vs_displacement` (y-axis) → `outputs/ex2GEN_force_vs_u.png`

### `example_3_bendingtest`
Three-point bending test.
- Mesh: `bendingtest.msh` · plane stress
- Material: E=20.8, ν=0.3, Gc=5.0e-4, l0=0.03, length=8.0
- BCs: tag `1` full clamp (left); tag `2` Dirichlet-y (right); tag `3` Dirichlet-y `yval=-0.08` (load)
- `dt=0.01`, `totaltime=1.0`, `maxiter=500`, `stagtol=1e-4`
- Graph: `force_vs_displacement` (y-axis) → `outputs/Bending_force_vs_u.png`

### `example_4_shear`
Notched plate under shear.
- Mesh: `shear.msh` · plane strain
- Material: E=210, ν=0.3, Gc=2.7e-3, l0=0.003, length=1.0
- BCs: tag `20` Dirichlet x&y `xval=0.04` (top); tag `10` full clamp (base)
- `dt=0.005`, `totaltime=0.9`, `maxiter=500`, `stagtol=1e-4`
- Graph: `force_vs_displacement` (x-axis) → `outputs/ex4_force_vs_u.png`



## How to Add a New Simulation

Create a new file inside `simulations/` and set all the parameters and boundary conditions that you want. `config_simulations.py` does not need to be modified. 

Remember that you need to generate the .msh file for your simulation using Gmsh and put it inside `simulations/`.

```python
# simulations/my_simulation.py
from config_simulations import (
    SimulationConfig, MaterialParameters,
    BoundaryCondition, GraphConfig, ReactionConfig
)

CONFIG = SimulationConfig(
    name='My Simulation',
    mesh_file='my_mesh.msh',
    mesh_type='gmsh',
    material=MaterialParameters(
        E=30.0, nu=0.2, Gc=1.2e-4, l0=10.0, length=200.0,
        material_type='plane_stress'
    ),
    boundary_conditions=[
        BoundaryCondition(name='base', node_filter=10, bc_type=0, xval=0.0, yval=0.0),
        BoundaryCondition(name='top',  node_filter=20, bc_type=1, xval=0.04, yval=0.0),
    ],
    simulation_params={
        'dt': 0.01,
        'totaltime': 1.0,
        'maxsteps': int(1e5),
        'maxiter': 500,
        'stagtol': 1e-4,
        # 'nthreads': 4,  # optional: override the system-detected thread count
    },
    output_base='outputs/my_sim_',
    imposed_displacement=0.04,
    graph_config=GraphConfig(
        graph_type='force_vs_displacement',
        output_file='outputs/my_sim_force_vs_u.png',
        x_label='u (mm)', y_label='Force (kN)',
        title='Force vs imposed u',
        displacement_axis='x',
    ),
    reaction_config=ReactionConfig(
        reaction_type='bottom_ids_x',
        reaction_dof=0,
        sign_factor=-1.0,
    ),
)
```

Then run it with:

```bash
python PhaseFieldJr.py my_simulation
```

The new simulation will also appear automatically in `list_available_simulations()`.

## Utility Functions

You can get the parameters of any example:

```python
if __name__ == "__main__":
    from config_simulations import get_simulation_config, list_available_simulations

    # Lists all available simulations 
    for name in list_available_simulations():
        print(name)
    # Load a specific configuration 
    config = get_simulation_config('example_1_simplebar')
    print(f"Simulation parameters of {config.name}:\n {config.simulation_params}")
```
