# MIT License
# © 2025 Nathan Shauer
# phasefield-jr

"""Configuration structures and loader for PhaseFieldJr simulations.

Each simulation is defined in ``simulations/<name>.py`` and exposes a module-
level ``CONFIG`` object. Materials use the unified orthotropic constants
``E1``, ``E2``, ``nu12`` and ``G12`` plus ``theta_deg``. FFracture toughness is defined 
by ``Gc1`` and ``Gc2`` for crack advance along material axes 1 and 2, respectively.
"""

from dataclasses import dataclass
import importlib.util
import math
import os
from typing import Any, Dict, List, Optional, Tuple, Union


# =============================================================================
# CONFIGURATION DATA
# =============================================================================

@dataclass
class BoundaryCondition:
  """Boundary condition applied to a node set.

  ``node_filter`` is either a numeric Gmsh physical-group tag or a physical-
  group name from the mesh.

  ``bc_type`` values are 0 (x and y Dirichlet), 1 (x Dirichlet), 2 (y
  Dirichlet), and 3 (x/y nodal Neumann load).

  ``distribute_total`` only applies to ``bc_type=3`` (Neumann). When True,
  ``xval``/``yval`` are the TOTAL force to apply across every node matched by
  ``node_filter``; the solver splits it evenly (xval/n, yval/n) so the user
  does not need to know or count how many nodes lie on that boundary. When
  False (default), ``xval``/``yval`` are already the per-node nodal force.
  """

  name: str
  node_filter: Union[str, int]
  bc_type: int
  xval: float
  yval: float
  distribute_total: bool = False

  def __post_init__(self):
    if not isinstance(self.node_filter, (str, int)):
      raise TypeError(
        f"BoundaryCondition '{self.name}' node_filter must be a Gmsh "
        "physical-group name (str) or tag (int)."
      )
    if self.bc_type not in (0, 1, 2, 3):
      raise ValueError(
        f"BoundaryCondition '{self.name}' has invalid bc_type "
        f"{self.bc_type}; use 0, 1, 2 or 3."
      )
    if self.distribute_total and self.bc_type != 3:
      raise ValueError(
        f"BoundaryCondition '{self.name}': distribute_total is only "
        "meaningful for bc_type=3 (Neumann)."
      )


@dataclass
class MaterialParameters:
  """Elastic and phase-field material properties.

  All materials use the orthotropic form::

    MaterialParameters(E1=..., E2=..., nu12=..., G12=...,
                       Gc1=..., Gc2=..., l0=...,
                       material_type='plane_stress', length=...,
                       theta_deg=...)

  An isotropic material is represented by ``E1 == E2`` and
  ``G12 == E1 / (2 * (1 + nu12))``.

    The material direction 1 forms ``theta_deg`` degrees with the global x-axis.
  ``Gc1`` and ``Gc2`` are the effective toughnesses for cracks propagating
  parallel to material directions 1 and 2, respectively. Equal values
  recover isotropic fracture.
  """

  E1: float
  E2: float
  nu12: float
  G12: float
  Gc1: float
  Gc2: float
  l0: float
  material_type: str = "plane_stress"
  length: float = 1.0
  theta_deg: float = 0.0

  def __post_init__(self):
    self._validate()

  def _validate(self):
    for name in ("E1", "E2", "G12", "l0", "length"):
      value = getattr(self, name)
      if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number.")
    for name in ("Gc1", "Gc2"):
      value = getattr(self, name)
      if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number.")
    if not math.isfinite(self.nu12):
      raise ValueError("nu12 must be finite.")
    if not math.isfinite(self.theta_deg):
      raise ValueError("theta_deg must be finite.")
    if self.material_type not in ("plane_stress", "plane_strain"):
      raise ValueError(
        f"Unknown material_type '{self.material_type}'. "
        "Use 'plane_stress' or 'plane_strain'."
      )

    nu21 = self.nu12 * self.E2 / self.E1
    if 1.0 - self.nu12 * nu21 <= 0.0:
      raise ValueError(
        "Invalid orthotropic constants: 1 - nu12*nu21 must be positive."
      )

@dataclass
class ReactionConfig:
  """Configuration for reaction-force monitoring."""

  reaction_type: str
  reaction_dof: int
  sign_factor: float
  support_coords: Optional[Tuple[float, float]] = None

  def __post_init__(self):
    if self.reaction_dof not in (0, 1):
      raise ValueError("reaction_dof must be 0 (x) or 1 (y).")


@dataclass
class DisplacementMonitorConfig:
  """Surface displacement to monitor in a force-controlled simulation."""

  node_filter: Optional[Union[str, int]] = None
  dof: int = 1
  aggregation: str = "average"

  def __post_init__(self):
    if self.node_filter is not None and not isinstance(
        self.node_filter,
        (str, int),
    ):
      raise TypeError(
        "Displacement monitor node_filter must be a Gmsh physical-group "
        "name (str) or tag (int)."
      )
    if self.dof not in (0, 1):
      raise ValueError("Displacement monitor dof must be 0 (x) or 1 (y).")
    if self.aggregation not in ("average", "max"):
      raise ValueError(
        "Displacement monitor aggregation must be 'average' or 'max'."
      )


@dataclass
class GraphConfig:
  """Configuration for automatic result plotting."""

  graph_type: str
  x_label: str
  y_label: str
  title: str
  output_file: str
  element_id: Optional[int] = None
  displacement_axis: str = "y"
  # Used by graph_type='stress_vs_time': 0 = sigma_x, 1 = sigma_y, 2 = tau_xy.
  stress_component: int = 0

  def __post_init__(self):
    if self.displacement_axis not in ("x", "y"):
      raise ValueError("displacement_axis must be 'x' or 'y'.")
    if self.stress_component not in (0, 1, 2):
      raise ValueError(
        "stress_component must be 0 (sigma_x), 1 (sigma_y) or 2 (tau_xy)."
      )


@dataclass
class SimulationConfig:
  """Complete configuration for one simulation run."""

  name: str
  mesh_file: Optional[str]
  mesh_type: str
  material: MaterialParameters
  boundary_conditions: List[BoundaryCondition]
  simulation_params: Dict[str, Any]
  output_base: str
  graph_config: Optional[GraphConfig] = None
  reaction_config: Optional[ReactionConfig] = None
  imposed_displacement: float = 0.04
  mesh_generation_params: Optional[Dict[str, Any]] = None
  # "auto" infers the control type from the non-zero boundary conditions.
  loading_type: str = "auto"
  displacement_monitor: Optional[DisplacementMonitorConfig] = None

  def __post_init__(self):
    if self.loading_type not in ("auto", "displacement", "force"):
      raise ValueError(
        f"Unknown loading_type '{self.loading_type}'; use "
        "'auto', 'displacement' or 'force'."
      )


# =============================================================================
# SIMULATION LOADER
# =============================================================================

SIMULATIONS_DIR = "simulations"


def list_available_simulations() -> List[str]:
  """Return sorted names of Python configurations in ``simulations/``."""
  if not os.path.isdir(SIMULATIONS_DIR):
    return []
  return sorted(
    entry.name[:-3]
    for entry in os.scandir(SIMULATIONS_DIR)
    if (
      entry.is_file()
      and entry.name.endswith(".py")
      and not entry.name.startswith("__")
    )
  )


def get_simulation_config(name: str) -> SimulationConfig:
  """Load ``simulations/<name>.py`` and return its ``CONFIG`` object."""
  filepath = os.path.join(SIMULATIONS_DIR, f"{name}.py")
  if not os.path.exists(filepath):
    raise ValueError(
      f"Simulation profile '{name}' not found.\n"
      f"Path checked: '{filepath}'.\n"
      f"Available simulations: {list_available_simulations()}"
    )

  spec = importlib.util.spec_from_file_location(
    f"dynamic_sim_{name}",
    filepath,
  )
  if spec is None or spec.loader is None:
    raise ImportError(f"Could not load module spec for '{filepath}'.")

  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  if not hasattr(module, "CONFIG"):
    raise AttributeError(
      f"Simulation file '{filepath}' does not define global CONFIG."
    )

  config = module.CONFIG
  if not isinstance(config, SimulationConfig):
    raise TypeError(
      f"CONFIG in '{filepath}' must be a SimulationConfig, found "
      f"{type(config).__name__}."
    )
  return config


if __name__ == "__main__":
  for simulation_name in list_available_simulations():
    print(simulation_name)
