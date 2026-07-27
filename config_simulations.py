# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Loader and structures for per-simulation configuration files

"""
This module centralizes all simulation configurations and structures.

Each simulation lives in its own file inside the 'simulations/' folder 
(e.g. simulations/example_1_simplebar.py), and each of those files exposes 
a module-level 'CONFIG' variable.
"""

import os
import importlib.util
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any

# ===============================================================================
# ======================== DATA STRUCTURES (DATA CLASSES) =======================
# ===============================================================================

@dataclass
class BoundaryCondition:
    """Define a boundary condition

    Attributes:
        name: identifier for the condition
        node_filter: lambda function to filter nodes (receives Node as argument)
                     OR a Gmsh physical-group tag (int/str)
        bc_type: type of BC (0: dirichlet x,y | 1: dirichlet x | 2: dirichlet y | 3: neumann)
        xval: displacement/force in x
        yval: displacement/force in y
    """
    name: str
    node_filter: callable  # callable or int/str
    bc_type: int
    xval: float
    yval: float

@dataclass
class MaterialParameters:
    """Material parameters

    Attributes:
        E: Young's modulus
        nu: Poisson's ratio
        Gc: Critical energy release rate
        l0: Length scale parameter
        material_type: type of material ('plane_stress', 'plane_strain')
        length: length of the bar (used for edge filters)
    """
    E: float
    nu: float
    Gc: float
    l0: float
    material_type: str = 'plane_stress'
    length: float = 1.0

@dataclass
class ReactionConfig:
    """Configuration for reaction force monitoring"""
    reaction_type: str
    reaction_dof: int
    sign_factor: float
    support_coords: Optional[Tuple[float, float]] = None

@dataclass
class GraphConfig:
    """Configuration for automatic plotting"""
    graph_type: str
    x_label: str
    y_label: str
    title: str
    output_file: str
    element_id: Optional[int] = None
    displacement_axis: str = 'y'

@dataclass
class SimulationConfig:
    """Global configuration configuration for a simulation run"""
    name: str
    mesh_file: Optional[str]
    mesh_type: str  # 'gmsh' or 'generated'
    material: MaterialParameters
    boundary_conditions: List[BoundaryCondition]
    simulation_params: Dict[str, Any]
    output_base: str
    graph_config: Optional[GraphConfig] = None
    reaction_config: Optional[ReactionConfig] = None
    imposed_displacement: float = 0.04
    mesh_generation_params: Optional[Dict[str, Any]] = None

# ===============================================================================
# =========================== NODE FILTER HELPERS ===============================
# ===============================================================================

def is_left_edge(length: float) -> Any:
    """Returns a function that identifies nodes on the left edge"""
    return lambda node: abs(node.x) < 1e-8

def is_right_edge(length: float) -> Any:
    """Returns a function that identifies nodes on the right edge"""
    return lambda node: abs(node.x - length) < 1e-8

def is_bottom_edge(height: float) -> Any:
    """Returns a function that identifies nodes on the bottom edge"""
    return lambda node: abs(node.y) < 1e-8

def is_top_edge(height: float) -> Any:
    """Returns a function that identifies nodes on the top edge"""
    return lambda node: abs(node.y - height) < 1e-8

def is_bottom_left_corner(length: float, height: float) -> Any:
    """Returns a function that identifies node on the bottom left corner"""
    return lambda node: (abs(node.x) < 1e-8 and abs(node.y) < 1e-8)

# ===============================================================================
# =================================== PATHS =====================================
# ===============================================================================

# Folder where each individual simulation configuration file is stored
SIMULATIONS_DIR = "simulations"

# ===============================================================================
# ============================= LOADER FUNCTIONS ================================
# ===============================================================================

def list_available_simulations() -> List[str]:
    """Scans the 'simulations/' folder and lists keys based on filenames.
    
    Example: 'simulations/tensile_test.py' becomes available as 'tensile_test'.
    """
    if not os.path.exists(SIMULATIONS_DIR):
        return []
    
    sim_names = []
    for f in os.listdir(SIMULATIONS_DIR):
        if f.endswith(".py") and not f.startswith("__"):
            sim_names.append(f[:-3])  # Remove the '.py' extension
    return sorted(sim_names)

def get_simulation_config(name: str) -> SimulationConfig:
    """Dynamically loads and returns the SimulationConfig from the simulation file.
    
    Looks for a file named simulations/{name}.py and reads its 'CONFIG' variable.
    """
    filepath = os.path.join(SIMULATIONS_DIR, f"{name}.py")
    
    if not os.path.exists(filepath):
        available = list_available_simulations()
        raise ValueError(
            f"Simulation profile '{name}' not found.\n"
            f"Path checked: '{filepath}'.\n"
            f"Available simulations in '{SIMULATIONS_DIR}/': {available}"
        )
        
    # Dynamic import magic via importlib
    spec = importlib.util.spec_from_file_location(f"dynamic_sim_{name}", filepath)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module spec for {filepath}")
        
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # Extract the CONFIG variable from the loaded module file
    if not hasattr(module, "CONFIG"):
        raise AttributeError(
            f"The simulation file '{filepath}' was loaded successfully, "
            f"but it is missing the required global 'CONFIG' variable."
        )
        
    config = getattr(module, "CONFIG")
    
    if not isinstance(config, SimulationConfig):
        raise TypeError(
            f"The global variable 'CONFIG' in '{filepath}' must be an instance "
            f"of 'SimulationConfig', but found '{type(config).__name__}' instead."
        )
        
    return config

# ===============================================================================
# =============================== EXAMPLE OF USE ================================
# ===============================================================================

if __name__ == "__main__":
    from config_simulations import get_simulation_config, list_available_simulations

    # Lists all available simulations (just scans filenames, does not import them)
    for name in list_available_simulations():
        print(name)
    # Load a specific configuration (only that file is imported)
    config = get_simulation_config('example_1_simplebar')
    print(f"Simulation parameters of {config.name}:\n {config.simulation_params}")