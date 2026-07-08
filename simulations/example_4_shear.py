# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Simulation: Simple Shear Test

"""
This file defines a single simulation configuration in the module-level
variable CONFIG. config_simulations.get_simulation_config('example_4_shear')
loads exactly this file (and no other) at runtime.
"""

from config_simulations import SimulationConfig, MaterialParameters, BoundaryCondition, GraphConfig, ReactionConfig

CONFIG = SimulationConfig(
    name='Example 4: Shear Test',
    mesh_file='shear.msh',
    mesh_type='gmsh',
    material=MaterialParameters(E=210, nu=0.3, Gc=2.7e-3, l0=0.003, length=1.0, material_type='plane_strain'),
    boundary_conditions=[
        BoundaryCondition(
            name='top_id',
            node_filter=20,  # Gmsh tag (top_ids)
            bc_type=0,  # Dirichlet in x and y
            xval=0.04,  # imposed_displacement_x
            yval=0.0
        ),
        BoundaryCondition(
            name='bottom_id',
            node_filter=10,  # Gmsh tag (bottom_ids)
            bc_type=0,  # Full clamp
            xval=0.0,
            yval=0.0
        ),
    ],
    simulation_params={
        'dt': 0.005,
        'totaltime': 0.9,
        'maxsteps': int(1e5),
        'maxiter': 500,
        'stagtol': 1e-4,
    },
    output_base='outputs/ex4_shear',
    imposed_displacement=0.04,
    graph_config=GraphConfig(
        graph_type='force_vs_displacement',
        output_file='outputs/ex4_shear_force_vs_u.png',
        x_label='u (mm)',
        y_label='Force (kN)',
        title='Force vs imposed u',
        displacement_axis='x',
    ),
    reaction_config=ReactionConfig(
        reaction_type='bottom_ids_x',
        reaction_dof=0,  # X direction
        sign_factor=-1.0,  # negative sign
    ),
)
