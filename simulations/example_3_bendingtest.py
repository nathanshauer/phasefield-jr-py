# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Simulation: Bending Test

"""
This file defines a single simulation configuration in the module-level
variable CONFIG. config_simulations.get_simulation_config('example_3_bendingtest')
loads exactly this file (and no other) at runtime.
"""

from config_simulations import SimulationConfig, MaterialParameters, BoundaryCondition, GraphConfig, ReactionConfig

CONFIG = SimulationConfig(
    name='Example 3: Bending Test',
    mesh_file='bendingtest.msh',
    mesh_type='gmsh',
    material=MaterialParameters(E=20.8, nu=0.3, Gc=5.0e-4, l0=0.03, length=8.0, material_type='plane_stress'),
    boundary_conditions=[
        BoundaryCondition(
            name='left_id',
            node_filter=1,  # Gmsh tag (left_ids)
            bc_type=0,  # Full clamp
            xval=0.0,
            yval=0.0
        ),
        BoundaryCondition(
            name='right_id',
            node_filter=2,  # Gmsh tag (right_ids)
            bc_type=2,  # Dirichlet in y
            xval=0.0,
            yval=0.0
        ),
        BoundaryCondition(
            name='top_id',
            node_filter=3,  # Gmsh tag (top_ids)
            bc_type=2,  # Displacement in y
            xval=0.0,
            yval=-0.08  # imposed_displacement_y
        ),
    ],
    simulation_params={
        'dt': 0.01,
        'totaltime': 1.0,
        'maxsteps': int(1e5),
        'maxiter': 500,
        'stagtol': 1e-4,
    },
    output_base='outputs/ex3_bendingtest',
    imposed_displacement=0.01,
    graph_config=GraphConfig(
        graph_type='force_vs_displacement',
        output_file='outputs/ex3_bendingtest_force_vs_u.png',
        x_label='u (mm)',
        y_label='Force (kN)',
        title='Force vs imposed u',
        displacement_axis='y',
    ),
    reaction_config=ReactionConfig(
        reaction_type='supports_y',
        reaction_dof=1,  # Y direction
        sign_factor=1.0,  # positive sign
        support_coords=(-4.0, 0.0),  # left and right supports at y=0
    ),
)
