# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Simulation: Notched Plate under Tension

from config_simulations import SimulationConfig, MaterialParameters, BoundaryCondition, GraphConfig, ReactionConfig

CONFIG = SimulationConfig(
    name='Example 2: Notched Plate under tension',
    mesh_file='notchedplate.msh',
    mesh_type='gmsh',
    material=MaterialParameters(E=210, nu=0.3, Gc=2.7e-3, l0=0.003, material_type='plane_strain'),
    boundary_conditions=[
        BoundaryCondition(
            name='top_id',
            node_filter='Top_ids',  # Gmsh physical group name
            bc_type=2,  # Dirichlet in y only
            xval=0.0,
            yval=0.01  # imposed_displacement_y
        ),
        BoundaryCondition(
            name='bottom_id',
            node_filter='Bottom_ids',  # Gmsh physical group name
            bc_type=2,  # Dirichlet in y
            xval=0.0,
            yval=0.0
        ),
        BoundaryCondition(
            name='fixed_x_id',
            node_filter='Fixed_x',  # Gmsh physical group name
            bc_type=1,  # Dirichlet in x
            xval=0.0,
            yval=0.0
        ),
    ],
    simulation_params={
        'dt': 0.01,
        'totaltime': 0.9,
        'maxsteps': int(1e5),
        'maxiter': 500,
        'stagtol': 1e-4,
    },
    output_base='outputs/ex2_notchedplate',
    imposed_displacement=0.01,
    graph_config=GraphConfig(
        graph_type='force_vs_displacement',
        output_file='outputs/ex2_notchedplate_force_vs_u.png',
        x_label='u (mm)',
        y_label='Force (kN)',
        title='Force vs imposed u',
        displacement_axis='y',
    ),
    reaction_config=ReactionConfig(
        reaction_type='Bottom_ids',
        reaction_dof=1,  # Y direction
        sign_factor=-1.0,  # negative sign
    ),
)
