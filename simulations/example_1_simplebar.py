# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Simulation: Simple Bar under Tension

from config_simulations import SimulationConfig, MaterialParameters, BoundaryCondition, GraphConfig

CONFIG = SimulationConfig(
    name='Example 1: Simple Bar under tension',
    mesh_file='simplebar.msh',
    mesh_type='gmsh',
    material=MaterialParameters(E=30, nu=0.2, Gc=1.2e-4, l0=10.0, length=200.0, material_type='plane_stress'),
    boundary_conditions=[
        BoundaryCondition(
            name='left_id',
            node_filter=100,  # Gmsh tag
            bc_type=0,  # Dirichlet in x and y
            xval=0.0,
            yval=0.0
        ),
        BoundaryCondition(
            name='right_id',
            node_filter=200,  # Gmsh tag
            bc_type=1,  # Dirichlet in x only
            xval=0.08,  # u_peak_at2 (approximate)
            yval=0.0
        ),
    ],
    simulation_params={
        'dt': 0.01,
        'totaltime': 1.5,
        'maxsteps': int(1e5),
        'maxiter': 500,
        'stagtol': 1e-4,
    },
    output_base='outputs/ex1_simplebar',
    imposed_displacement=0.08,
    graph_config=GraphConfig(
        graph_type='stress_vs_time',
        output_file='outputs/ex1_stress_vs_timeSimpleBar.png',
        x_label='Pseudo Time',
        y_label='Stress/Stress_peak',
        title='Stress vs Pseudo Time',
        element_id=450,
    ),
)
