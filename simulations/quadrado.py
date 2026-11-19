# MIT License
# © 2025 Nathan Shauer
# phasefield-jr
# Simulation: Square under vertical tension

"""
This file defines a single simulation configuration in the module-level
variable CONFIG. config_simulations.get_simulation_config('example_1_simplebar')
loads exactly this file (and no other) at runtime.
"""

from config_simulations import SimulationConfig, MaterialParameters, BoundaryCondition, GraphConfig

CONFIG = SimulationConfig(
    name='Quadrado',
    mesh_file='quadrado.msh',
    mesh_type='gmsh',
    # Material equivalente ao isotropico original (E=15, nu=0.2), agora
    # descrito pelas constantes ortotropicas: E1=E2=E, nu12=nu, G12=E/2(1+nu).
    # Para um material ortotropico "de verdade", basta usar E1 != E2 e/ou
    # theta_deg != 0 (direcao 1 nao alinhada com o eixo x global).
    material=MaterialParameters(E1=14000, E2=1000, nu12=0.3,
                                 G12=500, theta_deg=90.0,
                                 Gc1=1.2e9, Gc2=1.2e9,
                                 l0=10.0, length=10.0,
                                 material_type='plane_stress'),
    boundary_conditions=[
        BoundaryCondition(
            name='bottom_id',
            node_filter=3,  # Gmsh tag: Bottom_ids
            bc_type=0,  # Dirichlet in x and y
            xval=0.0,
            yval=0.0
        ),
        BoundaryCondition(
            name='top_force',
            node_filter=4,  # Gmsh tag: Top_ids
            bc_type=3,  # Neumann: total force in x/y
            xval=0.0,
            yval=10.0,  # Total force; the solver counts the nodes
            distribute_total=True,
        ),
    ],
    simulation_params={
        'dt': 0.5,
        'totaltime': 1.0,
        'maxsteps': int(1e5),
        'maxiter': 600,
        'stagtol': 1e-4,
    },
    output_base='outputs/quadrado',
    imposed_displacement=1.0,
    graph_config=GraphConfig(
        graph_type='stress_vs_time',
        output_file='outputs/ex1_stress_vs_timeQuadrado.png',
        x_label='Pseudo Time',
        y_label='Stress/Stress_peak',
        title='Stress vs Pseudo Time',
        element_id=450,
    ),
)
