# MIT License
# © 2025 Nathan Shauer
# phasefield-jr

import numpy as np
import time
import matplotlib.pyplot as plt
import os
import sys
from scipy import sparse
import meshio  
from concurrent.futures import ProcessPoolExecutor

# =============================== DATA STRUCTURES ===============================
# ===============================================================================
class QuadraturePoint:
  def __init__(self, xi, eta, weight):
    self.xi = xi
    self.eta = eta
    self.weight = weight

class MaterialParameters:
  def __init__(self, E, nu, Gc, l0, material_type='plane_stress'):
    self.E = E  # Young's modulus
    self.nu = nu  # Poisson's ratio
    self.Gc = Gc  # Critical strain energy release rate
    self.l0 = l0  # Length scale parameter
    self.material_type = material_type  # 'plane_stress' ou 'plane_strain'
  
  def get_factor(self):
    """Calculates the factor for the D matrix based on the material type"""
    if self.material_type == 'plane_stress':
      
      return self.E / (1 - self.nu * self.nu)
    elif self.material_type == 'plane_strain':
      
      return self.E / ((1 + self.nu) * (1 - 2 * self.nu))
    else:
      # Default: plane stress
      return self.E / (1 - self.nu * self.nu)
  
  def get_D_matrix(self):
    """Calculate the 3x3 elasticity matrix D based on the material type."""
    factor = self.get_factor()
    D = np.zeros((3, 3))
    
    if self.material_type == 'plane_stress':
      
      D[0, 0] = factor
      D[0, 1] = factor * self.nu
      D[1, 0] = factor * self.nu
      D[1, 1] = factor
      D[2, 2] = factor * (1 - self.nu) / 2.0
    elif self.material_type == 'plane_strain':
      
      D[0, 0] = factor * (1 - self.nu)
      D[0, 1] = factor * self.nu
      D[1, 0] = factor * self.nu
      D[1, 1] = factor * (1 - self.nu)
      D[2, 2] = factor * (1 - 2 * self.nu) / 2.0
    
    return D

def create2x2QuadratureRule():
  points = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
  weights = [1.0, 1.0]
  rule = [QuadraturePoint(xi, eta, w1 * w2) for xi, w1 in zip(points, weights) for eta, w2 in zip(points, weights)]
  return rule

class Element:
  def __init__(self, node_ids):
    self.node_ids = node_ids

class Node:
  def __init__(self, x, y):
    self.x = x
    self.y = y

class BC:
  def __init__(self, node, bc_type, xval, yval):
    self.node = node
    self.bctype = bc_type  # 0 dirichlet in x and y, 1 dirichlet in x, 2 dirichlet in y, 3 neumann
    self.xval = xval
    self.yval = yval

class Timer:
  def __init__(self):
    self.start = time.time()

  def elapsed(self, message=""):
    self.end = time.time()
    duration = self.end - self.start
    print(f"Timer for {message}: {duration:.1f} seconds")

# =============================== GLOBAL VARIABLES ==============================
# ===============================================================================
global Uelas, Upf, D, pseudotime, basefilename, vtkextension, intrule, isAlignedMesh
Uelas = np.zeros(1) 
Upf = np.zeros(1) 
D = np.zeros((3, 3)) 
pseudotime = 0.0 
vtkextension = ".vtk" 
intrule = create2x2QuadratureRule() 



isAlignedMesh = False 
 
def get_default_nthreads(reserve=2):
  """
  Determine the default number of threads/processes to use based on the
  own system, instead of a fixed value in the code.
 
  Uses os.cpu_count(), which is portable (Linux/Windows/Mac) and reads the number
  of logical CPUs available directly from the operating system --
  equivalent to what 'lscpu' would show in "CPU(s)", but without depending
  on a specific command or text parsing.
 
  reserve: how many cores to leave free for the system/other tasks
           (default = 2). Garante pelo menos 2 thread.
  """
  ncpus = os.cpu_count()
  if ncpus is None:
    
    return 4
  return max(1, ncpus - reserve)
 

DEFAULT_NTHREADS = get_default_nthreads()

# Worker-side shared state for process pool execution
WORKER_NODES = None
WORKER_MAT = None

# =============================== FUNCTION IMPLEMENTATIONS ======================
# ===============================================================================


def readGmshMesh(filename):
    
    mesh = meshio.read(filename)
    
    nodes = []
    for pt in mesh.points:
        nodes.append(Node(pt[0], pt[1]))
        
    elements = []
    if "quad" in mesh.cells_dict:
        quad_cells = mesh.cells_dict["quad"]
        for cell_nodes in quad_cells:
            elements.append(Element(cell_nodes.tolist()))
    else:
        raise ValueError("The mesh does not contain quadrilaterals.!")
        
    # Process Physical Groups (boundary conditions tags)
    physical_groups = {}
    
    
    if "gmsh:physical" in mesh.cell_data_dict:
        
        for cell_type in mesh.cells_dict:
            if cell_type in mesh.cell_data_dict["gmsh:physical"]:
                cells = mesh.cells_dict[cell_type]
                tags = mesh.cell_data_dict["gmsh:physical"][cell_type]
                
                
                for i, tag in enumerate(tags):
                    if tag not in physical_groups:
                        physical_groups[tag] = set()
                   
                    physical_groups[tag].update(np.atleast_1d(cells[i]))
                    
    
    for tag in physical_groups:
        physical_groups[tag] = list(physical_groups[tag])

    # Build a name -> tag map from the physical names defined in the .geo file
    # (e.g. Physical Curve("Left_ids") = {4};), so BCs in the config can use
    # the string name instead of having to know the numeric Gmsh tag.
    physical_name_to_tag = {}
    for name, data in mesh.field_data.items():
      physical_name_to_tag[name] = data[0]

    print(f"Mesh loaded: {len(nodes)} nodes and {len(elements)} elements.")
       
    return nodes, elements, physical_groups, physical_name_to_tag


def createSparseStructure(K, elements, nstate):
  for element in elements:
    for i in range(4):
      row = nstate * element.node_ids[i]
      for j in range(4):
        col = nstate * element.node_ids[j]
        for k in range(nstate):
          for l in range(nstate):            
            K[row + k, col + l] = 1. 

def splitElementsAmongThreads(elements, nthreads):
  if nthreads <= 1:
    return [elements]
  chunk_size = max(1, (len(elements) + nthreads - 1) // nthreads)
  return [elements[i:i + chunk_size] for i in range(0, len(elements), chunk_size)]


def initProcessWorker(nodes, mat, upf, uelas, d, aligned_mesh):
  global WORKER_NODES, WORKER_MAT
  global Upf, Uelas, D, isAlignedMesh

  WORKER_NODES = nodes
  WORKER_MAT = mat
  Upf = upf
  Uelas = uelas
  D = d
  isAlignedMesh = aligned_mesh


def computeElementContribution(element, nodes, mat, nstate):
  nnodesel = len(element.node_ids)
  ndofel = nstate * nnodesel
  Ke = np.zeros((ndofel, ndofel))
  Fe = np.zeros(ndofel)
  computeElementStiffness(Ke, Fe, nodes, element, mat, nstate)
  return element.node_ids, Ke, Fe


def computeElementChunk(chunk_nstate):
  chunk, nstate = chunk_nstate
  local_updates = []
  for element in chunk:
    local_updates.append(computeElementContribution(element, WORKER_NODES, WORKER_MAT, nstate))
  return local_updates

def assembleGlobalStiffness(K, F, elements, nodes, mat, nstate, nthreads=None):
  
  if nthreads is None:
    nthreads = DEFAULT_NTHREADS
  for row in K.data: 
    for i in range(len(row)):
      row[i] = 0.0
  F.fill(0)
  if nthreads < 2:
    element_chunks = [elements]
  else:
    chunks = splitElementsAmongThreads(elements, nthreads)
    element_chunks = chunks

  if nthreads < 2:
    local_updates = []
    for element in elements:
      local_updates.append(computeElementContribution(element, nodes, mat, nstate))
    chunk_results = [local_updates]
  else:
    with ProcessPoolExecutor(
    max_workers=nthreads,
    initializer=initProcessWorker,
    initargs=(nodes, mat, Upf, Uelas, D, isAlignedMesh)
) as executor:
      chunk_results = list(executor.map(computeElementChunk, [(chunk, nstate) for chunk in element_chunks]))

  for local_updates in chunk_results:
    for node_ids, Ke, Fe in local_updates:
      nnodesel = len(node_ids)
      for i in range(nnodesel):
        row = nstate * node_ids[i]
        for k in range(nstate):
          F[row + k] += Fe[nstate * i + k]
        for j in range(nnodesel):
          col = nstate * node_ids[j]
          for k in range(nstate):
            for l in range(nstate):
              index = K.rows[row+k].index(col+l)
              K.data[row+k][index] += Ke[nstate * i + k, nstate * j + l]
  
  


def computeReaction(K, F, nodes, elements, mat, target_ids, dof_offset=0, sign_factor=-1.0, nthreads=None):
    assembleGlobalStiffness(K, F, elements, nodes, mat, 2, nthreads=nthreads)
    residual = K.dot(Uelas)  # F is zero
    reaction = 0.0
    for i in target_ids:
        # 2*i is the X direction, 2*i + 1 is the Y direction.
        reaction += sign_factor * residual[2 * i + dof_offset]
    return reaction
  
def jacobian(nodes, dN):
  coords = np.array([[n.x, n.y] for n in nodes]).T  # shape (2, 4)
  J = coords @ dN.T  # shape (2, 2)
  detjac = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
  if abs(detjac) < 1e-8:
    raise ValueError("Jacobian determinant is too small, check the mesh or element shape")
  J_inv = np.array([[ J[1, 1], -J[0, 1]],
                    [-J[1, 0],  J[0, 0]]]) / detjac
  return J_inv, detjac

def constantJacobian(elnodes):
  n1, n2, n3, n4 = elnodes
  base = n2.x - n1.x
  height = n4.y - n1.y
  area = base * height
  detjac = area / 4.0
  dqsidx = 2.0 / base
  detady = 2.0 / height
  J_inv = np.diag([dqsidx, detady])
  return J_inv, detjac

def computeElementStiffness(Ke, Fe, nodes, element, mat, nstate):
  nnodes = len(element.node_ids)
  elnodes = [nodes[i] for i in element.node_ids] 
  if isAlignedMesh:
    J_inv, detjac = constantJacobian(elnodes)

  if nstate == 2: # compute elasticity stiffness
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      if not isAlignedMesh:
        J_inv, detjac = jacobian(elnodes, dN)
      dN_xy = J_inv.T @ dN
      B = createB(dN_xy)
      phase_field = sum(N[0, nstate * i] * Upf[element.node_ids[i]] for i in range(nnodes))
      Ddeteriorated = D.copy()
      Ddeteriorated *= (1 - phase_field) ** 2
      Ke += B.T @ Ddeteriorated @ B * qp.weight * abs(detjac)
  elif nstate == 1: # compute phase field stiffness
    Gc, l0 = mat.Gc, mat.l0
    c0 = 2.0
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      if not isAlignedMesh:
        J_inv, detjac = jacobian(elnodes, dN)   
      dN_xy = J_inv.T @ dN # Same as B_phi
      sigmaDotEps = calculateTensileSigmaDotEps(element, dN_xy)
      Ke += abs(detjac) * qp.weight * (Gc * l0 / c0 * (dN_xy.T @ dN_xy) + (Gc / (l0 * c0) + 0.5 * sigmaDotEps) * N.T @ N)
      Fe += abs(detjac) * qp.weight * 0.5 * sigmaDotEps * N.flatten()
  else:
    raise Exception("Invalid nstate")

def computeSigmaAtCenter(element, nodes, stress_vec):
  qsi, eta = 0.0, 0.0
  n1, n2, n3, n4 = [nodes[i] for i in element.node_ids]
  base = n2.x - n1.x
  height = n4.y - n1.y
  dqsidx = 2.0 / base
  dqsidy = 2.0 / height
  J_inv = np.diag([dqsidx, dqsidy])
  N, dN = shapeFunctions(qsi, eta, 2)
  dN_xy = J_inv.T @ dN
  dU = np.zeros((2, 2))
  for i in range(4):
    index = 2 * element.node_ids[i]
    dU[0, 0] += dN_xy[0, i] * Uelas[index]
    dU[0, 1] += dN_xy[1, i] * Uelas[index]
    dU[1, 0] += dN_xy[0, i] * Uelas[index + 1]
    dU[1, 1] += dN_xy[1, i] * Uelas[index + 1]
  strain = 0.5 * (dU + dU.T)
  strain_vec = np.array([strain[0, 0], strain[1, 1], 2 * strain[0, 1]])
  phase_field = sum(N[0, 2 * i] * Upf[element.node_ids[i]] for i in range(4))
  g = (1.0 - phase_field) ** 2
  stress_vec[:] = g * D @ strain_vec
  return phase_field

def calculateTensileSigmaDotEps(element, dN):
  dU = np.zeros((2, 2))
  for i in range(4):
    index = 2 * element.node_ids[i]
    dU[0, 0] += dN[0, i] * Uelas[index]
    dU[0, 1] += dN[1, i] * Uelas[index]
    dU[1, 0] += dN[0, i] * Uelas[index + 1]
    dU[1, 1] += dN[1, i] * Uelas[index + 1]
  strain = 0.5 * (dU + dU.T)
  strain_vec = np.array([strain[0, 0], strain[1, 1], 2 * strain[0, 1]])
  sigma_vec = D @ strain_vec
  sigma = np.array([[sigma_vec[0], sigma_vec[2]],
                    [sigma_vec[2], sigma_vec[1]]])

  eigsig = np.linalg.eigvals(sigma)
  eigstrain = np.linalg.eigvals(strain)
  sigmaDotEps = sum(eigsig[i] * eigstrain[i] if eigstrain[i] > 0 else 0.0 for i in range(len(eigsig)))

  return sigmaDotEps

def shapeFunctions(qsi, eta, nstate):
  phi1qsi = (1 + qsi) / 2.0
  phi0eta = (1 - eta) / 2.0
  phi1eta = (1 + eta) / 2.0
  phi0qsi = (1 - qsi) / 2.0
  shape = np.array([phi0qsi * phi0eta, phi1qsi * phi0eta, phi1qsi * phi1eta, phi0qsi * phi1eta])
  N = np.zeros((nstate, nstate * 4))
  if nstate == 1:
    N[0, :4] = shape
  else:
    for i in range(4):
      N[0, 2 * i] = shape[i]
      N[1, 2 * i + 1] = shape[i]
  dN = np.array([
    [0.25 * (-1 + eta), 0.25 * (1 - eta), 0.25 * (1 + eta), 0.25 * (-1 - eta)],
    [0.25 * (-1 + qsi), 0.25 * (-1 - qsi), 0.25 * (1 + qsi), 0.25 * (1 - qsi)]
  ])
  return N, dN

def createB(dN):
  B = np.zeros((3, 8))
  for i in range(4):
    B[0, 2 * i] = dN[0, i]
    B[1, 2 * i + 1] = dN[1, i]
    B[2, 2 * i] = dN[1, i]
    B[2, 2 * i + 1] = dN[0, i]
  return B

def zeroRowAndColumnOfSparseMatrix(K, target_row):
  K.data[target_row] = [0] * len(K.data[target_row])
  for i in range(K.shape[0]):    
    if target_row in K.rows[i]:
        col_idx = K.rows[i].index(target_row)
        K.data[i][col_idx] = 0

def applyBoundaryConditions(K, F, bc_nodes):
  
  isPenalty = True
  if isPenalty:
    B = 1e10
    for bc in bc_nodes:
      row = 2 * bc.node
      xval = bc.xval * pseudotime
      yval = bc.yval * pseudotime
      if bc.bctype == 0: 
        K[row, row] += B
        K[row + 1, row + 1] += B
        F[row] += B * xval
        F[row + 1] += B * yval
      elif bc.bctype == 1: 
        K[row, row] += B  
        F[row] += B * xval
      elif bc.bctype == 2: 
        K[row + 1, row + 1] += B
        F[row + 1] += B * yval
      elif bc.bctype == 3: 
        F[row] += xval
        F[row + 1] += yval
  else:
    for bc in bc_nodes:
      row = 2 * bc.node
      xval = bc.xval * pseudotime
      yval = bc.yval * pseudotime
      if bc.bctype == 0: 
        F -= K[row,:].toarray().flatten() * xval
        F -= K[row + 1,:].toarray().flatten() * yval
        zeroRowAndColumnOfSparseMatrix(K, row)
        zeroRowAndColumnOfSparseMatrix(K, row+1)
        K[row, row] = 1.0
        K[row + 1, row + 1] = 1.0
        F[row] = xval
        F[row + 1] = yval
      elif bc.bctype == 1: 
        F -= K[row,:].toarray().flatten() * xval
        zeroRowAndColumnOfSparseMatrix(K, row)
        K[row, row] = 1.0
        F[row] = xval
      elif bc.bctype == 2: 
        F -= K[row + 1,:].toarray().flatten() * yval
        zeroRowAndColumnOfSparseMatrix(K, row+1)
        K[row + 1, row + 1] = 1.0
        F[row + 1] = yval
      elif bc.bctype == 3: 
        F[row] += xval
        F[row + 1] += yval
  

def solveSystem(K, F, U):
  
  U[:] = sparse.linalg.splu(K.tocsc()).solve(F)
  

def generateVTKLegacyFile(nodes, elements, filename):
  with open(filename, 'w') as vtkFile:
    vtkFile.write("# vtk DataFile Version 2.0\n")
    vtkFile.write("FEM results\n")
    vtkFile.write("ASCII\n")
    vtkFile.write("DATASET UNSTRUCTURED_GRID\n")
    vtkFile.write(f"POINTS {len(nodes)} float\n")
    for node in nodes:
      vtkFile.write(f"{node.x} {node.y} 0.0\n")
    vtkFile.write(f"CELLS {len(elements)} {len(elements) * 5}\n")
    for element in elements:
      vtkFile.write(f"4 {' '.join(map(str, element.node_ids))}\n")
    vtkFile.write(f"CELL_TYPES {len(elements)}\n")
    for _ in elements:
      vtkFile.write("9\n")
    vtkFile.write(f"POINT_DATA {len(nodes)}\n")
    vtkFile.write("VECTORS displacements float\n")
    for i in range(len(nodes)):
      vtkFile.write(f"{Uelas[2 * i]} {Uelas[2 * i + 1]} 0.0\n")
    vtkFile.write("SCALARS phasefield float 1\n")
    vtkFile.write("LOOKUP_TABLE default\n")
    for i in range(len(nodes)):
      vtkFile.write(f"{Upf[i]}\n")


def create_bc_nodes_from_config(nodes, boundary_conditions, physical_groups=None, physical_name_to_tag=None):
  """
  Creates BC objects from the configuration.
  
  If node_filter is a function: uses the function to filter nodes
  If node_filter is an int: treated as the numeric Gmsh Physical Group tag
  If node_filter is a str: resolved via physical_name_to_tag to the matching
                            Physical Group tag (e.g. "Left_ids" -> 1), using
                            the names defined in the .geo file
                            (Physical Curve("Left_ids") = {...};)
  """
  bc_nodes = []
  
  for bc_config in boundary_conditions:

    if physical_groups and isinstance(bc_config.node_filter, (str, int)):

      node_filter = bc_config.node_filter

      # Resolve string names (e.g. "Left_ids") to the numeric Gmsh tag
      if isinstance(node_filter, str):
        if not physical_name_to_tag or node_filter not in physical_name_to_tag:
          available = list(physical_name_to_tag.keys()) if physical_name_to_tag else []
          raise ValueError(
              f"Boundary condition '{bc_config.name}' references physical group "
              f"name '{node_filter}', but it was not found in the mesh.\n"
              f"Available physical group names: {available}"
          )
        node_filter = physical_name_to_tag[node_filter]

      node_ids = physical_groups.get(node_filter, [])
      for node_id in node_ids:
        bc_nodes.append(BC(node_id, bc_config.bc_type, bc_config.xval, bc_config.yval))
  
      
  return bc_nodes


def resolve_reaction_target_ids(physical_groups, physical_name_to_tag, reaction_config):
  """Resolve target node IDs for reaction-force monitoring.

  Preference order:
  1. Exact physical-group name from the mesh (e.g. 'Bottom_ids').
  2. Legacy aliases used by older configs (e.g. 'bottom_ids_y', 'supports_y').
  """
  if not reaction_config:
    return []

  reaction_type = (reaction_config.reaction_type or "").strip()
  reaction_type_lower = reaction_type.lower()
  name_lookup = {name.lower(): name for name in (physical_name_to_tag or {})}

  def ids_for_name(group_name):
    actual_name = name_lookup.get(group_name.lower())
    if actual_name is None:
      return []
    tag = physical_name_to_tag[actual_name]
    return physical_groups.get(tag, [])

  base_name = reaction_type.rsplit("_", 1)[0] if reaction_type_lower.endswith(("_x", "_y")) else reaction_type
  target_ids = ids_for_name(base_name)
  if target_ids:
    return target_ids

  if "bottom" in reaction_type_lower:
    return ids_for_name("Bottom_ids")
  if "top" in reaction_type_lower:
    return ids_for_name("Top_ids")
  if "left" in reaction_type_lower:
    return ids_for_name("Left_ids")
  if "right" in reaction_type_lower:
    return ids_for_name("Right_ids")
  if "support" in reaction_type_lower:
    target_ids = ids_for_name("Left_ids") + ids_for_name("Right_ids")
    if target_ids:
      return target_ids

  return []

# =============================== MAIN ==========================================
# ===============================================================================


def main(config_name='default'):
  if not os.path.exists("outputs"):
    os.makedirs("outputs")
  try:
    from config_simulations import get_simulation_config
    print(f"\n Loading configuration: {config_name}")
    config = get_simulation_config(config_name)
    print(f" {config.name}\n")
    
    # Material
    E = config.material.E
    nu = config.material.nu
    Gc = config.material.Gc
    l0 = config.material.l0
    length = config.material.length
    
    # Simulation parameters
    dt = config.simulation_params['dt']
    totaltime = config.simulation_params['totaltime']
    maxsteps = int(config.simulation_params['maxsteps'])
    maxiter = config.simulation_params['maxiter']
    stagtol = config.simulation_params['stagtol']
    
    # Output file
    global basefilename
    basefilename = config.output_base
    
    if config.mesh_type == 'gmsh':
      mesh_file = os.path.join("simulations", config.mesh_file)
      print(f" Used mesh: {mesh_file}")
      nodes, elements, physical_groups, physical_name_to_tag = readGmshMesh(mesh_file)
    else:
      raise ValueError(f"Mesh type '{config.mesh_type}' not supported. Please use 'gmsh'.")
    
  except Exception as e:
    print(f"\n Error loading configuration: {e}")
    print("Could not read simulation parameters.")
    print(f"Requested configuration: {config_name}")
    print("Please check the 'config_simulations.py' file and try again.\n")
    return

  simulation_time = Timer()

  # Plot data
  imposed_displacement = config.imposed_displacement if config else 0.04
  if 'simplebar.msh' in mesh_file:
    sigma_peak_at2 = np.sqrt(27.0 * E * Gc / (256.0 * l0))
    u_peak_at2 = 16.0 / 9.0 * sigma_peak_at2* length / E
    print(f"Sigma peak: {sigma_peak_at2}")
    print(f"U peak: {u_peak_at2}\n")

  print(f"Imposed displacement: {imposed_displacement}")

  material = MaterialParameters(E, nu, Gc, l0, material_type=config.material.material_type)
  global D
  D = material.get_D_matrix()

  # Apply boundary conditions from the configuration
  print(f" Applying {len(config.boundary_conditions)} boundary conditions from the configuration...")
  bc_nodes = create_bc_nodes_from_config(nodes, config.boundary_conditions, physical_groups, physical_name_to_tag)

  nstate_elas = 2
  nstate_pf = 1
  nnodes = len(nodes)
  ndofs_elas = nstate_elas * nnodes
  ndofs_pf = nstate_pf * nnodes
  Kelas = sparse.lil_matrix((ndofs_elas, ndofs_elas))
  createSparseStructure(Kelas, elements, nstate_elas)
  Felas = np.zeros(ndofs_elas)
  global Uelas
  Uelas = np.zeros(ndofs_elas)

  Kpf = sparse.lil_matrix((ndofs_pf, ndofs_pf))
  createSparseStructure(Kpf, elements, nstate_pf)
  Fpf = np.zeros(ndofs_pf)  
  global Upf
  Upf = np.zeros(ndofs_pf)

  u_data = []
  force_data = []
  stress_data = []
  time_data = []

  global pseudotime
  mythread = config.simulation_params.get('nthreads', DEFAULT_NTHREADS)
  print(f" Using {mythread} threads (CPUs detected: {os.cpu_count()})")
  pseudotime = 0.0
  for step in range(maxsteps):
    pseudotime += dt
    if pseudotime > totaltime:
      break
    print(f"******************** Time Step {step} | Pseudo time = {pseudotime:.6f} | Time step = {dt} ********************")
    for iter in range(maxiter):      
      print(f"------ Staggered Iteration {iter} ------")
      assembleGlobalStiffness(Kelas, Felas, elements, nodes, material, nstate_elas, nthreads=mythread)
      applyBoundaryConditions(Kelas, Felas, bc_nodes)
      if iter != 0:
        residual = Kelas.tocsr() @ Uelas - Felas
        norm = np.linalg.norm(residual)
        print(f"Residual Elasticity Norm: {norm:.2e}")
        if norm < stagtol:
          print(f"------> Staggered scheme converged in {iter} iterations.")
          break
      solveSystem(Kelas, Felas, Uelas)
      assembleGlobalStiffness(Kpf, Fpf, elements, nodes, material, nstate_pf, nthreads=mythread)
      solveSystem(Kpf, Fpf, Upf)
    if iter == maxiter:
      print(f"------> Staggered scheme did not converge in {maxiter} iterations.\nAccepting current solution and continuing")
    filename = f"{basefilename}{step}{vtkextension}"
    generateVTKLegacyFile(nodes, elements, filename)

    # Collect data for graphs based on the configuration
    if config.graph_config:
      if config.graph_config.graph_type == 'stress_vs_time':
        
        element_id = config.graph_config.element_id or 450
        if element_id < len(elements):
          sig = np.zeros(3)
          pfmid = computeSigmaAtCenter(elements[element_id], nodes, sig)
          stress_data.append(sig[0] / sigma_peak_at2)
          time_data.append(pseudotime)
      elif config.graph_config.graph_type == 'force_vs_displacement':

        target_ids = resolve_reaction_target_ids(physical_groups, physical_name_to_tag, config.reaction_config)
        dof_offset = 0
        sign_factor = -1.0
        
        if config.reaction_config:
            dof_offset = config.reaction_config.reaction_dof
            sign_factor = config.reaction_config.sign_factor
            
        # Collect reaction force data
        reaction = computeReaction(Kelas, Felas, nodes, elements, material, target_ids, dof_offset, sign_factor)
        force_data.append(reaction)
        
        # Determine which displacement to use based on the configuration
        displacement_value = 0.0
        if config.graph_config.displacement_axis == 'x':
         
          for bc in config.boundary_conditions:
            if bc.xval != 0:
              displacement_value = bc.xval
              break
        else:  # 'y'

          for bc in config.boundary_conditions:
            if bc.yval != 0:
              displacement_value = abs(bc.yval)
              break
        
        u_data.append(pseudotime * displacement_value)
  
  # Generate graph if requested in the configuration
  if config.graph_config:
    plt.figure()
    if config.graph_config.graph_type == 'stress_vs_time':
      plt.plot(time_data, stress_data, 'o')
    elif config.graph_config.graph_type == 'force_vs_displacement':
      plt.plot(u_data, force_data, 'o')
    
    plt.xlabel(config.graph_config.x_label)
    plt.ylabel(config.graph_config.y_label)
    plt.title(config.graph_config.title)
    plt.grid(True)
    plt.savefig(config.graph_config.output_file)
    print(f"✓ Graph saved to: {config.graph_config.output_file}")

  print("\n================> Simulation completed!")
  simulation_time.elapsed("complete simulation")

if __name__ == "__main__":
  config_name = sys.argv[1] if len(sys.argv) > 1 else 'default'
  
  print("\n" + "="*80)
  print(f"phasefieldjr- Simulation: {config_name}")
  print("="*80 + "\n")
  
  main(config_name)