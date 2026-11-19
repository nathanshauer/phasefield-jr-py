# MIT License
# © 2025 Nathan Shauer
# phasefield-jr

"""Two-dimensional orthotropic phase-field finite-element solver.

This version preserves the original numerical formulation while reducing
repeated work during element and sparse-matrix assembly. The model assumes
unit out-of-plane thickness.
"""

from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import os
import sys
import time

import numpy as np


# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# DATA STRUCTURES
# =============================================================================

class QuadraturePoint:
  __slots__ = ("xi", "eta", "weight")

  def __init__(self, xi, eta, weight):
    self.xi = xi
    self.eta = eta
    self.weight = weight


class MaterialParameters:
  """Orthotropic elastic and phase-field material parameters."""

  __slots__ = (
    "E1",
    "E2",
    "nu12",
    "G12",
    "Gc1",
    "Gc2",
    "l0",
    "material_type",
    "theta_deg",
    "_Gc_reference",
    "_fracture_tensor",
  )

  def __init__(
      self,
      E1,
      E2,
      nu12,
      G12,
      Gc1,
      Gc2,
      l0,
      material_type="plane_stress",
      theta_deg=0.0,
  ):
    self.E1 = E1
    self.E2 = E2
    self.nu12 = nu12
    self.G12 = G12
    self.Gc1 = Gc1
    self.Gc2 = Gc2
    self.l0 = l0
    self.material_type = material_type
    self.theta_deg = theta_deg
    self._validate()
    self._Gc_reference = np.sqrt(self.Gc1 * self.Gc2)
    self._fracture_tensor = self._build_fracture_tensor()

  def _validate(self):
    for name in ("E1", "E2", "G12", "l0"):
      value = getattr(self, name)
      if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number.")
    for name in ("Gc1", "Gc2"):
      value = getattr(self, name)
      if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number.")
    if not np.isfinite(self.nu12):
      raise ValueError("nu12 must be finite.")
    if not np.isfinite(self.theta_deg):
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

  def is_isotropic_equivalent(self):
    return (
      np.isclose(self.E1, self.E2)
      and np.isclose(
        self.G12,
        self.E1 / (2.0 * (1.0 + self.nu12)),
      )
    )

  def is_fracture_isotropic(self):
    return np.isclose(self.Gc1, self.Gc2)

  def get_reference_toughness(self):
    """Return the symmetric normalization sqrt(Gc1*Gc2)."""
    return self._Gc_reference

  def _build_fracture_tensor(self):
    """Build the crack-gradient tensor from propagation toughnesses.

    Gc1 and Gc2 refer to crack advance along material axes 1 and 2.
    The gradient is normal to the crack, so Gc2 weights normal axis 1
    and Gc1 weights normal axis 2.
    """
    theta = np.radians(self.theta_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    rotation = np.array([[c, -s], [s, c]])

    principal = np.diag([
      self.Gc2 ** 2 / self._Gc_reference,
      self.Gc1 ** 2 / self._Gc_reference,
    ])
    return rotation @ principal @ rotation.T

  def get_fracture_tensor(self):
    """Return a copy of the crack-gradient tensor in global axes."""
    return self._fracture_tensor.copy()

  def get_fracture_toughness(self, normal):
    """Evaluate effective toughness for a crack unit-normal direction."""
    normal = np.asarray(normal, dtype=float)
    if normal.shape != (2,) or not np.all(np.isfinite(normal)):
      raise ValueError("normal must contain two finite components.")
    magnitude = np.linalg.norm(normal)
    if magnitude == 0.0:
      raise ValueError("normal must be non-zero.")
    unit_normal = normal / magnitude
    value = (
      self._Gc_reference
      * unit_normal @ self._fracture_tensor @ unit_normal
    )
    return np.sqrt(max(0.0, value))

  @staticmethod
  def _get_inplane_Q(E1, E2, nu12, G12):
    """Return Q11, Q12, Q22 and Q66 in the material 1-2 axes."""
    nu21 = nu12 * E2 / E1
    denominator = 1.0 - nu12 * nu21
    return (
      E1 / denominator,
      nu12 * E2 / denominator,
      E2 / denominator,
      G12,
    )

  def get_D_matrix(self):
    """Build the in-plane constitutive matrix in the global x-y axes.

    Genuine orthotropic plane strain requires out-of-plane constants that are
    not present in this model. Consequently, plane strain remains available
    only for an isotropic-equivalent material, as in the original solver.
    """
    if self.material_type == "plane_stress":
      E1, E2 = self.E1, self.E2
      nu12, G12 = self.nu12, self.G12
    elif self.material_type == "plane_strain":
      if not self.is_isotropic_equivalent():
        raise NotImplementedError(
          "material_type='plane_strain' is only supported for materials "
          "that behave isotropically in-plane. Genuine orthotropic plane "
          "strain requires E3, nu13 and nu23."
        )

      E_star = self.E1 / (1.0 - self.nu12 ** 2)
      nu_star = self.nu12 / (1.0 - self.nu12)
      E1, E2, nu12, G12 = E_star, E_star, nu_star, self.G12
    else:
      raise ValueError(
        f"Unknown material_type '{self.material_type}'. "
        "Use 'plane_stress' or 'plane_strain'."
      )

    Q11, Q12, Q22, Q66 = self._get_inplane_Q(E1, E2, nu12, G12)

    theta = np.radians(self.theta_deg)
    c = np.cos(theta)
    s = np.sin(theta)

    # Keep the original operation order: this matrix is built only once, so
    # changing powers into reused products offers no meaningful speed benefit
    # and can introduce tiny round-off differences for off-axis materials.
    Qbar11 = Q11 * c**4 + 2 * (Q12 + 2 * Q66) * s**2 * c**2 + Q22 * s**4
    Qbar12 = (
      (Q11 + Q22 - 4 * Q66) * s**2 * c**2
      + Q12 * (s**4 + c**4)
    )
    Qbar22 = Q11 * s**4 + 2 * (Q12 + 2 * Q66) * s**2 * c**2 + Q22 * c**4
    Qbar16 = (
      (Q11 - Q12 - 2 * Q66) * s * c**3
      + (Q12 - Q22 + 2 * Q66) * s**3 * c
    )
    Qbar26 = (
      (Q11 - Q12 - 2 * Q66) * s**3 * c
      + (Q12 - Q22 + 2 * Q66) * s * c**3
    )
    Qbar66 = (
      (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s**2 * c**2
      + Q66 * (s**4 + c**4)
    )

    return np.array([
      [Qbar11, Qbar12, Qbar16],
      [Qbar12, Qbar22, Qbar26],
      [Qbar16, Qbar26, Qbar66],
    ])


class Element:
  __slots__ = (
    "node_ids",
    "_integration_cache_key",
    "_integration_data",
    "_assembly_data",
  )

  def __init__(self, node_ids):
    self.node_ids = node_ids
    self._integration_cache_key = None
    self._integration_data = None
    self._assembly_data = {}


class Node:
  __slots__ = ("x", "y")

  def __init__(self, x, y):
    self.x = x
    self.y = y


class BC:
  __slots__ = ("node", "bctype", "xval", "yval")

  def __init__(self, node, bc_type, xval, yval):
    self.node = node
    self.bctype = bc_type
    self.xval = xval
    self.yval = yval


class Timer:
  __slots__ = ("start", "end")

  def __init__(self):
    self.start = time.perf_counter()

  def elapsed(self, message=""):
    self.end = time.perf_counter()
    duration = self.end - self.start
    print(f"Timer for {message}: {duration:.1f} seconds")


# =============================================================================
# NUMERICAL STATE AND WORKER STATE
# =============================================================================

elastic_state_size = 2
phase_field_state_size = 1

Uelas = np.zeros(1)
Upf = np.zeros(1)
D = np.zeros((3, 3))
pseudotime = 0.0
isAlignedMesh = False

WORKER_NODES = None
WORKER_MATERIAL = None
WORKER_UELAS = None
WORKER_UPF = None
WORKER_D = None


def create2x2QuadratureRule():
  point = 1.0 / np.sqrt(3.0)
  return tuple(
    QuadraturePoint(xi, eta, 1.0)
    for xi in (-point, point)
    for eta in (-point, point)
  )


INTEGRATION_RULE = create2x2QuadratureRule()


def get_default_nthreads(reserve=2):
  """Return available logical CPUs while leaving ``reserve`` CPUs free."""
  cpu_count = os.cpu_count()
  if cpu_count is None:
    return 4
  return max(1, cpu_count - reserve)


DEFAULT_NTHREADS = get_default_nthreads()


# =============================================================================
# MESH AND ELEMENT GEOMETRY
# =============================================================================

def readGmshMesh(filename):
  import meshio

  mesh = meshio.read(filename)

  nodes = [Node(point[0], point[1]) for point in mesh.points]

  if "quad" not in mesh.cells_dict:
    raise ValueError("The mesh does not contain quadrilateral elements.")
  elements = [Element(cell.tolist()) for cell in mesh.cells_dict["quad"]]

  physical_groups = {}
  physical_data = mesh.cell_data_dict.get("gmsh:physical", {})
  for cell_type, cells in mesh.cells_dict.items():
    if cell_type not in physical_data:
      continue
    for cell, tag in zip(cells, physical_data[cell_type]):
      physical_groups.setdefault(tag, set()).update(np.atleast_1d(cell))

  for tag, node_ids in physical_groups.items():
    physical_groups[tag] = list(node_ids)

  physical_name_to_tag = {
    name: data[0]
    for name, data in mesh.field_data.items()
  }

  print(f"Mesh loaded: {len(nodes)} nodes and {len(elements)} elements.")
  return nodes, elements, physical_groups, physical_name_to_tag


@lru_cache(maxsize=8)
def shapeFunctions(qsi, eta, nstate):
  """Return Q4 interpolation and parent-coordinate derivatives.

  Only eight combinations are used (four quadrature points and two state
  sizes), so caching removes thousands of identical small allocations.
  Returned arrays are treated as read-only by the solver.
  """
  phi1qsi = (1 + qsi) / 2.0
  phi0eta = (1 - eta) / 2.0
  phi1eta = (1 + eta) / 2.0
  phi0qsi = (1 - qsi) / 2.0
  shape = np.array([
    phi0qsi * phi0eta,
    phi1qsi * phi0eta,
    phi1qsi * phi1eta,
    phi0qsi * phi1eta,
  ])

  N = np.zeros((nstate, nstate * 4))
  if nstate == phase_field_state_size:
    N[0, :4] = shape
  else:
    for i in range(4):
      N[0, 2 * i] = shape[i]
      N[1, 2 * i + 1] = shape[i]

  dN = np.array([
    [
      0.25 * (-1 + eta),
      0.25 * (1 - eta),
      0.25 * (1 + eta),
      0.25 * (-1 - eta),
    ],
    [
      0.25 * (-1 + qsi),
      0.25 * (-1 - qsi),
      0.25 * (1 + qsi),
      0.25 * (1 - qsi),
    ],
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


def jacobian(nodes, dN):
  coordinates = np.array([[node.x, node.y] for node in nodes]).T
  J = coordinates @ dN.T
  determinant = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
  if abs(determinant) < 1.0e-8:
    raise ValueError(
      "Jacobian determinant is too small; check the mesh element."
    )
  inverse = np.array([
    [J[1, 1], -J[0, 1]],
    [-J[1, 0], J[0, 0]],
  ]) / determinant
  return inverse, determinant


def constantJacobian(element_nodes):
  n1, n2, _, n4 = element_nodes
  base = n2.x - n1.x
  height = n4.y - n1.y
  determinant = base * height / 4.0
  return np.diag([2.0 / base, 2.0 / height]), determinant


def _get_element_integration_data(element, nodes):
  """Return cached geometry and interpolation data for one element."""
  cache_key = (id(nodes), isAlignedMesh)
  if element._integration_cache_key == cache_key:
    return element._integration_data

  element_nodes = [nodes[node_id] for node_id in element.node_ids]
  if isAlignedMesh:
    constant_inverse, constant_determinant = constantJacobian(element_nodes)

  integration_data = []
  for quadrature_point in INTEGRATION_RULE:
    N_elastic, dN = shapeFunctions(
      quadrature_point.xi,
      quadrature_point.eta,
      elastic_state_size,
    )
    N_phase, _ = shapeFunctions(
      quadrature_point.xi,
      quadrature_point.eta,
      phase_field_state_size,
    )

    if isAlignedMesh:
      inverse, determinant = constant_inverse, constant_determinant
    else:
      inverse, determinant = jacobian(element_nodes, dN)

    dN_xy = inverse.T @ dN
    integration_data.append((
      N_elastic,
      N_phase,
      dN_xy,
      createB(dN_xy),
      quadrature_point.weight * abs(determinant),
    ))

  element._integration_cache_key = cache_key
  element._integration_data = tuple(integration_data)
  return element._integration_data


def prepareElementIntegrationData(nodes, elements):
  """Precompute geometry used in every staggered iteration."""
  for element in elements:
    _get_element_integration_data(element, nodes)


# =============================================================================
# SPARSE STRUCTURE AND ASSEMBLY
# =============================================================================

# Starting processes costs more than serial assembly for small meshes. This
# value can be overridden with simulation_params['parallel_min_elements'].
# Elastic assembly is much cheaper per element than phase-field assembly, so
# it uses a higher threshold before process startup becomes worthwhile.
default_parallel_min_elements = 2_000
elastic_parallel_threshold_multiplier = 20

def _element_dofs(node_ids, nstate):
  return tuple(
    nstate * node_id + local_dof
    for node_id in node_ids
    for local_dof in range(nstate)
  )


def createSparseStructure(K, elements, nodes_or_nstate, nstate=None):
  """Create the LIL pattern and cache local insertion indices.

  Both ``(K, elements, nstate)`` from the attached solver and the previous
  optimized ``(K, elements, nodes, nstate)`` call are accepted. Node data are
  not needed to construct the connectivity pattern.
  """
  if nstate is None:
    nstate = nodes_or_nstate
  row_columns = [set() for _ in range(K.shape[0])]

  for element in elements:
    dofs = _element_dofs(element.node_ids, nstate)
    for row in dofs:
      row_columns[row].update(dofs)

  for row, columns in enumerate(row_columns):
    sorted_columns = sorted(columns)
    K.rows[row] = sorted_columns
    K.data[row] = [1.0] * len(sorted_columns)

  row_positions = [
    {column: position for position, column in enumerate(columns)}
    for columns in K.rows
  ]

  for element in elements:
    dofs = _element_dofs(element.node_ids, nstate)
    positions = tuple(
      tuple(row_positions[row][column] for column in dofs)
      for row in dofs
    )
    element._assembly_data[nstate] = (dofs, positions)


def splitElementsAmongThreads(elements, nthreads):
  if nthreads <= 1:
    return [elements]
  chunk_size = max(1, (len(elements) + nthreads - 1) // nthreads)
  return [
    elements[start:start + chunk_size]
    for start in range(0, len(elements), chunk_size)
  ]


def initProcessWorker(nodes, mat, Uelas_data, Upf_data, D_mat):
  global WORKER_NODES, WORKER_MATERIAL, WORKER_UELAS, WORKER_UPF, WORKER_D
  WORKER_NODES = nodes
  WORKER_MATERIAL = mat
  WORKER_UELAS = Uelas_data
  WORKER_UPF = Upf_data
  WORKER_D = D_mat


def computeElementContribution(
    element,
    nodes,
    mat,
    nstate,
    Uelas_data=None,
    Upf_data=None,
    D_mat=None,
    tensile_values=None,
):
  """Compute one element contribution without touching the global matrix."""
  if Uelas_data is None:
    Uelas_data = Uelas
  if Upf_data is None:
    Upf_data = Upf
  if D_mat is None:
    D_mat = D

  local_size = nstate * len(element.node_ids)
  Ke = np.zeros((local_size, local_size))
  Fe = np.zeros(local_size)
  computeElementStiffness(
    Ke,
    Fe,
    nodes,
    element,
    mat,
    nstate,
    Uelas_data,
    Upf_data,
    D_mat,
    tensile_values,
  )
  return element.node_ids, Ke, Fe


def computeElementChunk(chunk_nstate):
  chunk, nstate = chunk_nstate
  return [
    computeElementContribution(
      element,
      WORKER_NODES,
      WORKER_MATERIAL,
      nstate,
      WORKER_UELAS,
      WORKER_UPF,
      WORKER_D,
    )
    for element in chunk
  ]


def _zero_lil_matrix(K):
  for row_data in K.data:
    row_data[:] = [0.0] * len(row_data)


def _assemble_local_contribution(K, F, element, Ke, Fe, nstate):
  dofs, positions = element._assembly_data[nstate]
  add_element_force = nstate == phase_field_state_size

  # The loop order matches the original assembly exactly: element, local row,
  # local column. Cached positions remove the expensive list.index operation.
  for local_row, global_row in enumerate(dofs):
    if add_element_force:
      F[global_row] += Fe[local_row]
    matrix_row = K.data[global_row]
    local_values = Ke[local_row]
    for local_column, position in enumerate(positions[local_row]):
      matrix_row[position] += local_values[local_column]


def _effective_worker_count(nthreads, number_of_elements,
                            parallel_min_elements):
  requested = max(1, int(nthreads))
  if requested < 2 or number_of_elements < parallel_min_elements:
    return 1
  return min(requested, number_of_elements)


def assembleGlobalStiffness(
    K,
    F,
    elements,
    nodes,
    mat,
    nstate,
    nthreads=None,
    parallel_min_elements=default_parallel_min_elements,
):
  """Assemble global stiffness and force arrays in deterministic order."""
  if nthreads is None:
    nthreads = DEFAULT_NTHREADS

  _zero_lil_matrix(K)
  F.fill(0.0)

  worker_count = _effective_worker_count(
    nthreads,
    len(elements),
    (
      parallel_min_elements * elastic_parallel_threshold_multiplier
      if nstate == elastic_state_size
      else parallel_min_elements
    ),
  )

  if worker_count == 1:
    tensile_values = None
    if nstate == phase_field_state_size:
      tensile_values = _calculate_all_tensile_values(
        elements,
        nodes,
        Uelas,
        D,
      )

    for element_index, element in enumerate(elements):
      _, Ke, Fe = computeElementContribution(
        element,
        nodes,
        mat,
        nstate,
        Uelas,
        Upf,
        D,
        None if tensile_values is None else tensile_values[element_index],
      )
      _assemble_local_contribution(K, F, element, Ke, Fe, nstate)
    return

  element_chunks = splitElementsAmongThreads(elements, worker_count)
  element_lookup = {
    tuple(element.node_ids): element
    for element in elements
  }
  with ProcessPoolExecutor(
      max_workers=worker_count,
      initializer=initProcessWorker,
      initargs=(nodes, mat, Uelas, Upf, D),
  ) as executor:
    chunk_results = executor.map(
      computeElementChunk,
      ((chunk, nstate) for chunk in element_chunks),
    )
    for local_updates in chunk_results:
      for node_ids, Ke, Fe in local_updates:
        # Elements are returned in input order; retrieve their cached map by
        # node IDs without changing floating-point accumulation order.
        element = element_lookup[tuple(node_ids)]
        _assemble_local_contribution(K, F, element, Ke, Fe, nstate)


# =============================================================================
# ELEMENT EQUATIONS
# =============================================================================

def computeElementStiffness(
    Ke,
    Fe,
    nodes,
    element,
    mat,
    nstate,
    Uelas_data=None,
    Upf_data=None,
    D_mat=None,
    tensile_values=None,
):
  if Uelas_data is None:
    Uelas_data = Uelas
  if Upf_data is None:
    Upf_data = Upf
  if D_mat is None:
    D_mat = D

  integration_data = _get_element_integration_data(element, nodes)
#Here it starts the iterative process between displacement and phase field
  if nstate == elastic_state_size:
    for N_elastic, _, _, B, integration_weight in integration_data:
      phase_field = sum(
        N_elastic[0, elastic_state_size * i]
        * Upf_data[element.node_ids[i]]
        for i in range(len(element.node_ids))
      )
      deteriorated_D = D_mat.copy()
      deteriorated_D *= (1.0 - phase_field) ** 2
      Ke += B.T @ deteriorated_D @ B * integration_weight 
    return

  if nstate == phase_field_state_size:
    Gc_reference, l0 = mat.get_reference_toughness(), mat.l0
    fracture_tensor = mat.get_fracture_tensor()
    c0 = 2.0
    if tensile_values is None:
      tensile_values = _calculate_element_tensile_values(
        element,
        integration_data,
        Uelas_data,
        D_mat,
      )
    for data, sigma_dot_epsilon in zip(integration_data, tensile_values):
      _, N_phase, dN_xy, _, integration_weight = data
      Ke += integration_weight * (
        l0 / c0 * (dN_xy.T @ fracture_tensor @ dN_xy) #dN_xy = grad phi
        + (Gc_reference / (l0 * c0) + 0.5 * sigma_dot_epsilon)
        * N_phase.T @ N_phase
      )
      Fe += (
        integration_weight
        * 0.5
        * sigma_dot_epsilon
        * N_phase.flatten()
      )
    return

  raise ValueError(f"Invalid state size: {nstate}")


def _strain_and_stress_at_point(element, dN, displacement, D_mat):
  displacement_gradient = np.zeros((2, 2))
  for i, node_id in enumerate(element.node_ids):
    index = elastic_state_size * node_id
    displacement_gradient[0, 0] += dN[0, i] * displacement[index]
    displacement_gradient[0, 1] += dN[1, i] * displacement[index]
    displacement_gradient[1, 0] += dN[0, i] * displacement[index + 1]
    displacement_gradient[1, 1] += dN[1, i] * displacement[index + 1]

  strain = 0.5 * (displacement_gradient + displacement_gradient.T)
  strain_vector = np.array([
    strain[0, 0],
    strain[1, 1],
    2.0 * strain[0, 1],
  ])
  stress_vector = D_mat @ strain_vector
  stress = np.array([
    [stress_vector[0], stress_vector[2]],
    [stress_vector[2], stress_vector[1]],
  ])
  return strain, stress


def _calculate_element_tensile_values(
    element,
    integration_data,
    displacement,
    D_mat,
):
  """Evaluate all element quadrature points in two batched eigensolves."""
  point_count = len(integration_data)
  strains = np.empty((point_count, 2, 2))
  stresses = np.empty((point_count, 2, 2))

  for point, data in enumerate(integration_data):
    dN_xy = data[2]
    strains[point], stresses[point] = _strain_and_stress_at_point(
      element,
      dN_xy,
      displacement,
      D_mat,
    )

  # Stress and strain are symmetric tensors. eigvalsh exploits that structure,
  # returns ordered real eigenvalues, and avoids tiny complex round-off parts
  # that cannot be accumulated into the real phase-field matrix.
  stress_eigenvalues = np.linalg.eigvalsh(stresses)
  strain_eigenvalues = np.linalg.eigvalsh(strains)
  return tuple(
    sum(
      stress_eigenvalues[point, i] * strain_eigenvalues[point, i]
      if strain_eigenvalues[point, i] > 0.0 else 0.0
      for i in range(stress_eigenvalues.shape[1])
    )
    for point in range(point_count)
  )


def _calculate_all_tensile_values(elements, nodes, displacement, D_mat):
  """Batch eigensolves for every quadrature point in a serial assembly."""
  element_count = len(elements)
  point_count = len(INTEGRATION_RULE)
  strains = np.empty((element_count, point_count, 2, 2))
  stresses = np.empty((element_count, point_count, 2, 2))

  for element_index, element in enumerate(elements):
    integration_data = _get_element_integration_data(element, nodes)
    for point, data in enumerate(integration_data):
      strains[element_index, point], stresses[element_index, point] = (
        _strain_and_stress_at_point(
          element,
          data[2],
          displacement,
          D_mat,
        )
      )

  stress_eigenvalues = np.linalg.eigvalsh(stresses)
  strain_eigenvalues = np.linalg.eigvalsh(strains)
  first = np.where(
    strain_eigenvalues[:, :, 0] > 0.0,
    stress_eigenvalues[:, :, 0] * strain_eigenvalues[:, :, 0],
    0.0,
  )
  second = np.where(
    strain_eigenvalues[:, :, 1] > 0.0,
    stress_eigenvalues[:, :, 1] * strain_eigenvalues[:, :, 1],
    0.0,
  )
  return first + second


def calculateTensileSigmaDotEps(
    element,
    dN,
    Uelas_data=None,
    D_mat=None,
):
  if Uelas_data is None:
    Uelas_data = Uelas
  if D_mat is None:
    D_mat = D

  strain, stress = _strain_and_stress_at_point(
    element,
    dN,
    Uelas_data,
    D_mat,
  )

  stress_eigenvalues = np.linalg.eigvalsh(stress)
  strain_eigenvalues = np.linalg.eigvalsh(strain)
  return sum(
    stress_eigenvalues[i] * strain_eigenvalues[i]
    if strain_eigenvalues[i] > 0.0 else 0.0
    for i in range(len(stress_eigenvalues))
  )


def computeSigmaAtCenter(
    element,
    nodes,
    stress_vector,
    Uelas_data=None,
    Upf_data=None,
    D_mat=None,
):
  """Compute degraded stress at an element center (legacy rectangle method)."""
  if Uelas_data is None:
    Uelas_data = Uelas
  if Upf_data is None:
    Upf_data = Upf
  if D_mat is None:
    D_mat = D

  qsi = eta = 0.0
  n1, n2, _, n4 = [nodes[i] for i in element.node_ids]
  inverse = np.diag([
    2.0 / (n2.x - n1.x),
    2.0 / (n4.y - n1.y),
  ])
  N, dN = shapeFunctions(qsi, eta, elastic_state_size)
  dN_xy = inverse.T @ dN

  displacement_gradient = np.zeros((2, 2))
  for i, node_id in enumerate(element.node_ids):
    index = elastic_state_size * node_id
    displacement_gradient[0, 0] += dN_xy[0, i] * Uelas_data[index]
    displacement_gradient[0, 1] += dN_xy[1, i] * Uelas_data[index]
    displacement_gradient[1, 0] += dN_xy[0, i] * Uelas_data[index + 1]
    displacement_gradient[1, 1] += dN_xy[1, i] * Uelas_data[index + 1]

  strain = 0.5 * (displacement_gradient + displacement_gradient.T)
  strain_vector = np.array([
    strain[0, 0],
    strain[1, 1],
    2.0 * strain[0, 1],
  ])
  phase_field = sum(
    N[0, elastic_state_size * i] * Upf_data[element.node_ids[i]]
    for i in range(4)
  )
  stress_vector[:] = (1.0 - phase_field) ** 2 * D_mat @ strain_vector
  return phase_field


# =============================================================================
# BOUNDARY CONDITIONS, SOLVER AND OUTPUT
# =============================================================================

penalty_factor = 1.0e10
vtk_extension = ".vtk"

def applyBoundaryConditions(K, F, bc_nodes):
  """Apply the original penalty constraints and nodal Neumann forces."""
  for bc in bc_nodes:
    row = elastic_state_size * bc.node
    x_value = bc.xval * pseudotime
    y_value = bc.yval * pseudotime

    if bc.bctype == 0:
      K[row, row] += penalty_factor
      K[row + 1, row + 1] += penalty_factor
      F[row] += penalty_factor * x_value
      F[row + 1] += penalty_factor * y_value
    elif bc.bctype == 1:
      K[row, row] += penalty_factor
      F[row] += penalty_factor * x_value
    elif bc.bctype == 2:
      K[row + 1, row + 1] += penalty_factor
      F[row + 1] += penalty_factor * y_value
    elif bc.bctype == 3:
      F[row] += x_value
      F[row + 1] += y_value
    else:
      raise ValueError(f"Unknown boundary-condition type: {bc.bctype}")


def solveSystem(K, F, U):
  from scipy import sparse

  U[:] = sparse.linalg.splu(K.tocsc()).solve(F)


def computeReaction(
    K,
    F,
    nodes,
    elements,
    mat,
    target_ids,
    dof_offset=0,
    sign_factor=-1.0,
    nthreads=None,
    parallel_min_elements=default_parallel_min_elements,
):
  assembleGlobalStiffness(
    K,
    F,
    elements,
    nodes,
    mat,
    elastic_state_size,
    nthreads=nthreads,
    parallel_min_elements=parallel_min_elements,
  )
  residual = K.dot(Uelas)
  return sum(
    sign_factor * residual[elastic_state_size * node_id + dof_offset]
    for node_id in target_ids
  )


def generateVTKLegacyFile(
    nodes,
    elements,
    filename,
    Uelas_data=None,
    Upf_data=None,
):
  if Uelas_data is None:
    Uelas_data = Uelas
  if Upf_data is None:
    Upf_data = Upf

  with open(filename, "w", encoding="utf-8") as vtk_file:
    vtk_file.write("# vtk DataFile Version 2.0\n")
    vtk_file.write("FEM results\n")
    vtk_file.write("ASCII\n")
    vtk_file.write("DATASET UNSTRUCTURED_GRID\n")
    vtk_file.write(f"POINTS {len(nodes)} float\n")
    vtk_file.writelines(
      f"{node.x} {node.y} 0.0\n"
      for node in nodes
    )
    vtk_file.write(f"CELLS {len(elements)} {len(elements) * 5}\n")
    vtk_file.writelines(
      f"4 {' '.join(map(str, element.node_ids))}\n"
      for element in elements
    )
    vtk_file.write(f"CELL_TYPES {len(elements)}\n")
    vtk_file.writelines("9\n" for _ in elements)
    vtk_file.write(f"POINT_DATA {len(nodes)}\n")
    vtk_file.write("VECTORS displacements float\n")
    vtk_file.writelines(
      f"{Uelas_data[2 * i]} {Uelas_data[2 * i + 1]} 0.0\n"
      for i in range(len(nodes))
    )
    vtk_file.write("SCALARS phasefield float 1\n")
    vtk_file.write("LOOKUP_TABLE default\n")
    vtk_file.writelines(f"{value}\n" for value in Upf_data)


def create_bc_nodes_from_config(
    boundary_conditions,
    physical_groups=None,
    physical_name_to_tag=None,
):
  """Expand named or numeric Gmsh physical groups into BC objects.

  When a condition has ``bc_type == 3`` (Neumann) and ``distribute_total``
  is set, ``xval``/``yval`` are treated as the TOTAL force over all nodes
  matched by ``node_filter`` and are divided evenly among them.
  """
  bc_nodes = []
  for condition in boundary_conditions:
    node_filter = condition.node_filter
    if not isinstance(node_filter, (str, int)):
      raise TypeError(
        f"Boundary condition '{condition.name}' node_filter must be a "
        "physical-group name (str) or tag (int)."
      )
    if physical_groups is None:
      raise ValueError(
        "Physical groups must be loaded before resolving boundary conditions."
      )
    if isinstance(node_filter, str):
      if (
          not physical_name_to_tag
          or node_filter not in physical_name_to_tag
      ):
        available = (
          sorted(physical_name_to_tag)
          if physical_name_to_tag else []
        )
        raise ValueError(
          f"Boundary condition '{condition.name}' references physical "
          f"group '{node_filter}', which was not found in the mesh. "
          f"Available names: {available}"
        )
      node_filter = physical_name_to_tag[node_filter]
    node_ids = physical_groups.get(node_filter, [])

    xval, yval = condition.xval, condition.yval
    if condition.bc_type == 3 and getattr(condition, "distribute_total", False):
      if not node_ids:
        raise ValueError(
          f"Force boundary condition '{condition.name}' found no nodes "
          f"in physical group {node_filter!r}."
        )
      xval = xval / len(node_ids)
      yval = yval / len(node_ids)

    bc_nodes.extend(
      BC(
        node_id,
        condition.bc_type,
        xval,
        yval,
      )
      for node_id in node_ids
    )
  return bc_nodes


def resolve_node_ids(
    nodes,
    node_filter,
    physical_groups=None,
    physical_name_to_tag=None,
):
  del nodes  # Kept for compatibility with the previous internal signature.
  if not isinstance(node_filter, (str, int)):
    raise TypeError(
      "node_filter must be a physical-group name (str) or tag (int)."
    )
  if physical_groups is None:
    return []
  if isinstance(node_filter, str):
    if not physical_name_to_tag or node_filter not in physical_name_to_tag:
      return []
    node_filter = physical_name_to_tag[node_filter]
  return list(physical_groups.get(node_filter, []))


def resolve_reaction_target_ids(
    physical_groups,
    physical_name_to_tag,
    reaction_config,
):
  """Resolve reaction-monitor nodes from current and legacy group names."""
  if not reaction_config:
    return []

  reaction_type = (reaction_config.reaction_type or "").strip()
  reaction_type_lower = reaction_type.lower()
  name_lookup = {
    name.lower(): name
    for name in (physical_name_to_tag or {})
  }

  def ids_for_name(group_name):
    actual_name = name_lookup.get(group_name.lower())
    if actual_name is None:
      return []
    tag = physical_name_to_tag[actual_name]
    return list(physical_groups.get(tag, []))

  base_name = (
    reaction_type.rsplit("_", 1)[0]
    if reaction_type_lower.endswith(("_x", "_y"))
    else reaction_type
  )
  target_ids = ids_for_name(base_name)
  if target_ids:
    return target_ids

  aliases = (
    ("bottom", "Bottom_ids"),
    ("top", "Top_ids"),
    ("left", "Left_ids"),
    ("right", "Right_ids"),
  )
  for keyword, group_name in aliases:
    if keyword in reaction_type_lower:
      return ids_for_name(group_name)

  if "support" in reaction_type_lower:
    return ids_for_name("Left_ids") + ids_for_name("Right_ids")
  return []


def get_loaded_surface_node_ids(bc_nodes):
  """Return unique nodes carrying Neumann loads, preserving input order."""
  return list(dict.fromkeys(
    bc.node for bc in bc_nodes if bc.bctype == 3
  ))


def computeMonitoredDisplacement(
    node_ids,
    dof,
    Uelas,
    aggregation="average",
):
  if not node_ids:
    return 0.0
  values = [
    Uelas[elastic_state_size * node_id + dof]
    for node_id in node_ids
  ]
  if aggregation == "max":
    return max(values, key=abs)
  return sum(values) / len(values)


def detectLoadingType(config):
  """Return ``force`` or ``displacement`` from the active loading BCs.

  Zero-valued Dirichlet conditions are supports and do not characterize the
  loading type. An explicit ``loading_type`` still overrides auto-detection,
  which is useful for intentionally mixed loading cases.
  """
  configured_type = getattr(config, "loading_type", "auto")
  if configured_type in ("force", "displacement"):
    return configured_type
  if configured_type not in (None, "auto"):
    raise ValueError(
      f"Unknown loading_type '{configured_type}'. Use 'auto', "
      "'force' or 'displacement'."
    )

  has_force_loading = False
  has_displacement_loading = False
  for condition in config.boundary_conditions:
    if condition.bc_type == 3:
      has_force_loading |= (
        condition.xval != 0.0 or condition.yval != 0.0
      )
    elif condition.bc_type == 0:
      has_displacement_loading |= (
        condition.xval != 0.0 or condition.yval != 0.0
      )
    elif condition.bc_type == 1:
      has_displacement_loading |= condition.xval != 0.0
    elif condition.bc_type == 2:
      has_displacement_loading |= condition.yval != 0.0

  if has_force_loading and not has_displacement_loading:
    return "force"
  if has_displacement_loading and not has_force_loading:
    return "displacement"
  if has_force_loading and has_displacement_loading:
    raise ValueError(
      "Automatic loading detection found both a non-zero Neumann force and "
      "a non-zero prescribed displacement. Set loading_type explicitly to "
      "'force' or 'displacement' to select the monitoring convention."
    )
  raise ValueError(
    "Automatic loading detection found no non-zero force or prescribed "
    "displacement in boundary_conditions."
  )


# =============================================================================
# MAIN SIMULATION
# =============================================================================

def _load_configuration(config_name):
  from config_simulations import get_simulation_config

  print(f"\n Loading configuration: {config_name}")
  config = get_simulation_config(config_name)
  print(f" {config.name}\n")

  if config.mesh_type != "gmsh":
    raise ValueError(
      f"Mesh type '{config.mesh_type}' not supported. Please use 'gmsh'."
    )

  mesh_file = os.path.join("GmshFiles", config.mesh_file)
  print(f" Used mesh: {mesh_file}")
  (
    nodes,
    elements,
    physical_groups,
    physical_name_to_tag,
  ) = readGmshMesh(mesh_file)
  return (
    config,
    mesh_file,
    nodes,
    elements,
    physical_groups,
    physical_name_to_tag,
  )


def _material_from_config(material_config):
  """Build the solver material from the unified orthotropic fields."""
  required = ("E1", "E2", "nu12", "G12", "Gc1", "Gc2", "l0")
  missing = [
    name for name in required
    if not hasattr(material_config, name)
  ]
  if missing:
    raise ValueError(
      "Material configuration must define E1, E2, nu12, G12, Gc1, Gc2 and "
      f"l0. Missing: {missing}"
    )

  return MaterialParameters(
    E1=material_config.E1,
    E2=material_config.E2,
    nu12=material_config.nu12,
    G12=material_config.G12,
    Gc1=material_config.Gc1,
    Gc2=material_config.Gc2,
    l0=material_config.l0,
    material_type=getattr(
      material_config,
      "material_type",
      "plane_stress",
    ),
    theta_deg=getattr(material_config, "theta_deg", 0.0),
  )


def _print_simplebar_reference(config, mesh_file, material):
  if "simplebar.msh" not in mesh_file:
    return None

  if (
      not material.is_isotropic_equivalent()
      or not material.is_fracture_isotropic()
  ):
    print(
      "Note: analytical sigma/u peak formula is isotropic-only; "
      "skipped for this orthotropic material.\n"
    )
    return None

  sigma_peak = np.sqrt(
    27.0 * material.E1 * material.Gc1 / (256.0 * material.l0)
  )
  u_peak = (
    16.0 / 9.0
    * sigma_peak
    * config.material.length
    / material.E1
  )
  print(f"Sigma peak: {sigma_peak}")
  print(f"U peak: {u_peak}\n")
  return sigma_peak


def _collect_graph_data(
    config,
    loading_type,
    elements,
    nodes,
    physical_groups,
    physical_name_to_tag,
    bc_nodes,
    Kelas,
    Felas,
    material,
    sigma_peak_at2,
    u_data,
    force_data,
    stress_data,
    time_data,
    nthreads,
    parallel_min_elements,
):
  if not config.graph_config:
    return

  graph = config.graph_config
  if graph.graph_type == "stress_vs_time":
    element_id = 450 if graph.element_id is None else graph.element_id
    if element_id < len(elements):
      stress = np.zeros(3)
      computeSigmaAtCenter(elements[element_id], nodes, stress)
      component = getattr(graph, "stress_component", 0)
      stress_data.append(
        stress[component] / sigma_peak_at2
        if sigma_peak_at2 is not None
        else stress[component]
      )
      time_data.append(pseudotime)
    return

  if graph.graph_type != "force_vs_displacement":
    return

  dof = 0 if graph.displacement_axis == "x" else 1
  if loading_type == "force":
    applied_force = sum(
      (bc.xval if dof == 0 else bc.yval) * pseudotime
      for bc in bc_nodes
      if bc.bctype == 3
    )
    force_data.append(applied_force)

    monitor = getattr(config, "displacement_monitor", None)
    if monitor and monitor.node_filter is not None:
      monitor_ids = resolve_node_ids(
        nodes,
        monitor.node_filter,
        physical_groups,
        physical_name_to_tag,
      )
      monitor_dof = monitor.dof
      aggregation = monitor.aggregation
    else:
      monitor_ids = get_loaded_surface_node_ids(bc_nodes)
      monitor_dof = dof
      aggregation = monitor.aggregation if monitor else "average"

    u_data.append(computeMonitoredDisplacement(
      monitor_ids,
      monitor_dof,
      Uelas,
      aggregation,
    ))
    return

  reaction_config = config.reaction_config
  target_ids = resolve_reaction_target_ids(
    physical_groups,
    physical_name_to_tag,
    reaction_config,
  )
  dof_offset = 0
  sign_factor = -1.0
  if reaction_config:
    dof_offset = reaction_config.reaction_dof
    sign_factor = reaction_config.sign_factor

  force_data.append(computeReaction(
    Kelas,
    Felas,
    nodes,
    elements,
    material,
    target_ids,
    dof_offset,
    sign_factor,
    nthreads,
    parallel_min_elements,
  ))

  displacement_value = 0.0
  if graph.displacement_axis == "x":
    for condition in config.boundary_conditions:
      if condition.xval != 0:
        displacement_value = condition.xval
        break
  else:
    for condition in config.boundary_conditions:
      if condition.yval != 0:
        displacement_value = abs(condition.yval)
        break
  u_data.append(pseudotime * displacement_value)


def _save_graph(config, u_data, force_data, stress_data, time_data):
  graph = config.graph_config
  if not graph:
    return

  # Import lazily so worker processes and simulations without graphs avoid the
  # comparatively expensive Matplotlib startup.
  import matplotlib.pyplot as plt

  figure, axes = plt.subplots()
  if graph.graph_type == "stress_vs_time":
    axes.plot(time_data, stress_data, "o")
  elif graph.graph_type == "force_vs_displacement":
    axes.plot(u_data, force_data, "o")

  axes.set_xlabel(graph.x_label)
  axes.set_ylabel(graph.y_label)
  axes.set_title(graph.title)
  axes.grid(True)
  figure.savefig(graph.output_file)
  plt.close(figure)
  print(f"✓ Graph saved to: {graph.output_file}")


def _final_monitored_node_ids(
    config,
    loading_type,
    nodes,
    bc_nodes,
    physical_groups,
    physical_name_to_tag,
):
  if loading_type == "force":
    monitor = getattr(config, "displacement_monitor", None)
    if monitor and monitor.node_filter is not None:
      return resolve_node_ids(
        nodes,
        monitor.node_filter,
        physical_groups,
        physical_name_to_tag,
      )
    return get_loaded_surface_node_ids(bc_nodes)
  return list(dict.fromkeys(
    bc.node
    for bc in bc_nodes
    if (
      (bc.bctype in (0, 1) and bc.xval != 0.0)
      or (bc.bctype in (0, 2) and bc.yval != 0.0)
    )
  ))


def main(config_name="default"):
  from scipy import sparse

  os.makedirs("outputs", exist_ok=True)

  try:
    (
      config,
      mesh_file,
      nodes,
      elements,
      physical_groups,
      physical_name_to_tag,
    ) = _load_configuration(config_name)
  except Exception as error:
    print(f"\n Error loading configuration: {error}")
    print("Could not read simulation parameters.")
    print(f"Requested configuration: {config_name}")
    print("Please check the configuration and try again.\n")
    return

  simulation_time = Timer()
  material_config = config.material
  simulation = config.simulation_params
  material = _material_from_config(material_config)

  sigma_peak_at2 = _print_simplebar_reference(
    config,
    mesh_file,
    material,
  )
  loading_type = detectLoadingType(config)
  automatically_detected = getattr(config, "loading_type", "auto") in (
    None,
    "auto",
  )
  if loading_type == "displacement":
    label = " detected automatically" if automatically_detected else ""
    print(f"Loading type{label}: imposed displacement (Dirichlet BC).")
  else:
    label = " detected automatically" if automatically_detected else ""
    print(
      f"Loading type{label}: force-controlled (Neumann BCs); "
      "displacement is a result of the solve."
    )
  print(
    f"Material: E1={material.E1}, E2={material.E2}, "
    f"nu12={material.nu12}, G12={material.G12}, "
    f"theta={material.theta_deg} deg, "
    f"Gc1={material.Gc1}, Gc2={material.Gc2}"
  )

  global D, Uelas, Upf, pseudotime
  D = material.get_D_matrix()

  print(
    f" Applying {len(config.boundary_conditions)} boundary conditions "
    "from the configuration..."
  )
  bc_nodes = create_bc_nodes_from_config(
    config.boundary_conditions,
    physical_groups,
    physical_name_to_tag,
  )

  # Geometry and sparse insertion positions are invariant during staggered
  # iterations and are therefore computed once.
  prepareElementIntegrationData(nodes, elements)

  number_of_nodes = len(nodes)
  elastic_dofs = elastic_state_size * number_of_nodes
  phase_field_dofs = phase_field_state_size * number_of_nodes

  Kelas = sparse.lil_matrix((elastic_dofs, elastic_dofs))
  createSparseStructure(Kelas, elements, elastic_state_size)
  Felas = np.zeros(elastic_dofs)
  Uelas = np.zeros(elastic_dofs)

  Kpf = sparse.lil_matrix((phase_field_dofs, phase_field_dofs))
  createSparseStructure(Kpf, elements, phase_field_state_size)
  Fpf = np.zeros(phase_field_dofs)
  Upf = np.zeros(phase_field_dofs)

  u_data = []
  force_data = []
  stress_data = []
  time_data = []

  dt = simulation["dt"]
  total_time = simulation["totaltime"]
  maximum_steps = int(simulation["maxsteps"])
  maximum_iterations = int(simulation["maxiter"])
  staggered_tolerance = simulation["stagtol"]
  nthreads = simulation.get("nthreads", DEFAULT_NTHREADS)
  parallel_min_elements = simulation.get(
    "parallel_min_elements",
    default_parallel_min_elements,
  )
  if maximum_iterations < 1:
    raise ValueError("simulation_params['maxiter'] must be at least 1.")
  if parallel_min_elements < 0:
    raise ValueError(
      "simulation_params['parallel_min_elements'] cannot be negative."
    )
  elastic_threads = _effective_worker_count(
    nthreads,
    len(elements),
    parallel_min_elements * elastic_parallel_threshold_multiplier,
  )
  phase_field_threads = _effective_worker_count(
    nthreads,
    len(elements),
    parallel_min_elements,
  )
  print(
    f" Using {elastic_threads} elastic / {phase_field_threads} phase-field "
    f"worker(s); {nthreads} requested, "
    f"CPUs detected: {os.cpu_count()}"
  )

  pseudotime = 0.0
  for step in range(maximum_steps):
    pseudotime += dt
    if pseudotime > total_time:
      break

    print(
      f"******************** Time Step {step} | "
      f"Pseudo time = {pseudotime:.6f} | Time step = {dt} "
      "********************"
    )

    converged = False
    for iteration in range(maximum_iterations):
      print(f"------ Staggered Iteration {iteration} ------")
      assembleGlobalStiffness(
        Kelas,
        Felas,
        elements,
        nodes,
        material,
        elastic_state_size,
        nthreads,
        parallel_min_elements,
      )
      applyBoundaryConditions(Kelas, Felas, bc_nodes)

      if iteration != 0:
        residual = Kelas.tocsr() @ Uelas - Felas
        residual_norm = np.linalg.norm(residual)
        print(f"Residual Elasticity Norm: {residual_norm:.2e}")
        if residual_norm < staggered_tolerance:
          print(
            "------> Staggered scheme converged in "
            f"{iteration} iterations."
          )
          converged = True
          break

      solveSystem(Kelas, Felas, Uelas)
      assembleGlobalStiffness(
        Kpf,
        Fpf,
        elements,
        nodes,
        material,
        phase_field_state_size,
        nthreads,
        parallel_min_elements,
      )
      solveSystem(Kpf, Fpf, Upf)

    if maximum_iterations > 0 and not converged:
      print(
        "------> Staggered scheme reached the configured iteration limit.\n"
        "Accepting current solution and continuing"
      )

    output_filename = f"{config.output_base}{step}{vtk_extension}"
    generateVTKLegacyFile(
      nodes,
      elements,
      output_filename,
      Uelas,
      Upf,
    )

    _collect_graph_data(
      config,
      loading_type,
      elements,
      nodes,
      physical_groups,
      physical_name_to_tag,
      bc_nodes,
      Kelas,
      Felas,
      material,
      sigma_peak_at2,
      u_data,
      force_data,
      stress_data,
      time_data,
      nthreads,
      parallel_min_elements,
    )

  _save_graph(config, u_data, force_data, stress_data, time_data)

  monitored_ids = _final_monitored_node_ids(
    config,
    loading_type,
    nodes,
    bc_nodes,
    physical_groups,
    physical_name_to_tag,
  )
  final_ux = computeMonitoredDisplacement(monitored_ids, 0, Uelas)
  final_uy = computeMonitoredDisplacement(monitored_ids, 1, Uelas)
  print("\nFinal mean displacements on loaded surface:")
  print(f"  X: {final_ux:.6e}")
  print(f"  Y: {final_uy:.6e}")
  print(f"  Nodes monitored: {len(monitored_ids)}")

  print("\n================> Simulation completed!")
  simulation_time.elapsed("complete simulation")


if __name__ == "__main__":
  selected_config = sys.argv[1] if len(sys.argv) > 1 else "default"
  print("\n" + "=" * 80)
  print(f"PhaseFieldJr - Simulation: {selected_config}")
  print("=" * 80 + "\n")
  main(selected_config)
