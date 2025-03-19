# MIT License

# © 2025 Nathan Shauer

# phasefield-jr

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
from numpy import dot
import time
import matplotlib.pyplot as plt
import os
from scipy import sparse

# =============================== DATA STRUCTURES ===============================
# ===============================================================================
class QuadraturePoint:
  def __init__(self, xi, eta, weight):
    self.xi = xi
    self.eta = eta
    self.weight = weight

class MaterialParameters:
  def __init__(self, E, nu, G, l):
    self.E = E  # Young's modulus
    self.nu = nu  # Poisson's ratio
    self.G = G  # Critical strain energy release rate
    self.l = l  # Length scale parameter

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
    self.type = bc_type  # 0 dirichlet in x and y, 1 dirichlet in x, 2 dirichlet in y, 3 neumann
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
# This is not ideal, but since it is a simple example there is no problem
global Uelas, Upf, D, pseudotime, basefilename, vtkextension, intrule
Uelas = np.zeros(1)  # Vector with one position valuing zero
Upf = np.zeros(1)  # Initialize with the correct size
D = np.zeros((3, 3)) # Assuming all elements have same materials
pseudotime = 0.0 # Is updated to impose displacement gradually
basefilename = "outputs/output_ex2_"
vtkextension = ".vtk"
intrule = create2x2QuadratureRule()  # Adopting 2x2 quadrature rule

# =============================== FUNCTION IMPLEMENTATIONS ======================
# ===============================================================================
def createDoubleNodeMesh(nodes, elements, num_elements_x, num_elements_y, length, height):
  xsize = ysize = 0.001
  num_elements_y //= 2
  num_elements_x_small = num_elements_y_small = 8
  num_elements_y_large = num_elements_y - num_elements_y_small
  num_elements_x_large = num_elements_x - num_elements_x_small
  y_small = ysize * num_elements_y_small
  y_largeel_size = (height / 2 - y_small) / num_elements_y_large
  x_small = xsize * num_elements_x_small
  x_largeel_size = (length - x_small) / num_elements_x_large

  # Generate nodes
  if num_elements_x % 2 != 0:
    raise ValueError("num_elements_x must be an even number")
  ndoublenodes = int(num_elements_x_large / 2 + num_elements_x_small / 2)
  nodes.resize((num_elements_x+1) * (2*num_elements_y+1) + ndoublenodes, refcheck=False)

  for j in range(num_elements_y + 1):
    xnow = -0.5
    for i in range(num_elements_x + 1):
      if j <= num_elements_y_small:
        y = j * ysize
      else:
        y = y_small + (j - num_elements_y_small) * y_largeel_size
      nodes[j * (num_elements_x + 1) + i] = Node(xnow, y)
      if i < num_elements_x_large/2 or i > num_elements_x_large/2 + num_elements_x_small-1:
        xnow += x_largeel_size
      else:
        xnow += xsize      

  # Generate elements
  elements.resize(num_elements_x * 2*num_elements_y, refcheck=False)
  for j in range(num_elements_y):
    for i in range(num_elements_x):
      n1 = j * (num_elements_x + 1) + i
      n2 = n1 + 1
      n3 = n1 + num_elements_x + 1
      n4 = n3 + 1
      elements[j * num_elements_x + i] = Element([n1, n2, n4, n3])

  # Mirror the mesh to negative y
  original_node_count = (num_elements_x+1) * (num_elements_y+1)
  node_map = {}

  count = original_node_count;
  for i in range(original_node_count):
    if abs(nodes[i].y) > 1.e-8:      
      mirrored_node = Node(nodes[i].x, -nodes[i].y)
      node_map[i] = count
      nodes[count] = mirrored_node
      count = count + 1
    else:
      if nodes[i].x < -1.e-8: # nodes until the tip
        duplicated_node = Node(nodes[i].x, nodes[i].y)
        node_map[i] = count
        nodes[count] = duplicated_node
        count = count + 1
      else:
        node_map[i] = i
  
  # Generate elements for the mirrored part
  original_element_count = num_elements_x * num_elements_y
  count = original_element_count
  for i in range(original_element_count):
    mirrored_element = Element([node_map[node_id] for node_id in elements[i].node_ids])
    mirrored_element.node_ids.reverse()
    elements[count] = mirrored_element
    count = count + 1    

def createSparseStructure(K, elements, nodes, nstate):
  for element in elements:
    for i in range(4):
      row = nstate * element.node_ids[i]
      for j in range(4):
        col = nstate * element.node_ids[j]
        for k in range(nstate):
          for l in range(nstate):            
            K[row + k, col + l] = 1. # we need a number to create sparse structure, later we zero it out

def assembleGlobalStiffness(K, F, elements, nodes, mat, nstate):
  time = Timer()
  for row in K.data:
    for i in range(len(row)):
      row[i] = 0.0
  F.fill(0)
  for element in elements:
    nquadnodes = 4
    ndofel = nstate * nquadnodes
    Ke = np.zeros((ndofel, ndofel))
    Fe = np.zeros(ndofel)
    computeElementStiffness(Ke, Fe, nodes, element, mat, nstate)
    for i in range(nquadnodes):
      row = nstate * element.node_ids[i]
      for k in range(nstate):
        F[row + k] += Fe[nstate * i + k]
      for j in range(nquadnodes):
        col = nstate * element.node_ids[j]
        for k in range(nstate):
          for l in range(nstate):
            index = K.rows[row+k].index(col+l)
            K.data[row+k][index] += Ke[nstate * i + k, nstate * j + l]
  # time.elapsed("assembly")

def computeReaction(K, F, nodes, elements, mat):
  assembleGlobalStiffness(K, F, elements, nodes, mat, 2)
  residual = K.dot(Uelas)  # F is zero
  reaction = 0.0
  for i in range(len(nodes)):
    if abs(nodes[i].y + 0.5) < 1.e-8:
      reaction += -residual[2 * i + 1]
  return reaction
  
def computeElementStiffness(Ke, Fe, nodes, element, mat, nstate):
  n1, n2, n3, n4 = [nodes[i] for i in element.node_ids]
  base = n2.x - n1.x
  height = n4.y - n1.y
  area = base * height
  detjac = area / 4.0
  dqsidx = 2.0 / base
  dqsidy = 2.0 / height
  J_inv = np.diag([dqsidx, dqsidy])

  if nstate == 2: # compute elasticity stiffness
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      dN_xy = J_inv.T @ dN.T
      B = createB(dN_xy.T)
      Ddeteriorated = D.copy()
      phase_field = sum(N[0, 2 * i] * Upf[element.node_ids[i]] for i in range(4))
      Ddeteriorated *= (1 - phase_field) ** 2
      Ke += B.T @ Ddeteriorated @ B * qp.weight * detjac
  elif nstate == 1: # compute phase field stiffness
    G, l = mat.G, mat.l
    c0 = 2.0
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      dN_xy = J_inv.T @ dN.T # Same as B_phi
      sigmaDotEps = calculateSigmaDotEps(element, dN_xy)
      Ke += detjac * qp.weight * (G * l / c0 * (dN_xy.T @ dN_xy) + (G / (l * c0) + 0.5 * sigmaDotEps) * N.T @ N)
      Fe += detjac * qp.weight * 0.5 * sigmaDotEps * N.flatten()
  else:
    raise Exception("Invalid state")

def calculateSigmaDotEps(element, dN):
  dU = np.zeros((2, 2))
  for i in range(4):
    index = 2 * element.node_ids[i]
    dU[0, 0] += dN[0, i] * Uelas[index]
    dU[0, 1] += dN[1, i] * Uelas[index]
    dU[1, 0] += dN[0, i] * Uelas[index + 1]
    dU[1, 1] += dN[1, i] * Uelas[index + 1]
  strain = 0.5 * (dU + dU.T)
  strain_vec = np.array([strain[0, 0], strain[1, 1], 2 * strain[0, 1]])
  stress_vec = D @ strain_vec
  sigmaDotEps = stress_vec @ strain_vec
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
    [0.25 * (-1 + eta), 0.25 * (-1 + qsi)],
    [0.25 * (1 - eta), 0.25 * (-1 - qsi)],
    [0.25 * (1 + eta), 0.25 * (1 + qsi)],
    [0.25 * (-1 - eta), 0.25 * (1 - qsi)]
  ])
  return N, dN

def createB(dN):
  B = np.zeros((3, 8))
  for i in range(4):
    B[0, 2 * i] = dN[i, 0]
    B[1, 2 * i + 1] = dN[i, 1]
    B[2, 2 * i] = dN[i, 1]
    B[2, 2 * i + 1] = dN[i, 0]
  return B

def zeroRowAndColumnOfSparseMatrix(K, target_row):
  K.data[target_row] = [0] * len(K.data[target_row])
  for i in range(K.shape[0]):    
    if target_row in K.rows[i]:
        col_idx = K.rows[i].index(target_row)
        K.data[i][col_idx] = 0

def applyBoundaryConditions(K, F, bc_nodes):
  for bc in bc_nodes:
    row = 2 * bc.node
    xval = bc.xval * pseudotime
    yval = bc.yval * pseudotime
    if bc.type == 0:
      F -= K[row,:].toarray().flatten() * xval
      F -= K[row + 1,:].toarray().flatten() * yval
      zeroRowAndColumnOfSparseMatrix(K, row)
      zeroRowAndColumnOfSparseMatrix(K, row+1)
      K[row, row] = 1.0
      K[row + 1, row + 1] = 1.0
      F[row] = xval
      F[row + 1] = yval
    elif bc.type == 1:
      F -= K[row,:].toarray().flatten() *  xval
      zeroRowAndColumnOfSparseMatrix(K, row)
      K[row, row] = 1.0
      F[row] = xval
    elif bc.type == 2:
      F -= K[row + 1,:].toarray().flatten() * yval
      zeroRowAndColumnOfSparseMatrix(K, row+1)
      K[row + 1, row + 1] = 1.0
      F[row + 1] = yval
    elif bc.type == 3:
      F[row] += xval
      F[row + 1] += yval

def solveSystem(K, F, U):
  time = Timer()
  U[:] = sparse.linalg.splu(K.tocsc()).solve(F)
  # time.elapsed("solve")

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

# =============================== MAIN ==========================================
# ===============================================================================
def main():
  # Create folder outputs if it does not exist
  if not os.path.exists("outputs"):
    os.makedirs("outputs")

  simulation_time = Timer()
  E = 210.0  # Young's modulus in Pascals
  nu = 0.3  # Poisson's ratio
  G = 2.7e-3  # Strain energy release rate
  l = 0.005  # Length scale parameter

  # Define mesh and time step parameters
  num_elements_x = 50 # has to be even number
  num_elements_y = 30  # has to be even number
  length = 1.0
  height = 1.0
  dt = 0.01
  totaltime = 0.8
  maxsteps = int(1e5)  # maximum number of time steps (in case using adaptive time step)
  maxiter = 1000  # maximum number of iterations for the staggered scheme
  stagtol = 1e-7  # tolerance to consider the staggered scheme converged

  # Boundary conditions
  imposed_displacement_y = 0.01  # such that we have nucleation at step 50

  nodes = np.array([], dtype=object)
  elements = np.array([], dtype=object)
  bc_nodes = np.array([], dtype=object)

  material = MaterialParameters(E, nu, G, l)
  factor = E / (1 - nu * nu)
  global D
  D = np.zeros((3, 3))
  D[0, 0] = factor
  D[0, 1] = factor * nu
  D[1, 0] = factor * nu
  D[1, 1] = factor
  D[2, 2] = factor * (1 - nu) / 2.0

  createDoubleNodeMesh(nodes, elements, num_elements_x, num_elements_y, length, height)
  
  firstnode = True
  for i in range(len(nodes)):
    if abs(nodes[i].y + 0.5) < 1.e-8:
      if firstnode:
        bc_nodes = np.append(bc_nodes, BC(i, 0, 0.0, 0.0))  # Fix x and y displacement
        firstnode = False
      else:
        bc_nodes = np.append(bc_nodes, BC(i, 2, 0.0, 0.0))  # Fix y displacement
    elif abs(nodes[i].y - 0.5) < 1.e-8:
      bc_nodes = np.append(bc_nodes, BC(i, 2, 0.0, imposed_displacement_y))  # Impose total y displacement on the top edge

  nstate_elas = 2
  nstate_pf = 1
  nnodes = len(nodes)
  ndofs_elas = nstate_elas * nnodes
  ndofs_pf = nstate_pf * nnodes
  Kelas =  sparse.lil_matrix((ndofs_elas, ndofs_elas))
  createSparseStructure(Kelas, elements, nodes, nstate_elas)
  Felas = np.zeros(ndofs_elas)
  global Uelas
  Uelas = np.zeros(ndofs_elas)

  Kpf = sparse.lil_matrix((ndofs_pf, ndofs_pf))
  createSparseStructure(Kpf, elements, nodes, nstate_pf)
  Fpf = np.zeros(ndofs_pf)  
  global Upf
  Upf = np.zeros(ndofs_pf)

  # Data structure to save the data
  u_data = []
  force_data = []

  global pseudotime
  pseudotime = 0.0
  for step in range(maxsteps):
    pseudotime += dt
    if pseudotime > totaltime:
      break
    print(f"******************** Time Step {step} | Pseudo time = {pseudotime:.6f} | Time step = {dt} ********************")
    for iter in range(maxiter):      
      print(f"------ Staggered Iteration {iter} ------")
      assembleGlobalStiffness(Kelas, Felas, elements, nodes, material, nstate_elas)
      applyBoundaryConditions(Kelas, Felas, bc_nodes)
      if iter != 0:
        residual = Kelas.tocsr().dot(Uelas) - Felas
        norm = np.linalg.norm(residual)
        print(f"Residual Elasticity Norm: {norm:.2e}")
        if norm < stagtol:
          print(f"------> Staggered scheme converged in {iter} iterations.")
          break
      solveSystem(Kelas, Felas, Uelas)
      assembleGlobalStiffness(Kpf, Fpf, elements, nodes, material, nstate_pf)
      solveSystem(Kpf, Fpf, Upf)
    if iter == maxiter:
      print(f"------> Staggered scheme did not converge in {maxiter} iterations.\nAccepting current solution and continuing")
    filename = f"{basefilename}{step}{vtkextension}"
    generateVTKLegacyFile(nodes, elements, filename)

    # Save the data
    reaction = computeReaction(Kelas, Felas, nodes, elements, material)
    u_data.append(pseudotime*imposed_displacement_y)
    force_data.append(reaction)

  # Plot the data using matplotlib
  plt.figure()
  plt.plot(u_data, force_data, 'o')
  plt.xlabel('u (mm)')
  plt.ylabel('Force (kN)')
  plt.title('Force vs imposed u')
  plt.grid(True)
  plt.savefig('outputs/force_vs_u.png')
  # plt.show()

  print("\n================> Simulation completed!")
  simulation_time.elapsed("complete simulation")

if __name__ == "__main__":
  main()
