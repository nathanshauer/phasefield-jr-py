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
import time
import matplotlib.pyplot as plt
import os

# =============================== DATA STRUCTURES ===============================
# ===============================================================================
class QuadraturePoint:
  def __init__(self, xi, eta, weight):
    self.xi = xi
    self.eta = eta
    self.weight = weight

class MaterialParameters:
  def __init__(self, E, nu, Gc, l0):
    self.E = E  # Young's modulus
    self.nu = nu  # Poisson's ratio
    self.Gc = Gc  # Critical strain energy release rate
    self.l0 = l0  # Length scale parameter

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
# This is not ideal, but since it is a simple example there is no problem
global Uelas, Upf, D, pseudotime, basefilename, vtkextension, intrule
Uelas = np.zeros(1) # Global vector with nodal values for the displacement approximation
Upf = np.zeros(1) # Global vector with nodal values for the phase field approximation
D = np.zeros((3, 3)) # Constitutive matrix of the 2D elasticity problem. 
pseudotime = 0.0 # Pseudo time used to control incremental displacement/load steps
basefilename = "outputs/output_ex1_" # Base name for the Paraview output files
vtkextension = ".vtk" # Extension for the Paraview output files
intrule = create2x2QuadratureRule() # Integration rule. Adopting 2x2 quadrature rule

# =============================== FUNCTION IMPLEMENTATIONS ======================
# ===============================================================================
def createRectangularMesh(nodes, elements, num_elements_x, num_elements_y, length, height):
  nodes.resize((num_elements_x+1) * (num_elements_y+1), refcheck=False)
  for j in range(num_elements_y+1):
    for i in range(num_elements_x+1):
      nodes[j * (num_elements_x+1) + i] = Node(i * length / num_elements_x, j * height / num_elements_y)

  elements.resize(num_elements_x * num_elements_y, refcheck=False)
  for j in range(num_elements_y):
    for i in range(num_elements_x):
      n1 = j * (num_elements_x + 1) + i
      n2 = n1 + 1
      n3 = n1 + num_elements_x + 1
      n4 = n3 + 1
      elements[j * num_elements_x + i] = Element([n1, n2, n4, n3])

def assembleGlobalStiffness(K, F, elements, nodes, mat, nstate):
  time = Timer()
  K.fill(0)
  F.fill(0)
  for element in elements:
    nnodesel = len(element.node_ids)
    ndofel = nstate * nnodesel
    Ke = np.zeros((ndofel, ndofel))
    Fe = np.zeros(ndofel)
    computeElementStiffness(Ke, Fe, nodes, element, mat, nstate)
    for i in range(nnodesel):
      row = nstate * element.node_ids[i]
      for k in range(nstate):
        F[row + k] += Fe[nstate * i + k]
      for j in range(nnodesel):
        col = nstate * element.node_ids[j]
        for k in range(nstate):
          for l in range(nstate):
            K[row + k, col + l] += Ke[nstate * i + k, nstate * j + l]
  # time.elapsed("assembly")

def computeElementStiffness(Ke, Fe, nodes, element, mat, nstate):
  nnodes = len(element.node_ids)
  n1, n2, n3, n4 = [nodes[i] for i in element.node_ids]
  base = n2.x - n1.x
  height = n4.y - n1.y
  area = base * height
  detjac = area / 4.0
  dqsidx = 2.0 / base
  detady = 2.0 / height
  J_inv = np.diag([dqsidx, detady])

  if nstate == 2: # compute elasticity stiffness
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      dN_xy = J_inv.T @ dN
      B = createB(dN_xy)
      phase_field = sum(N[0, nstate * i] * Upf[element.node_ids[i]] for i in range(nnodes))
      Ddeteriorated = D.copy()
      Ddeteriorated *= (1 - phase_field) ** 2
      Ke += B.T @ Ddeteriorated @ B * qp.weight * detjac
  elif nstate == 1: # compute phase field stiffness
    Gc, l0 = mat.Gc, mat.l0
    c0 = 2.0
    for qp in intrule:
      N, dN = shapeFunctions(qp.xi, qp.eta, nstate)
      dN_xy = J_inv.T @ dN # Same as B_phi
      sigmaDotEps = calculateSigmaDotEps(element, dN_xy)
      Ke += detjac * qp.weight * (Gc * l0 / c0 * (dN_xy.T @ dN_xy) + (Gc / (l0 * c0) + 0.5 * sigmaDotEps) * N.T @ N)
      Fe += detjac * qp.weight * 0.5 * sigmaDotEps * N.flatten()
  else:
    raise Exception("Invalid nstate")

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
  return stress_vec

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

def applyBoundaryConditions(K, F, bc_nodes):
  for bc in bc_nodes:
    row = 2 * bc.node
    xval = bc.xval * pseudotime
    yval = bc.yval * pseudotime
    if bc.bctype == 0: # displacement in x and y
      F -= K[:, row] * xval
      F -= K[:, row + 1] * yval
      K[row, :] = 0
      K[:, row] = 0
      K[row + 1, :] = 0
      K[:, row + 1] = 0
      K[row, row] = 1.0
      K[row + 1, row + 1] = 1.0
      F[row] = xval
      F[row + 1] = yval
    elif bc.bctype == 1: # displacement in x
      F -= K[:, row] * xval
      K[row, :] = 0
      K[:, row] = 0
      K[row, row] = 1.0
      F[row] = xval
    elif bc.bctype == 2: # displacement in y
      F -= K[:, row + 1] * yval
      K[row + 1, :] = 0
      K[:, row + 1] = 0
      K[row + 1, row + 1] = 1.0
      F[row + 1] = yval
    elif bc.bctype == 3: # nodal load in x and y
      F[row] += xval
      F[row + 1] += yval

def solveSystem(K, F, U):
  time = Timer()
  U[:] = np.linalg.solve(K, F)
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
  E = 30
  nu = 0.2
  Gc = 1.2e-4
  l0 = 10.0

  num_elements_x = 100
  num_elements_y = 10
  length = 200.0
  height = 20.0
  dt = 0.02
  totaltime = 1.5
  maxsteps = int(1e5)
  maxiter = 600
  stagtol = 1e-6

  sigma_peak_at2 = np.sqrt(27.0 * E * Gc / (256.0 * l0))
  u_peak_at2 = 16.0 / 9.0 * sigma_peak_at2 * length / E
  print(f"Sigma peak: {sigma_peak_at2}")
  print(f"U peak: {u_peak_at2}")
  imposed_displacement_x = u_peak_at2

  nodes = np.array([], dtype=object)
  elements = np.array([], dtype=object)
  bc_nodes = np.array([], dtype=object)

  material = MaterialParameters(E, nu, Gc, l0)
  factor = E / (1 - nu * nu)
  global D
  D = np.zeros((3, 3))
  D[0, 0] = factor
  D[0, 1] = factor * nu
  D[1, 0] = factor * nu
  D[1, 1] = factor
  D[2, 2] = factor * (1 - nu) / 2.0

  createRectangularMesh(nodes, elements, num_elements_x, num_elements_y, length, height)
  
  for i in range(len(nodes)):
    if abs(nodes[i].x) < 1e-8:
      bc_nodes = np.append(bc_nodes, BC(i, 0, 0.0, 0.0))
    elif abs(nodes[i].x - length) < 1e-8:
      bc_nodes = np.append(bc_nodes, BC(i, 1, imposed_displacement_x, 0.0))

  nstate_elas = 2
  nstate_pf = 1
  nnodes = len(nodes)
  ndofs_elas = nstate_elas * nnodes
  ndofs_pf = nstate_pf * nnodes
  Kelas = np.zeros((ndofs_elas, ndofs_elas))
  Felas = np.zeros(ndofs_elas)
  global Uelas
  Uelas = np.zeros(ndofs_elas)

  Kpf = np.zeros((ndofs_pf, ndofs_pf))
  Fpf = np.zeros(ndofs_pf)  
  global Upf
  Upf = np.zeros(ndofs_pf)

  # Data structure to save the data
  time_data = []
  stress_data = []

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
        residual = Kelas @ Uelas - Felas
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

    sig = np.zeros(3)
    computeSigmaAtCenter(elements[50], nodes, sig) # element 50 is in the middle of the domain

    # Save the data
    time_data.append(pseudotime)
    stress_data.append(sig[0]/sigma_peak_at2)

  # Plot the data using matplotlib
  plt.figure()
  plt.plot(time_data, stress_data, 'o')
  plt.xlabel('Pseudo Time')
  plt.ylabel('Stress/Stress_peak')
  plt.title('Stress vs Pseudo Time')
  plt.grid(True)
  plt.savefig('outputs/stress_vs_time.png')
  # plt.show()

  print("\n================> Simulation completed!")
  simulation_time.elapsed("complete simulation")

if __name__ == "__main__":
  main()
