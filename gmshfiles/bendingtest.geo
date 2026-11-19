// =================================================================
// 1. PARÂMETROS DA GEOMETRIA (Em mm)
// =================================================================
L = 8.0;      // Comprimento da viga
H = 2.0;      // Altura da viga
a = 0.4;      // Comprimento da trinca inicial
c = 0.001;    // Espessura da trinca (muito fina para emular fenda)

// Parâmetros de Malha

h_grosso = 0.3;   
h_fino   = 0.02;   

// =================================================================
// 2. PONTOS GEOMÉTRICOS
// =================================================================
Point(1) = {-L/2, 0, 0};
Point(2) = {-c,   0, 0};

// A Trinca (Sobe pelo lábio esquerdo, desce pelo lábio direito)
Point(3) = {-c,   a, 0}; // Ponta da trinca (esquerda)
Point(4) = { c,   a, 0}; // Ponta da trinca (direita)
Point(5) = { c,   0, 0}; // Base da trinca (direita)

// Base inferior (Trinca até à Direita)
Point(6) = { L/2, 0, 0};

// Topo
Point(7) = { L/2, H, 0};
Point(8) = { 0,   H, 0};   // Ponto onde a carga será aplicada
Point(9) = {-L/2, H, 0};

// =================================================================
// 3. LINHAS E CONTORNO
// =================================================================
Line(1) = {1, 2}; // Base esquerda
Line(2) = {2, 3}; // Lábio esquerdo da trinca
Line(3) = {3, 4}; // Topo da trinca
Line(4) = {4, 5}; // Lábio direito da trinca
Line(5) = {5, 6}; // Base direita
Line(6) = {6, 7}; // Borda direita
Line(7) = {7, 8}; // Topo direito
Line(8) = {8, 9}; // Topo esquerdo
Line(9) = {9, 1}; // Borda esquerda

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Plane Surface(1) = {1};

// =================================================================
// 4. GRUPOS FÍSICOS (Para identificar no Python)
// =================================================================
Physical Point("Left_ids") = {1};
Physical Point("Right_ids")  = {6};
Physical Point("Top_ids")    = {8};
Physical Surface("Dominio")      = {1};

// =================================================================
// 5. CAIXA DE REFINAMENTO (Mesh Field)
// =================================================================
Field[1] = Box;
Field[1].VIn = h_fino;
Field[1].VOut = h_grosso;
Field[1].XMin = -0.14;
Field[1].XMax = 0.14;
Field[1].YMin = 0.3;  
Field[1].YMax = H - 0.23;

Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// =================================================================
// 6. MALHA 2D APENAS COM QUADRILÁTEROS (CORRIGIDO)
// =================================================================
// Usa o algoritmo Frontal 2D (altamente robusto e estável)
// Recombina triângulos em quadriláteros
Recombine Surface {1};

// Algoritmo de recombinação (Blossom é o melhor para quads)
Mesh.RecombinationAlgorithm = 1; 

// Força a malha a ser estritamente de quadriláteros
Mesh.RecombineAll = 1;

// SE houver triângulos remanescentes, este comando os divide em quads:
Mesh.SubdivisionAlgorithm = 1;