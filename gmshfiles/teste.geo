// =================================================================
// 1. PARÂMETROS DA GEOMETRIA (mm)
// =================================================================

L = 8.0;      // Comprimento da viga
H = 2.0;      // Altura da viga
a = 0.4;      // Comprimento da trinca inicial
c = 0.001;    // Metade da espessura da trinca

// Parâmetros da malha
h_grosso = 0.30;
h_fino   = 0.03;

// Dimensões do triângulo refinado
x_triangulo = -0.75;
y_triangulo =  1.5;


// =================================================================
// 2. PONTOS DA GEOMETRIA PRINCIPAL
// =================================================================

// Base esquerda
Point(1) = {-L/2, 0, 0};
Point(2) = {-c,   0, 0};

// Trinca
Point(3) = {-c, a, 0};  // Ponta esquerda da trinca
Point(4) = { c, a, 0};  // Ponta direita da trinca
Point(5) = { c, 0, 0};

// Base direita
Point(6) = {L/2, 0, 0};

// Parte superior
Point(7) = { L/2, H, 0};
Point(8) = { 0,   H, 0};  // Ponto de aplicação da carga
Point(9) = {-L/2, H, 0};


// =================================================================
// 3. PONTOS DO TRIÂNGULO REFINADO
// =================================================================
//
// O terceiro vértice é o Point(3), isto é, a ponta da trinca.
//
//        Point(10) ---------------- Point(11)
//               \                    |
//                \                   |
//                 \                  |
//                  \                |
//                   \               |
//                    \               |
//                     Point(3) ------
//
// Ajuste x_triangulo e y_triangulo para alterar seu tamanho.
//

Point(10) = {x_triangulo, y_triangulo, 0};
Point(11) = {-c,          y_triangulo, 0};


// =================================================================
// 4. LINHAS DO CONTORNO PRINCIPAL
// =================================================================

Line(1) = {1, 2};  // Base esquerda

Line(2) = {2, 3};  // Lábio esquerdo da trinca
Line(3) = {3, 4};  // Topo da trinca
Line(4) = {4, 5};  // Lábio direito da trinca

Line(5) = {5, 6};  // Base direita
Line(6) = {6, 7};  // Borda direita
Line(7) = {7, 8};  // Topo direito
Line(8) = {8, 9};  // Topo esquerdo
Line(9) = {9, 1};  // Borda esquerda


// =================================================================
// 5. LINHAS DO TRIÂNGULO REFINADO
// =================================================================

Line(10) = {10, 11}; // Lado horizontal superior
Line(11) = {11, 3};  // Lado vertical direito
Line(12) = {3, 10};  // Lado inclinado


// =================================================================
// 6. SUPERFÍCIES
// =================================================================

// Contorno externo completo
Curve Loop(1) = {
    1, 2, 3, 4, 5, 6, 7, 8, 9
};

// Contorno triangular
Curve Loop(2) = {
    10, 11, 12
};

// Região externa ao triângulo
Plane Surface(1) = {1, 2};

// Região interna do triângulo
Plane Surface(2) = {2};


// =================================================================
// 7. GRUPOS FÍSICOS
// =================================================================

Physical Point("Left_ids")  = {1};
Physical Point("Right_ids") = {6};
Physical Point("Top_ids")   = {8};

// As duas superfícies pertencem ao mesmo domínio físico
Physical Surface("Dominio") = {1, 2};


// =================================================================
// 8. CAMPO DE REFINAMENTO
// =================================================================
//
// Surface(2): interior do triângulo
// Surface(1): restante da viga
//

Field[1] = Constant;

// Superfície que receberá h_fino
Field[1].SurfacesList = {2};

// Tamanho dentro do triângulo
Field[1].VIn = h_fino;

// Tamanho no restante da viga
Field[1].VOut = h_grosso;

// Inclui as linhas que delimitam o triângulo
Field[1].IncludeBoundary = 1;

Background Field = 1;


// =================================================================
// 9. CONTROLE DAS FONTES DE TAMANHO
// =================================================================

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;


// =================================================================
// 10. MALHA 2D EM QUADRILÁTEROS
// =================================================================

// Recombinam-se as duas superfícies
Recombine Surface {1, 2};

// Blossom
Mesh.RecombinationAlgorithm = 1;

// Recombinação global
Mesh.RecombineAll = 1;

// Divide eventuais triângulos remanescentes em quadriláteros
Mesh.SubdivisionAlgorithm = 1;