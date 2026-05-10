function A = assembleStiffness2D(node, elem)
% ASSEMBLESTIFFNESS2D  Assemble the P1 stiffness matrix on a 2D triangular mesh.
%
%   A_ij = \int_\Omega \nabla \phi_i \cdot \nabla \phi_j  dx
%
%   A = ASSEMBLESTIFFNESS2D(node, elem)
%
%   Input:
%     node - N x 2  vertex coordinates
%     elem - NT x 3 element connectivity (1-indexed)
%   Output:
%     A    - N x N sparse stiffness matrix
%
%   Fully vectorised over elements — no per-element loop.
%   For P1 elements, basis gradients are constant on each triangle, so no
%   quadrature loop is needed.

N = size(node, 1);
NT = size(elem, 1);

% Element vertices:  [x1,y1], [x2,y2], [x3,y3]
x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);

% Signed area * 2:  2|T| = (x2-x1)(y3-y1) - (x3-x1)(y2-y1)
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);

% Gradient of barycentric coordinates (constant per element):
%   \nabla \lambda_1 = [ y2 - y3;  x3 - x2 ] / (2|T|)
%   \nabla \lambda_2 = [ y3 - y1;  x1 - x3 ] / (2|T|)
%   \nabla \lambda_3 = [ y1 - y2;  x2 - x1 ] / (2|T|)
g1x = (y2 - y3) ./ area2;    g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2;    g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2;    g3y = (x2 - x1) ./ area2;

% Local stiffness:  K_loc(i,j) = |T| * (\nabla\lambda_i \cdot \nabla\lambda_j)
% where |T| = |area2| / 2
area = abs(area2) / 2;

k11 = area .* (g1x.^2 + g1y.^2);
k22 = area .* (g2x.^2 + g2y.^2);
k33 = area .* (g3x.^2 + g3y.^2);
k12 = area .* (g1x .* g2x + g1y .* g2y);
k13 = area .* (g1x .* g3x + g1y .* g3y);
k23 = area .* (g2x .* g3x + g2y .* g3y);

% ---- Sparse assembly using (i,j,value) triplets ---------------------------
% 3x3 symmetric local matrix → 9 entries per element
% Diagonal (3 entries per element)
ii = [elem(:,1);  elem(:,2);  elem(:,3)];
jj = [elem(:,1);  elem(:,2);  elem(:,3)];
ss = [k11;        k22;        k33];

% Off-diagonal (6 entries per element, symmetric)
ii = [ii;  elem(:,1);  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,2);  elem(:,3)];
jj = [jj;  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,1);  elem(:,3);  elem(:,2)];
ss = [ss;  k12;        k12;        k13;        k13;        k23;        k23];

A = sparse(ii, jj, ss, N, N);
end
