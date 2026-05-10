function A = assembleStiffness3D(node, elem)
% ASSEMBLESTIFFNESS3D  Assemble the P1 stiffness matrix on a 3D tetrahedral mesh.
%
%   A_ij = \int_\Omega \nabla \phi_i \cdot \nabla \phi_j  dx
%
%   A = ASSEMBLESTIFFNESS3D(node, elem)
%
%   Input:
%     node - N x 3  vertex coordinates
%     elem - NT x 4 element connectivity (1-indexed)
%   Output:
%     A    - N x N sparse stiffness matrix
%
%   Gradient computation via Cramer's rule using vectorised cross products
%   — no per-element loops, no matrix inversion.

N = size(node, 1);
NT = size(elem, 1);

% Element vertices as 3D coordinates (NT x 3)
v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

% Edge vectors from vertex 1
e12 = v2 - v1;                            % NT x 3
e13 = v3 - v1;
e14 = v4 - v1;

% Jacobian determinant:  det(J) = e12 · (e13 × e14) = 6 * signed_volume
detJ = e12(:,1) .* (e13(:,2).*e14(:,3) - e13(:,3).*e14(:,2)) ...
     + e12(:,2) .* (e13(:,3).*e14(:,1) - e13(:,1).*e14(:,3)) ...
     + e12(:,3) .* (e13(:,1).*e14(:,2) - e13(:,2).*e14(:,1));

volume = abs(detJ) / 6;                   % |T|

% Gradient of barycentric coordinates via Cramer's rule:
%   J = [e12, e13, e14]   (3x3, columns)
%   J^{-T} rows:  cross(e13,e14)/detJ,  cross(e14,e12)/detJ,  cross(e12,e13)/detJ
%
%   \nabla \lambda_2 = cross(e13, e14) / detJ
%   \nabla \lambda_3 = cross(e14, e12) / detJ
%   \nabla \lambda_4 = cross(e12, e13) / detJ
%   \nabla \lambda_1 = -(\nabla\lambda_2 + \nabla\lambda_3 + \nabla\lambda_4)

c2 = cross(e13, e14);                     % = detJ * \nabla\lambda_2,  NT x 3
c3 = cross(e14, e12);                     % = detJ * \nabla\lambda_3
c4 = cross(e12, e13);                     % = detJ * \nabla\lambda_4
c1 = -(c2 + c3 + c4);                     % = detJ * \nabla\lambda_1

invDetJ = 1 ./ detJ;

g1x = c1(:,1) .* invDetJ;  g1y = c1(:,2) .* invDetJ;  g1z = c1(:,3) .* invDetJ;
g2x = c2(:,1) .* invDetJ;  g2y = c2(:,2) .* invDetJ;  g2z = c2(:,3) .* invDetJ;
g3x = c3(:,1) .* invDetJ;  g3y = c3(:,2) .* invDetJ;  g3z = c3(:,3) .* invDetJ;
g4x = c4(:,1) .* invDetJ;  g4y = c4(:,2) .* invDetJ;  g4z = c4(:,3) .* invDetJ;

% Local stiffness:  K_loc(i,j) = |T| * (\nabla\lambda_i · \nabla\lambda_j)
k11 = volume .* (g1x.^2 + g1y.^2 + g1z.^2);
k22 = volume .* (g2x.^2 + g2y.^2 + g2z.^2);
k33 = volume .* (g3x.^2 + g3y.^2 + g3z.^2);
k44 = volume .* (g4x.^2 + g4y.^2 + g4z.^2);
k12 = volume .* (g1x.*g2x + g1y.*g2y + g1z.*g2z);
k13 = volume .* (g1x.*g3x + g1y.*g3y + g1z.*g3z);
k14 = volume .* (g1x.*g4x + g1y.*g4y + g1z.*g4z);
k23 = volume .* (g2x.*g3x + g2y.*g3y + g2z.*g3z);
k24 = volume .* (g2x.*g4x + g2y.*g4y + g2z.*g4z);
k34 = volume .* (g3x.*g4x + g3y.*g4y + g3z.*g4z);

% ---- Sparse assembly: 4x4 symmetric → 16 entries per element --------------
% Diagonal (4 entries)
ii = [elem(:,1);  elem(:,2);  elem(:,3);  elem(:,4)];
jj = [elem(:,1);  elem(:,2);  elem(:,3);  elem(:,4)];
ss = [k11;        k22;        k33;        k44];

% (1,2) & (2,1)
ii = [ii; elem(:,1); elem(:,2)];
jj = [jj; elem(:,2); elem(:,1)];
ss = [ss; k12; k12];
% (1,3) & (3,1)
ii = [ii; elem(:,1); elem(:,3)];
jj = [jj; elem(:,3); elem(:,1)];
ss = [ss; k13; k13];
% (1,4) & (4,1)
ii = [ii; elem(:,1); elem(:,4)];
jj = [jj; elem(:,4); elem(:,1)];
ss = [ss; k14; k14];
% (2,3) & (3,2)
ii = [ii; elem(:,2); elem(:,3)];
jj = [jj; elem(:,3); elem(:,2)];
ss = [ss; k23; k23];
% (2,4) & (4,2)
ii = [ii; elem(:,2); elem(:,4)];
jj = [jj; elem(:,4); elem(:,2)];
ss = [ss; k24; k24];
% (3,4) & (4,3)
ii = [ii; elem(:,3); elem(:,4)];
jj = [jj; elem(:,4); elem(:,3)];
ss = [ss; k34; k34];

A = sparse(ii, jj, ss, N, N);
end
