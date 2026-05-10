function M = assembleMass3D(node, elem)
% ASSEMBLEMASS3D  Assemble the P1 mass matrix on a 3D tetrahedral mesh.
%
%   M_ij = \int_\Omega \phi_i \phi_j  dx
%
%   M = ASSEMBLEMASS3D(node, elem)
%
%   Uses the closed-form local mass matrix for P1 tetrahedra:
%     M_loc = |T|/20 * [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2]

N = size(node, 1);

v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

e12 = v2 - v1;  e13 = v3 - v1;  e14 = v4 - v1;

% 6 * |T| = |det([e12, e13, e14])|
detJ = e12(:,1) .* (e13(:,2).*e14(:,3) - e13(:,3).*e14(:,2)) ...
     + e12(:,2) .* (e13(:,3).*e14(:,1) - e13(:,1).*e14(:,3)) ...
     + e12(:,3) .* (e13(:,1).*e14(:,2) - e13(:,2).*e14(:,1));

volume = abs(detJ) / 6;

% Local mass:  diag = 2*|T|/20 = |T|/10,   off-diag = |T|/20
diag_val = volume / 10;
off_val  = volume / 20;

% ---- Sparse assembly: 16 entries per 4x4 symmetric matrix -----------------
ii = [elem(:,1);  elem(:,2);  elem(:,3);  elem(:,4)];
jj = [elem(:,1);  elem(:,2);  elem(:,3);  elem(:,4)];
ss = [diag_val;   diag_val;   diag_val;   diag_val];

% (1,2) & (2,1)
ii = [ii; elem(:,1); elem(:,2)];
jj = [jj; elem(:,2); elem(:,1)];
ss = [ss; off_val; off_val];
% (1,3) & (3,1)
ii = [ii; elem(:,1); elem(:,3)];
jj = [jj; elem(:,3); elem(:,1)];
ss = [ss; off_val; off_val];
% (1,4) & (4,1)
ii = [ii; elem(:,1); elem(:,4)];
jj = [jj; elem(:,4); elem(:,1)];
ss = [ss; off_val; off_val];
% (2,3) & (3,2)
ii = [ii; elem(:,2); elem(:,3)];
jj = [jj; elem(:,3); elem(:,2)];
ss = [ss; off_val; off_val];
% (2,4) & (4,2)
ii = [ii; elem(:,2); elem(:,4)];
jj = [jj; elem(:,4); elem(:,2)];
ss = [ss; off_val; off_val];
% (3,4) & (4,3)
ii = [ii; elem(:,3); elem(:,4)];
jj = [jj; elem(:,4); elem(:,3)];
ss = [ss; off_val; off_val];

M = sparse(ii, jj, ss, N, N);
end
