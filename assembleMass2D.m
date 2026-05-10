function M = assembleMass2D(node, elem)
% ASSEMBLEMASS2D  Assemble the P1 mass matrix on a 2D triangular mesh.
%
%   M_ij = \int_\Omega \phi_i \phi_j  dx
%
%   M = ASSEMBLEMASS2D(node, elem)
%
%   Uses the closed-form local mass matrix for P1 triangles:
%     M_loc = |T|/12 * [2 1 1; 1 2 1; 1 1 2]
%   No quadrature loop needed.

N = size(node, 1);

x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);

area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));

% Local mass:  diag = 2*|T|/12 = |T|/6,   off-diag = |T|/12
diag_val = area / 6;
off_val  = area / 12;

% ---- Sparse assembly: 9 entries per 3x3 symmetric matrix ------------------
ii = [elem(:,1);  elem(:,2);  elem(:,3)];
jj = [elem(:,1);  elem(:,2);  elem(:,3)];
ss = [diag_val;   diag_val;   diag_val];

ii = [ii;  elem(:,1);  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,2);  elem(:,3)];
jj = [jj;  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,1);  elem(:,3);  elem(:,2)];
ss = [ss;  off_val;    off_val;    off_val;    off_val;    off_val;    off_val];

M = sparse(ii, jj, ss, N, N);
end
