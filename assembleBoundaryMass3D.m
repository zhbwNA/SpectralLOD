function Mb = assembleBoundaryMass3D(node, elem, bdFlag)
% ASSEMBLEBOUNDARYMASS3D  Assemble the P1 boundary mass matrix on a 3D mesh.
%
%   Mb_ij = \int_{\partial\Omega} \phi_i \phi_j  ds
%
%   Mb = ASSEMBLEBOUNDARYMASS3D(node, elem, bdFlag)
%
%   bdFlag(t,f) = 1 if face f of tet t lies on the Dirichlet boundary.
%   (Face f is opposite vertex f.)
%
%   Uses the closed-form 2D boundary mass on triangular faces:
%     Mb_loc = area/12 * [2 1 1; 1 2 1; 1 1 2]

N = size(node, 1);
NT = size(elem, 1);

% Vertex lists for each of the 4 faces (cyclic ordering, opposite vertex f)
faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};

maxBd = 4 * NT;                          % upper bound
ii = zeros(maxBd * 9, 1);                % 9 entries per face (3x3 symmetric)
jj = zeros(maxBd * 9, 1);
ss = zeros(maxBd * 9, 1);
idx = 0;

for f = 1:4                               % 4 faces per tet (small loop)
    bdFaces = (bdFlag(:,f) == 1);
    if ~any(bdFaces), continue; end

    fv = faceVerts{f};
    e = elem(bdFaces, :);
    vA = e(:, fv(1));
    vB = e(:, fv(2));
    vC = e(:, fv(3));

    % Face area:  0.5 * |(vB - vA) × (vC - vA)|
    AB = node(vB, :) - node(vA, :);       % nBd x 3
    AC = node(vC, :) - node(vA, :);
    cr = cross(AB, AC);
    area = 0.5 * sqrt(cr(:,1).^2 + cr(:,2).^2 + cr(:,3).^2);

    % P1 local boundary mass (closed form):
    %   Mb_loc = area/12 * [2 1 1; 1 2 1; 1 1 2]
    nBd = length(area);
    nxt = idx + 1;
    idx = idx + 9 * nBd;

    % Diagonal: area * 2/12 = area/6
    ii(nxt:9:idx)   = vA;  jj(nxt:9:idx)   = vA;  ss(nxt:9:idx)   = area / 6;
    ii(nxt+1:9:idx) = vB;  jj(nxt+1:9:idx) = vB;  ss(nxt+1:9:idx) = area / 6;
    ii(nxt+2:9:idx) = vC;  jj(nxt+2:9:idx) = vC;  ss(nxt+2:9:idx) = area / 6;

    % Off-diagonal: area/12
    ii(nxt+3:9:idx) = vA;  jj(nxt+3:9:idx) = vB;  ss(nxt+3:9:idx) = area / 12;
    ii(nxt+4:9:idx) = vB;  jj(nxt+4:9:idx) = vA;  ss(nxt+4:9:idx) = area / 12;
    ii(nxt+5:9:idx) = vA;  jj(nxt+5:9:idx) = vC;  ss(nxt+5:9:idx) = area / 12;
    ii(nxt+6:9:idx) = vC;  jj(nxt+6:9:idx) = vA;  ss(nxt+6:9:idx) = area / 12;
    ii(nxt+7:9:idx) = vB;  jj(nxt+7:9:idx) = vC;  ss(nxt+7:9:idx) = area / 12;
    ii(nxt+8:9:idx) = vC;  jj(nxt+8:9:idx) = vB;  ss(nxt+8:9:idx) = area / 12;
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
