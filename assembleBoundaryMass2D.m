function Mb = assembleBoundaryMass2D(node, elem, bdFlag)
% ASSEMBLEBOUNDARYMASS2D  Assemble the P1 boundary mass matrix on a 2D mesh.
%
%   Mb_ij = \int_{\partial\Omega} \phi_i \phi_j  ds
%
%   Mb = ASSEMBLEBOUNDARYMASS2D(node, elem, bdFlag)
%
%   bdFlag(t,k) = 1 if edge k of element t lies on the Dirichlet boundary.
%   (Edge k is the edge opposite vertex k.)
%
%   Uses the closed-form 1D boundary mass on edges:
%     Mb_loc = L/6 * [2 1; 1 2]    where L = edge length

N = size(node, 1);

% Pre-allocate for all possible boundary edges (3 edges * NT)
maxBd = 3 * size(elem, 1);
ii = zeros(maxBd * 4, 1);                % 4 entries per edge (2x2 symmetric)
jj = zeros(maxBd * 4, 1);
ss = zeros(maxBd * 4, 1);
idx = 0;

% For each of the 3 edges: extract those flagged as boundary and assemble
edgeVertex = [2 3;   % edge 1 (opposite v1): vertices v2, v3
              3 1;   % edge 2 (opposite v2): vertices v3, v1
              1 2];  % edge 3 (opposite v3): vertices v1, v2

for k = 1:3                               % 3 edges, small loop — not over elements
    bdEdges = (bdFlag(:,k) == 1);
    if ~any(bdEdges), continue; end

    e = elem(bdEdges, :);                 % elements contributing to this boundary edge
    vA = e(:, edgeVertex(k,1));
    vB = e(:, edgeVertex(k,2));

    % Edge length  L = |vB - vA|
    L = sqrt((node(vB,1) - node(vA,1)).^2 + ...
             (node(vB,2) - node(vA,2)).^2);

    % Local boundary mass (closed form):  L/6 * [2 1; 1 2]
    nBd = length(L);
    nxt = idx + 1;
    idx = idx + 4 * nBd;

    % Diagonal: L*2/6 = L/3
    ii(nxt:4:idx)   = vA;  jj(nxt:4:idx)   = vA;  ss(nxt:4:idx)   = L / 3;
    ii(nxt+1:4:idx) = vB;  jj(nxt+1:4:idx) = vB;  ss(nxt+1:4:idx) = L / 3;
    % Off-diagonal: L/6
    ii(nxt+2:4:idx) = vA;  jj(nxt+2:4:idx) = vB;  ss(nxt+2:4:idx) = L / 6;
    ii(nxt+3:4:idx) = vB;  jj(nxt+3:4:idx) = vA;  ss(nxt+3:4:idx) = L / 6;
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
