function P = prolongateNestedP1(coarseNode, coarseElem, fineNode)
% PROLONGATENESTEDP1  P1 prolongation matrix from a coarse mesh to fine nodes.
%
%   P maps coarse P1 nodal values to their values at fineNode.
%   The meshes must be nested geometrically.

Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);
dim = size(coarseNode, 2);
nv = dim + 1;

[elemId, lambda] = locateSimplexP1(coarseNode, coarseElem, fineNode, 1e-10);
if any(elemId == 0)
    bad = find(elemId == 0, 1);
    error('prolongateNestedP1:notNested', ...
        'Fine node %d was not found in the coarse mesh.', bad);
end

ii = zeros(nv * Nf, 1);
jj = zeros(nv * Nf, 1);
ss = zeros(nv * Nf, 1);
idx = 0;
for q = 1:Nf
    cverts = coarseElem(elemId(q), :);
    rows = idx + (1:nv);
    ii(rows) = q;
    jj(rows) = cverts(:);
    ss(rows) = lambda(q, :).';
    idx = idx + nv;
end

P = sparse(ii, jj, ss, Nf, Nc);
end
