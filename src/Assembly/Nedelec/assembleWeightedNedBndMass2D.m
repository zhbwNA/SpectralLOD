function Mb = assembleWeightedNedBndMass2D(node, elem, bdFlag, coef)
% ASSEMBLEWEIGHTEDNEDBNDMASS2D  Assemble int_Gamma coef (w_i.t)(w_j.t) ds.

if nargin < 4 || isempty(coef), coef = 1; end
if isnumeric(coef) && isscalar(coef)
    Mb = coef * assembleNedBndMass2D(node, elem, bdFlag);
    return;
end

[~, edgeIdx] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
bdFlagToEdgeIdx = [2, 3, 1];
edgeVerts = [2 3; 3 1; 1 2];
[xi, w] = localGauss1D01(3);
rows = [];
vals = [];

for e = 1:3
    isBd = bdFlag(:, e) == 1;
    if ~any(isBd), continue; end
    tri = elem(isBd, :);
    vA = tri(:, edgeVerts(e,1));
    vB = tri(:, edgeVerts(e,2));
    xA = node(vA, 1); yA = node(vA, 2);
    xB = node(vB, 1); yB = node(vB, 2);
    L = sqrt((xB - xA).^2 + (yB - yA).^2);
    intCoef = zeros(size(L));
    for q = 1:numel(w)
        xq = (1 - xi(q)) * xA + xi(q) * xB;
        yq = (1 - xi(q)) * yA + xi(q) * yB;
        intCoef = intCoef + w(q) * L .* evalPDECoefficient(coef, xq, yq, [], []);
    end
    eid = edgeIdx(isBd, bdFlagToEdgeIdx(e));
    rows = [rows; eid(:)]; %#ok<AGROW>
    vals = [vals; intCoef(:) ./ (L(:).^2)]; %#ok<AGROW>
end

Mb = sparse(rows, rows, vals, NE, NE);
end


function [x, w] = localGauss1D01(n)
[x, w] = gauss1D(n);
end
