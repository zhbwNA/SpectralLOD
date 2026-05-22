function Mb = assembleWeightedBoundaryMass2D(node, elem, bdFlag, degree, coef)
% ASSEMBLEWEIGHTEDBOUNDARYMASS2D  Assemble int_Gamma coef phi_i phi_j ds.

if nargin < 4 || isempty(degree), degree = 1; end
if nargin < 5 || isempty(coef), coef = 1; end

if isnumeric(coef) && isscalar(coef)
    Mb = coef * assembleBoundaryMass2D(node, elem, bdFlag, degree);
    return;
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
[xi, w1d] = localGauss1D01(max(2 * degree, 2));
phi1d = localLagrange1D(degree, xi);
nEdgeDof = degree + 1;
edgeVerts = [2 3; 3 1; 1 2];
[aa, bb] = ndgrid(1:nEdgeDof, 1:nEdgeDof);
aa = aa(:).'; bb = bb(:).';

nEntries = nnz(bdFlag) * nEdgeDof^2;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for e = 1:3
    bdEdges = bdFlag(:, e) == 1;
    if ~any(bdEdges), continue; end
    eBd = elem(bdEdges, :);
    vA = eBd(:, edgeVerts(e,1));
    vB = eBd(:, edgeVerts(e,2));
    xA = node(vA, 1); yA = node(vA, 2);
    xB = node(vB, 1); yB = node(vB, 2);
    L = sqrt((xB - xA).^2 + (yB - yA).^2);
    edgeDofs = eBd(:, localEdgeDofs2D(e, degree));
    S = zeros(size(edgeDofs, 1), nEdgeDof^2);
    for q = 1:numel(w1d)
        xq = (1 - xi(q)) * xA + xi(q) * xB;
        yq = (1 - xi(q)) * yA + xi(q) * yB;
        cq = evalPDECoefficient(coef, xq, yq, [], []);
        S = S + (w1d(q) * L) .* cq .* phi1d(q, aa) .* phi1d(q, bb);
    end
    nNew = numel(S);
    nxt = idx + 1; idx = idx + nNew;
    ii(nxt:idx) = reshape(edgeDofs(:, aa), [], 1);
    jj(nxt:idx) = reshape(edgeDofs(:, bb), [], 1);
    ss(nxt:idx) = reshape(S, [], 1);
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function dofs = localEdgeDofs2D(edgeK, degree)
switch degree
    case 1
        map = {[2, 3], [3, 1], [1, 2]};
    case 2
        map = {[2, 5, 3], [3, 6, 1], [1, 4, 2]};
    case 3
        map = {[2, 6, 7, 3], [3, 8, 9, 1], [1, 4, 5, 2]};
    otherwise
        error('localEdgeDofs2D:degree', 'Degree %d is not supported.', degree);
end
dofs = map{edgeK};
end


function [x, w] = localGauss1D01(n)
[x, w] = gauss1D(n);
end


function phi = localLagrange1D(degree, x)
switch degree
    case 1
        phi = [1 - x, x];
    case 2
        phi = [(2*x - 1).*(x - 1), 4*x.*(1 - x), x.*(2*x - 1)];
    case 3
        phi = zeros(numel(x), 4);
        phi(:,1) = -4.5 * (x - 1/3) .* (x - 2/3) .* (x - 1);
        phi(:,2) = 13.5 * x .* (x - 2/3) .* (x - 1);
        phi(:,3) = -13.5 * x .* (x - 1/3) .* (x - 1);
        phi(:,4) = 4.5 * x .* (x - 1/3) .* (x - 2/3);
    otherwise
        error('localLagrange1D:degree', 'Degree %d is not supported.', degree);
end
end
