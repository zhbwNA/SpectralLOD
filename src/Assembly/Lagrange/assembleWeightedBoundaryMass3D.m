function Mb = assembleWeightedBoundaryMass3D(node, elem, bdFlag, degree, coef)
% ASSEMBLEWEIGHTEDBOUNDARYMASS3D  Assemble int_Gamma coef phi_i phi_j ds.

if nargin < 4 || isempty(degree), degree = 1; end
if nargin < 5 || isempty(coef), coef = 1; end

if isnumeric(coef) && isscalar(coef)
    Mb = coef * assembleBoundaryMass3D(node, elem, bdFlag, degree);
    return;
end

if degree > 1 && size(elem, 2) == 4
    [node, elem] = extendMesh3D(node, elem, degree);
end

N = size(node, 1);
nFaceDof = (degree + 1) * (degree + 2) / 2;
[lambda2, w2] = quadtriangle(max(2 * degree, 2));
[phiFace, ~] = lagrange2D(degree, lambda2);
faceVerts = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
[aa, bb] = ndgrid(1:nFaceDof, 1:nFaceDof);
aa = aa(:).'; bb = bb(:).';

nEntries = nnz(bdFlag) * nFaceDof^2;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for f = 1:4
    bdFaces = bdFlag(:, f) == 1;
    if ~any(bdFaces), continue; end
    tri = elem(bdFaces, :);
    fv = faceVerts{f};
    vA = tri(:, fv(1)); vB = tri(:, fv(2)); vC = tri(:, fv(3));
    xA = node(vA, :); xB = node(vB, :); xC = node(vC, :);
    cr = cross(xB - xA, xC - xA);
    area = 0.5 * sqrt(sum(cr.^2, 2));
    faceDofs = tri(:, localFaceDofs3D(f, degree));
    S = zeros(size(faceDofs, 1), nFaceDof^2);
    for q = 1:numel(w2)
        lq = lambda2(q, :);
        xq = lq(1) * xA(:,1) + lq(2) * xB(:,1) + lq(3) * xC(:,1);
        yq = lq(1) * xA(:,2) + lq(2) * xB(:,2) + lq(3) * xC(:,2);
        zq = lq(1) * xA(:,3) + lq(2) * xB(:,3) + lq(3) * xC(:,3);
        cq = evalPDECoefficient(coef, xq, yq, zq, []);
        S = S + (2 * w2(q) * area) .* cq .* phiFace(q, aa) .* phiFace(q, bb);
    end
    nNew = numel(S);
    nxt = idx + 1; idx = idx + nNew;
    ii(nxt:idx) = reshape(faceDofs(:, aa), [], 1);
    jj(nxt:idx) = reshape(faceDofs(:, bb), [], 1);
    ss(nxt:idx) = reshape(S, [], 1);
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function dofs = localFaceDofs3D(faceK, degree)
switch degree
    case 1
        map = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
    case 2
        map = {[2,3,4, 8,10,9], [1,3,4, 6,10,7], ...
               [1,2,4, 5,9,7], [1,2,3, 5,8,6]};
    case 3
        map = {[2,3,4, 11,12, 15,16, 14,13, 20], ...
               [1,3,4, 7,8, 15,16, 10,9, 19], ...
               [1,2,4, 5,6, 13,14, 9,10, 18], ...
               [1,2,3, 5,6, 11,12, 8,7, 17]};
    otherwise
        error('localFaceDofs3D:degree', 'Degree %d is not supported.', degree);
end
dofs = map{faceK};
end
