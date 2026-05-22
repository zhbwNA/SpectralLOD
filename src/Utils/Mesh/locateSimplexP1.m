function [elemId, lambda] = locateSimplexP1(node, elem, points, tol)
% LOCATESIMPLEXP1  Locate points in a P1 triangle/tetrahedron mesh.
%
%   [elemId,lambda] = LOCATESIMPLEXP1(node,elem,points)
%
%   elemId(q) is the containing element index, or 0 if not found.
%   lambda(q,:) are barycentric coordinates in that element.

if nargin < 4 || isempty(tol)
    tol = 1e-11;
end

dim = size(node, 2);
np = size(points, 1);
nv = dim + 1;
elemId = zeros(np, 1);
lambda = zeros(np, nv);

emin = zeros(size(elem, 1), dim);
emax = zeros(size(elem, 1), dim);
for d = 1:dim
    vals = reshape(node(elem, d), size(elem));
    emin(:, d) = min(vals, [], 2) - tol;
    emax(:, d) = max(vals, [], 2) + tol;
end

for p = 1:np
    x = points(p, :);
    cand = true(size(elem, 1), 1);
    for d = 1:dim
        cand = cand & x(d) >= emin(:, d) & x(d) <= emax(:, d);
    end
    ids = find(cand).';
    for t = ids
        lam = barycentricPoint(node(elem(t, :), :), x);
        if all(lam >= -tol) && all(lam <= 1 + tol)
            lam(abs(lam) < tol) = 0;
            lam(abs(lam - 1) < tol) = 1;
            elemId(p) = t;
            lambda(p, :) = lam / sum(lam);
            break;
        end
    end
end
end


function lam = barycentricPoint(v, x)
dim = size(v, 2);
if dim == 2
    B = [v(1,:) - v(3,:); v(2,:) - v(3,:)].';
    rhs = (x - v(3,:)).';
    a = B \ rhs;
    lam = [a(1), a(2), 1 - a(1) - a(2)];
elseif dim == 3
    B = [v(1,:) - v(4,:); v(2,:) - v(4,:); v(3,:) - v(4,:)].';
    rhs = (x - v(4,:)).';
    a = B \ rhs;
    lam = [a(1), a(2), a(3), 1 - sum(a)];
else
    error('locateSimplexP1:dim', 'Only 2D and 3D meshes are supported.');
end
end
