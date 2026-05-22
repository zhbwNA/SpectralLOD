function Q = weightedClementP1(fineNode, fineElem, coarseNode, coarseElem, weightFun, quadOrder)
% WEIGHTEDCLEMENTP1  Matrix weighted Clement quasi-interpolation, P1 fine to P1 coarse.
%
%   Q = WEIGHTEDCLEMENTP1(fineNode,fineElem,coarseNode,coarseElem)
%   maps fine nodal values uf to coarse nodal values uc = Q*uf.
%
%   The row for a coarse node z is the weighted patch average
%       int w u_h Phi_z / int w Phi_z
%   assembled over the fine mesh, where Phi_z is the coarse P1 hat.

if nargin < 5 || isempty(weightFun)
    weightFun = [];
end
if nargin < 6 || isempty(quadOrder)
    quadOrder = 4;
end

dim = size(fineNode, 2);
nv = dim + 1;
Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);
NTf = size(fineElem, 1);

if dim == 2
    [lambdaF, wRef] = quadtriangle(quadOrder);
elseif dim == 3
    [lambdaF, wRef] = quadtet(quadOrder);
else
    error('weightedClementP1:dim', 'Only 2D and 3D meshes are supported.');
end
nQuad = length(wRef);

maxNnz = NTf * nQuad * nv * nv;
ii = zeros(maxNnz, 1);
jj = zeros(maxNnz, 1);
ss = zeros(maxNnz, 1);
denom = zeros(Nc, 1);
idx = 0;

for t = 1:NTf
    fv = fineElem(t, :);
    vFine = fineNode(fv, :);
    jacScale = simplexJacobianScale(vFine);
    xq = lambdaF * vFine;
    [ct, lambdaC] = locateSimplexP1(coarseNode, coarseElem, xq, 1e-10);
    if any(ct == 0)
        error('weightedClementP1:notNested', ...
            'A fine quadrature point was not found in the coarse mesh.');
    end

    for q = 1:nQuad
        wq = jacScale * wRef(q) * evalWeight(weightFun, xq(q, :));
        cv = coarseElem(ct(q), :);
        phiC = lambdaC(q, :);
        phiF = lambdaF(q, :);
        for a = 1:nv
            row = cv(a);
            denom(row) = denom(row) + wq * phiC(a);
            cols = fv(:);
            rows = idx + (1:nv);
            ii(rows) = row;
            jj(rows) = cols;
            ss(rows) = wq * phiC(a) * phiF(:);
            idx = idx + nv;
        end
    end
end

B = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Nc, Nf);
good = denom > 100 * eps(max(1, max(abs(denom))));
if any(~good)
    warning('weightedClementP1:emptyPatch', ...
        '%d coarse patch denominator(s) are zero.', nnz(~good));
    denom(~good) = 1;
end
Q = spdiags(1 ./ denom, 0, Nc, Nc) * B;
end


function s = simplexJacobianScale(v)
dim = size(v, 2);
if dim == 2
    area = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:)])) / 2;
    s = 2 * area;
else
    vol = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:); v(4,:) - v(1,:)])) / 6;
    s = 6 * vol;
end
end


function w = evalWeight(weightFun, x)
if isempty(weightFun)
    w = 1;
else
    if numel(x) == 2
        w = weightFun(x(1), x(2));
    else
        w = weightFun(x(1), x(2), x(3));
    end
end
end
