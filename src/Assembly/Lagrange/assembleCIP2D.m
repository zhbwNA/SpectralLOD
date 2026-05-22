function C = assembleCIP2D(node, elem, degree, gamma, opts)
% ASSEMBLECIP2D  Continuous interior penalty matrix for 2D Lagrange FEM.
%
%   C = ASSEMBLECIP2D(node,elem,degree,gamma)
%
%   C contains the interior-edge terms
%       sum_e sum_{j=1}^p i*gamma_j*h_e^(2j-1)
%           int_e [d_n^j u] [d_n^j v] ds.
%
%   degree is 1, 2, or 3. gamma may be empty, scalar, vector, or a function
%   handle gamma(j, edgeMidpoint, hEdge). By default real gamma values are
%   multiplied by 1i. Set opts.gammaIsCoefficient = true to pass the final
%   complex coefficients directly.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4, gamma = []; end
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'gammaIsCoefficient') || isempty(opts.gammaIsCoefficient)
    opts.gammaIsCoefficient = false;
end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2, degree + 1);
end

if degree < 1 || degree > 3
    error('assembleCIP2D:degree', 'degree must be 1, 2, or 3.');
end

baseElem = elem(:, 1:3);
if degree == 1
    nodeH = node;
    elemH = baseElem;
else
    if size(elem, 2) == 3
        [nodeH, elemH] = extendMesh2D(node, baseElem, degree);
    else
        nodeH = node;
        elemH = elem;
    end
end

N = size(nodeH, 1);
nLB = size(elemH, 2);
[edge, edgeIdx, ~, edge2elem] = edgeMesh2D(baseElem);
interior = find(edge2elem(:,2) > 0);

[tq, wq] = gauss1D(min(6, opts.quadOrder));
nq = length(wq);

maxNnz = numel(interior) * degree * nq * (2*nLB)^2;
ii = zeros(maxNnz, 1);
jj = zeros(maxNnz, 1);
ss = zeros(maxNnz, 1);
idx = 0;

monos = lagrangeMonomials2D(degree);
localEdgePairs = [1 2; 2 3; 3 1];

for ee = reshape(interior, 1, [])
    tPlus = edge2elem(ee, 1);
    tMinus = edge2elem(ee, 2);
    locPlus = find(edgeIdx(tPlus, :) == ee, 1);
    locMinus = find(edgeIdx(tMinus, :) == ee, 1);

    xA = nodeH(edge(ee,1), :);
    xB = nodeH(edge(ee,2), :);
    tau = xB - xA;
    hEdge = norm(tau);
    normal = [tau(2), -tau(1)] / hEdge;
    midpoint = 0.5 * (xA + xB);

    gradPlus = baryGradients2D(nodeH(baseElem(tPlus, :), :));
    gradMinus = baryGradients2D(nodeH(baseElem(tMinus, :), :));
    alphaPlus = gradPlus * normal(:);
    alphaMinus = gradMinus * normal(:);

    dofs = [elemH(tPlus, :), elemH(tMinus, :)];
    for q = 1:nq
        lamPlus = edgeLambda(baseElem(tPlus, :), localEdgePairs(locPlus, :), edge(ee, :), tq(q));
        lamMinus = edgeLambda(baseElem(tMinus, :), localEdgePairs(locMinus, :), edge(ee, :), tq(q));
        weightEdge = hEdge * wq(q);

        for j = 1:degree
            coeff = cipCoefficient(gamma, j, midpoint, hEdge, opts) * hEdge^(2*j - 1);
            dPlus = directionalDerivativeMonomials(monos, lamPlus, alphaPlus(:).', j);
            dMinus = directionalDerivativeMonomials(monos, lamMinus, alphaMinus(:).', j);
            jump = [dPlus, -dMinus];

            rows = idx + (1:numel(dofs)^2);
            [aa, bb] = ndgrid(1:numel(dofs), 1:numel(dofs));
            ii(rows) = dofs(aa(:));
            jj(rows) = dofs(bb(:));
            ss(rows) = coeff * weightEdge * (jump(aa(:)) .* jump(bb(:)));
            idx = rows(end);
        end
    end
end

C = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function coeff = cipCoefficient(gamma, j, midpoint, hEdge, opts)
if isempty(gamma)
    defaults = [0.1, 0.01, 0.001];
    g = defaults(j);
elseif isa(gamma, 'function_handle')
    g = gamma(j, midpoint, hEdge);
elseif isscalar(gamma)
    g = gamma;
else
    g = gamma(j);
end

if opts.gammaIsCoefficient
    coeff = g;
else
    coeff = 1i * g;
end
end


function lam = edgeLambda(elemVerts, pair, globalEdge, t)
lam = zeros(1, 3);
if elemVerts(pair(1)) == globalEdge(1)
    lam(pair(1)) = 1 - t;
    lam(pair(2)) = t;
else
    lam(pair(1)) = t;
    lam(pair(2)) = 1 - t;
end
end


function G = baryGradients2D(v)
x1 = v(1,1); y1 = v(1,2);
x2 = v(2,1); y2 = v(2,2);
x3 = v(3,1); y3 = v(3,2);
area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
G = [(y2 - y3) / area2, (x3 - x2) / area2;
     (y3 - y1) / area2, (x1 - x3) / area2;
     (y1 - y2) / area2, (x2 - x1) / area2];
end


function monos = lagrangeMonomials2D(degree)
switch degree
    case 1
        monos = {
            [1 0 0 1]
            [0 1 0 1]
            [0 0 1 1]
        };
    case 2
        monos = {
            [2 0 0 2; 1 0 0 -1]
            [0 2 0 2; 0 1 0 -1]
            [0 0 2 2; 0 0 1 -1]
            [1 1 0 4]
            [0 1 1 4]
            [1 0 1 4]
        };
    case 3
        monos = {
            [3 0 0 4.5; 2 0 0 -4.5; 1 0 0 1]
            [0 3 0 4.5; 0 2 0 -4.5; 0 1 0 1]
            [0 0 3 4.5; 0 0 2 -4.5; 0 0 1 1]
            [2 1 0 13.5; 1 1 0 -4.5]
            [1 2 0 13.5; 1 1 0 -4.5]
            [0 2 1 13.5; 0 1 1 -4.5]
            [0 1 2 13.5; 0 1 1 -4.5]
            [1 0 2 13.5; 1 0 1 -4.5]
            [2 0 1 13.5; 1 0 1 -4.5]
            [1 1 1 27]
        };
    otherwise
        error('lagrangeMonomials2D:degree', 'degree must be 1, 2, or 3.');
end
end


function vals = directionalDerivativeMonomials(monos, lambda, alpha, order)
nLB = numel(monos);
vals = zeros(1, nLB);
betas = derivativeMultiIndices(order);
multiCoef = multinomialRows(betas);
for a = 1:nLB
    terms = monos{a};
    val = 0;
    for r = 1:size(terms, 1)
        pow = terms(r, 1:3);
        c = terms(r, 4);
        for b = 1:size(betas, 1)
            beta = betas(b, :);
            if any(beta > pow)
                continue;
            end
            fall = prod(fallingFactorial(pow, beta));
            remPow = pow - beta;
            val = val + c * multiCoef(b) * fall * prod(alpha.^beta) * prod(lambda.^remPow);
        end
    end
    vals(a) = val;
end
end


function betas = derivativeMultiIndices(order)
switch order
    case 1
        betas = [1 0 0; 0 1 0; 0 0 1];
    case 2
        betas = [2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; 0 1 1];
    case 3
        betas = [3 0 0; 0 3 0; 0 0 3; 2 1 0; 2 0 1; ...
                 1 2 0; 0 2 1; 1 0 2; 0 1 2; 1 1 1];
    otherwise
        error('derivativeMultiIndices:order', 'order must be 1, 2, or 3.');
end
end


function c = multinomialRows(betas)
ord = sum(betas, 2);
c = factorial(ord) ./ prod(factorial(betas), 2);
end


function f = fallingFactorial(pow, beta)
f = ones(size(pow));
for k = 1:numel(pow)
    if beta(k) > 0
        f(k) = prod(pow(k) - (0:beta(k)-1));
    end
end
end
