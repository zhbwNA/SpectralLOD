function jumpData = normalDerivativeJump2D(node, elem, degree, order, opts)
% NORMALDERIVATIVEJUMP2D  Evaluate Lagrange normal-derivative jumps on interior edges.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(order), order = 1; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2, degree + 1);
end
if degree < 1 || degree > 3 || order < 1 || order > degree
    error('normalDerivativeJump2D:order', 'degree and order must satisfy 1 <= order <= degree <= 3.');
end

baseElem = elem(:, 1:3);
if degree == 1
    nodeH = node;
    elemH = baseElem;
elseif size(elem, 2) == 3
    [nodeH, elemH] = extendMesh2D(node, baseElem, degree);
else
    nodeH = node;
    elemH = elem;
end

[edge, edgeIdx, ~, edge2elem] = edgeMesh2D(baseElem);
interior = find(edge2elem(:,2) > 0);
[tq, wq] = gauss1D(min(6, opts.quadOrder));
nQuad = numel(wq);
nLB = size(elemH, 2);
nTraceDof = 2 * nLB;
nEdge = numel(interior);

jumpData = struct();
jumpData.node = nodeH;
jumpData.elem = elemH;
jumpData.numDof = size(nodeH, 1);
jumpData.interiorEdge = interior;
jumpData.quadPoint = tq;
jumpData.quadWeight = wq(:).';
jumpData.order = order;
jumpData.degree = degree;
jumpData.dof = zeros(nEdge, nTraceDof);
jumpData.jump = zeros(nEdge, nQuad, nTraceDof);
jumpData.hEdge = zeros(nEdge, 1);
jumpData.midpoint = zeros(nEdge, 2);
if nEdge == 0
    return;
end

tPlus = edge2elem(interior, 1);
tMinus = edge2elem(interior, 2);
locPlus = sum((edgeIdx(tPlus, :) == interior) .* (1:3), 2);
locMinus = sum((edgeIdx(tMinus, :) == interior) .* (1:3), 2);

xA = nodeH(edge(interior, 1), :);
xB = nodeH(edge(interior, 2), :);
tau = xB - xA;
hEdge = sqrt(sum(tau.^2, 2));
normal = [tau(:,2), -tau(:,1)] ./ hEdge;
midpoint = 0.5 * (xA + xB);

gradAll = baryGradients2DAll(nodeH, baseElem);
gradPlus = gradAll(tPlus, :, :);
gradMinus = gradAll(tMinus, :, :);
alphaPlus = gradPlus(:,:,1) .* normal(:,1) + gradPlus(:,:,2) .* normal(:,2);
alphaMinus = gradMinus(:,:,1) .* normal(:,1) + gradMinus(:,:,2) .* normal(:,2);

localEdgePairs = [1 2; 2 3; 3 1];
pairPlus = localEdgePairs(locPlus, :);
pairMinus = localEdgePairs(locMinus, :);
globalEdge = edge(interior, :);
monos = lagrangeMonomials2D(degree);

jump = zeros(nEdge, nQuad, nTraceDof);
for q = 1:nQuad
    lamPlus = edgeLambdaVector(baseElem(tPlus, :), pairPlus, globalEdge, tq(q));
    lamMinus = edgeLambdaVector(baseElem(tMinus, :), pairMinus, globalEdge, tq(q));
    dPlus = directionalDerivativeMonomialsVector(monos, lamPlus, alphaPlus, order);
    dMinus = directionalDerivativeMonomialsVector(monos, lamMinus, alphaMinus, order);
    jump(:, q, :) = [dPlus, -dMinus];
end

jumpData.dof = [elemH(tPlus, :), elemH(tMinus, :)];
jumpData.jump = jump;
jumpData.hEdge = hEdge;
jumpData.midpoint = midpoint;
end


function G = baryGradients2DAll(node, elem)
x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
G = zeros(size(elem, 1), 3, 2);
G(:,1,1) = (y2 - y3) ./ area2; G(:,1,2) = (x3 - x2) ./ area2;
G(:,2,1) = (y3 - y1) ./ area2; G(:,2,2) = (x1 - x3) ./ area2;
G(:,3,1) = (y1 - y2) ./ area2; G(:,3,2) = (x2 - x1) ./ area2;
end


function lam = edgeLambdaVector(elemVerts, pair, globalEdge, t)
nEdge = size(elemVerts, 1);
lam = zeros(nEdge, 3);
row = (1:nEdge).';
p1 = pair(:,1);
p2 = pair(:,2);
localFirst = elemVerts(sub2ind([nEdge, 3], row, p1));
sameOrientation = localFirst == globalEdge(:,1);
v1 = sameOrientation * (1 - t) + (~sameOrientation) * t;
v2 = sameOrientation * t + (~sameOrientation) * (1 - t);
lam(sub2ind([nEdge, 3], row, p1)) = v1;
lam(sub2ind([nEdge, 3], row, p2)) = v2;
end


function monos = lagrangeMonomials2D(degree)
switch degree
    case 1
        monos = {[1 0 0 1]; [0 1 0 1]; [0 0 1 1]};
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
end
end


function vals = directionalDerivativeMonomialsVector(monos, lambda, alpha, order)
nEdge = size(lambda, 1);
nLB = numel(monos);
vals = zeros(nEdge, nLB);
betas = derivativeMultiIndices(order);
multiCoef = multinomialRows(betas);
for a = 1:nLB
    terms = monos{a};
    val = zeros(nEdge, 1);
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
            val = val + c * multiCoef(b) * fall .* ...
                prod(alpha.^beta, 2) .* prod(lambda.^remPow, 2);
        end
    end
    vals(:, a) = val;
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
