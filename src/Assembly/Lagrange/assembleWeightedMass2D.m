function M = assembleWeightedMass2D(node, elem, degree, coef)
% ASSEMBLEWEIGHTEDMASS2D  Assemble int coef phi_i phi_j on triangles.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end

if isnumeric(coef) && isscalar(coef)
    M = coef * assembleMass2D(node, elem, degree);
    return;
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);
quadOrder = max(2 * degree, 2);
[lambda, weight] = quadtriangle(min(6, quadOrder));
[phi, ~] = lagrange2D(degree, lambda);

x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));

[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));
for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    cq = evalPDECoefficient(coef, xq, yq, [], []);
    S = S + (2 * weight(q) * area) .* cq .* phi(q, aa) .* phi(q, bb);
end

ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(S, [], 1);
M = sparse(ii, jj, ss, N, N);
end
