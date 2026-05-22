function M = assembleWeightedMass3D(node, elem, degree, coef)
% ASSEMBLEWEIGHTEDMASS3D  Assemble int coef phi_i phi_j on tetrahedra.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end

if isnumeric(coef) && isscalar(coef)
    M = coef * assembleMass3D(node, elem, degree);
    return;
end

if degree > 1 && size(elem, 2) == 4
    [node, elem] = extendMesh3D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);
[lambda, weight] = quadtet(max(2 * degree, 2));
[phi, ~] = lagrange3D(degree, lambda);

v1 = node(elem(:,1), :); v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :); v4 = node(elem(:,4), :);
e12 = v2 - v1; e13 = v3 - v1; e14 = v4 - v1;
detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;

[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));
for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1) * v1(:,1) + lq(2) * v2(:,1) + lq(3) * v3(:,1) + lq(4) * v4(:,1);
    yq = lq(1) * v1(:,2) + lq(2) * v2(:,2) + lq(3) * v3(:,2) + lq(4) * v4(:,2);
    zq = lq(1) * v1(:,3) + lq(2) * v2(:,3) + lq(3) * v3(:,3) + lq(4) * v4(:,3);
    cq = evalPDECoefficient(coef, xq, yq, zq, []);
    S = S + (6 * weight(q) * volume) .* cq .* phi(q, aa) .* phi(q, bb);
end

ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(S, [], 1);
M = sparse(ii, jj, ss, N, N);
end
