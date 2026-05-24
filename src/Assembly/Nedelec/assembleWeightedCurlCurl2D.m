function A = assembleWeightedCurlCurl2D(node, elem, coef, opts)
% ASSEMBLEWEIGHTEDCURLCURL2D  Assemble int coef curl w_i curl w_j for NE_1 2D.

if nargin < 3 || isempty(coef), coef = 1; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = 2;
end
if isnumeric(coef) && isscalar(coef)
    A = coef * assembleCurlCurl2D(node, elem);
    return;
end

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

[lambda, weight] = quadtriangle(min(6, opts.quadOrder));
x1 = node(elem(:,1), :); x2 = node(elem(:,2), :); x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;
invArea2 = 1 ./ area2;

g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];

curlVal = zeros(NT, 3);
curlVal(:,1) = 2 * (g2(:,1).*g3(:,2) - g2(:,2).*g3(:,1));
curlVal(:,2) = 2 * (g3(:,1).*g1(:,2) - g3(:,2).*g1(:,1));
curlVal(:,3) = 2 * (g1(:,1).*g2(:,2) - g1(:,2).*g2(:,1));

bc = [2, 3, 1];
eid = edgeIdx(:, bc);
sig = edgeSign(:, bc);

[aa, bb] = ndgrid(1:3, 1:3);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));

for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1)*x1(:,1) + lq(2)*x2(:,1) + lq(3)*x3(:,1);
    yq = lq(1)*x1(:,2) + lq(2)*x2(:,2) + lq(3)*x3(:,2);
    cq = evalPDECoefficient(coef, xq, yq, [], []);
    S = S + (2 * weight(q) * area) .* cq .* ...
        sig(:,aa) .* sig(:,bb) .* curlVal(:,aa) .* curlVal(:,bb);
end

ii = reshape(eid(:, aa), [], 1);
jj = reshape(eid(:, bb), [], 1);
ss = reshape(S, [], 1);
A = sparse(ii, jj, ss, NE, NE);
end
