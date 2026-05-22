function M = assembleWeightedNedMass2D(node, elem, coef)
% ASSEMBLEWEIGHTEDNEDMASS2D  Assemble int coef w_i.w_j for NE_1 2D.

if nargin < 3 || isempty(coef), coef = 1; end
if isnumeric(coef) && isscalar(coef)
    M = coef * assembleNedMass2D(node, elem);
    return;
end

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

[lambda_q, weight] = quadtriangle(2);
x1 = node(elem(:,1), :); x2 = node(elem(:,2), :); x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;
invArea2 = 1 ./ area2;

g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];
bc = [2, 3, 1];
eid = edgeIdx(:, bc);
sig = edgeSign(:, bc);

nEntries = NT * 9 * numel(weight);
ii = zeros(nEntries, 1); jj = zeros(nEntries, 1); ss = zeros(nEntries, 1);
idx = 0;
for q = 1:numel(weight)
    l = lambda_q(q, :);
    xq = l(1)*x1(:,1) + l(2)*x2(:,1) + l(3)*x3(:,1);
    yq = l(1)*x1(:,2) + l(2)*x2(:,2) + l(3)*x3(:,2);
    cq = evalPDECoefficient(coef, xq, yq, [], []);
    phix = [l(2)*g3(:,1) - l(3)*g2(:,1), ...
            l(3)*g1(:,1) - l(1)*g3(:,1), ...
            l(1)*g2(:,1) - l(2)*g1(:,1)];
    phiy = [l(2)*g3(:,2) - l(3)*g2(:,2), ...
            l(3)*g1(:,2) - l(1)*g3(:,2), ...
            l(1)*g2(:,2) - l(2)*g1(:,2)];
    w = 2 * weight(q) * area .* cq;
    for a = 1:3
        for b = 1:3
            dotab = phix(:,a).*phix(:,b) + phiy(:,a).*phiy(:,b);
            nxt = idx + 1; idx = idx + NT;
            ii(nxt:idx) = eid(:,a);
            jj(nxt:idx) = eid(:,b);
            ss(nxt:idx) = sig(:,a) .* sig(:,b) .* w .* dotab;
        end
    end
end
M = sparse(ii(1:idx), jj(1:idx), ss(1:idx), NE, NE);
end
