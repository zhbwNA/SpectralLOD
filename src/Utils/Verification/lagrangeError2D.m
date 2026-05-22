function err = lagrangeError2D(node, elem, degree, uh, uExact, gradExact, k, quadOrder)
% LAGRANGEERROR2D  Compute L2, H1, and k-weighted errors for P1-P3 fields.

if nargin < 8 || isempty(quadOrder), quadOrder = max(4, 2 * degree + 2); end
if nargin < 7 || isempty(k), k = 1; end
if size(elem, 2) == 3 && degree > 1
    [node, elem] = extendMesh2D(node, elem, degree);
end
if numel(uh) ~= size(node, 1)
    error('lagrangeError2D:size', 'uh must have one value per active node.');
end
uh = uh(:);

[lambda, weight] = quadtriangle(min(6, quadOrder));
[phi, DphiRef] = lagrange2D(degree, lambda);
NT = size(elem, 1);
nQuad = numel(weight);

x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area = abs(area2) / 2;

g1x = (y2 - y3) ./ area2; g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2; g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2; g3y = (x2 - x1) ./ area2;

elU = uh(elem);
l2 = 0;
h1 = 0;
for q = 1:nQuad
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    uhq = elU * phi(q, :).';
    uq = uExact(xq, yq);
    du = uhq - uq;

    Dq = squeeze(DphiRef(q, :, :));
    phix = g1x * Dq(:,1).' + g2x * Dq(:,2).' + g3x * Dq(:,3).';
    phiy = g1y * Dq(:,1).' + g2y * Dq(:,2).' + g3y * Dq(:,3).';
    uxq = sum(elU .* phix, 2);
    uyq = sum(elU .* phiy, 2);
    [uexX, uexY] = gradExact(xq, yq);
    dg2 = abs(uxq - uexX).^2 + abs(uyq - uexY).^2;

    wphys = 2 * weight(q) * area;
    l2 = l2 + sum(wphys .* abs(du).^2);
    h1 = h1 + sum(wphys .* dg2);
end

err = struct();
err.L2 = sqrt(real(l2));
err.H1semi = sqrt(real(h1));
err.energy = sqrt(real(h1 + k^2 * l2));
end
