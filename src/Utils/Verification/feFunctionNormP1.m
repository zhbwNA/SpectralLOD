function out = feFunctionNormP1(node, elem, u, s, p, opts)
% FEFUNCTIONNORMP1  Compute sampled W^{s,p} diagnostics for P1 FE functions.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder), opts.quadOrder = 4; end
if ~isfield(opts, 'maxSamples') || isempty(opts.maxSamples), opts.maxSamples = 240; end

dim = size(node, 2);
if dim ~= 2 && dim ~= 3
    error('feFunctionNormP1:dim', 'Only 2D and 3D P1 meshes are supported.');
end
if size(elem, 2) ~= dim + 1
    error('feFunctionNormP1:p1Only', 'Input connectivity must be P1.');
end
if numel(u) ~= size(node, 1)
    error('feFunctionNormP1:size', 'u must have one value per node.');
end
u = u(:);

[lp, samples] = lpNormP1(node, elem, u, p, opts.quadOrder);
[gradLp, gradSamples] = gradNormP1(node, elem, u, p);
frac = 0;

if s > 0 && s < 1
    frac = fractionalPointSeminorm(samples.x, samples.w, samples.u, s, p, dim, opts.maxSamples);
elseif s > 1
    frac = fractionalPointSeminorm(gradSamples.x, gradSamples.w, gradSamples.g, s - 1, p, dim, opts.maxSamples);
end

if s == 0
    total = lp;
elseif s <= 1
    total = combineNorms([lp, gradLp * double(s == 1), frac], p);
else
    total = combineNorms([lp, gradLp, frac], p);
end

out = struct('total', total, 'lp', lp, 'gradLp', gradLp, ...
    'frac', frac, 's', s, 'p', p, 'sampled', s ~= 0 && s ~= 1);
end


function [val, samples] = lpNormP1(node, elem, u, p, quadOrder)
dim = size(node, 2);
if dim == 2
    [lambda, weight] = quadtriangle(quadOrder);
else
    [lambda, weight] = quadtet(quadOrder);
end
nQuad = numel(weight);
NT = size(elem, 1);

if dim == 2
    jac = simplexMeasureScale2D(node, elem);
else
    jac = simplexMeasureScale3D(node, elem);
end

xAll = zeros(NT * nQuad, dim);
uAll = zeros(NT * nQuad, 1);
wAll = zeros(NT * nQuad, 1);
idx = 0;
acc = 0;
maxv = 0;
for q = 1:nQuad
    ids = idx + (1:NT);
    xq = zeros(NT, dim);
    for d = 1:dim
        coord = reshape(node(elem(:), d), size(elem));
        xq(:, d) = coord * lambda(q, :).';
    end
    uq = u(elem) * lambda(q, :).';
    wq = jac * weight(q);
    xAll(ids, :) = xq;
    uAll(ids) = uq;
    wAll(ids) = wq;
    if isinf(p)
        maxv = max(maxv, max(abs(uq)));
    else
        acc = acc + sum(wq .* abs(uq).^p);
    end
    idx = ids(end);
end

if isinf(p)
    val = maxv;
else
    val = acc^(1 / p);
end
samples = struct('x', xAll, 'u', uAll, 'w', wAll);
end


function [val, samples] = gradNormP1(node, elem, u, p)
dim = size(node, 2);
NT = size(elem, 1);
if dim == 2
    [measure, gradLambda, centroid] = p1Geometry2D(node, elem);
else
    [measure, gradLambda, centroid] = p1Geometry3D(node, elem);
end

grad = zeros(NT, dim);
for a = 1:(dim + 1)
    grad = grad + u(elem(:, a)) .* gradLambda(:, :, a);
end
gmag = sqrt(sum(abs(grad).^2, 2));
if isinf(p)
    val = max(gmag);
else
    val = sum(measure .* gmag.^p)^(1 / p);
end
samples = struct('x', centroid, 'g', grad, 'w', measure);
end


function val = fractionalPointSeminorm(x, w, y, s, p, dim, maxSamples)
n = size(x, 1);
if n > maxSamples
    totalWeight = sum(w);
    pick = unique(round(linspace(1, n, maxSamples)));
    x = x(pick, :);
    w = w(pick);
    y = y(pick, :);
    scale = totalWeight / max(eps, sum(w));
    w = w * scale;
    n = numel(pick);
end

if n < 2
    val = 0;
    return;
end

if isinf(p)
    acc = 0;
    for i = 1:n-1
        dx = x(i+1:n, :) - x(i, :);
        dy = y(i+1:n, :) - y(i, :);
        r = sqrt(sum(dx.^2, 2));
        jump = sqrt(sum(abs(dy).^2, 2));
        acc = max(acc, max(jump ./ max(r.^s, eps)));
    end
    val = acc;
    return;
end

acc = 0;
for i = 1:n-1
    dx = x(i+1:n, :) - x(i, :);
    dy = y(i+1:n, :) - y(i, :);
    r = sqrt(sum(dx.^2, 2));
    jump = sqrt(sum(abs(dy).^2, 2));
    kern = jump.^p ./ max(r.^(dim + s * p), eps);
    acc = acc + 2 * w(i) * sum(w(i+1:n) .* kern);
end
val = acc^(1 / p);
end


function v = combineNorms(parts, p)
parts = parts(isfinite(parts) & parts > 0);
if isempty(parts)
    v = 0;
elseif isinf(p)
    v = max(parts);
else
    v = sum(parts.^p)^(1 / p);
end
end


function jac = simplexMeasureScale2D(node, elem)
x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
jac = abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));
end


function jac = simplexMeasureScale3D(node, elem)
NT = size(elem, 1);
jac = zeros(NT, 1);
for t = 1:NT
    v = node(elem(t, :), :);
    jac(t) = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:); v(4,:) - v(1,:)]));
end
end


function [area, gradLambda, centroid] = p1Geometry2D(node, elem)
NT = size(elem, 1);
x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area = abs(area2) / 2;
gradLambda = zeros(NT, 2, 3);
gradLambda(:, :, 1) = [(y2 - y3) ./ area2, (x3 - x2) ./ area2];
gradLambda(:, :, 2) = [(y3 - y1) ./ area2, (x1 - x3) ./ area2];
gradLambda(:, :, 3) = [(y1 - y2) ./ area2, (x2 - x1) ./ area2];
centroid = (node(elem(:,1), :) + node(elem(:,2), :) + node(elem(:,3), :)) / 3;
end


function [vol, gradLambda, centroid] = p1Geometry3D(node, elem)
NT = size(elem, 1);
vol = zeros(NT, 1);
gradLambda = zeros(NT, 3, 4);
centroid = zeros(NT, 3);
for t = 1:NT
    v = node(elem(t, :), :);
    B = [ones(4,1), v];
    C = inv(B);
    gradLambda(t, :, :) = reshape(C(2:4, :), 1, 3, 4);
    vol(t) = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:); v(4,:) - v(1,:)])) / 6;
    centroid(t, :) = mean(v, 1);
end
end
