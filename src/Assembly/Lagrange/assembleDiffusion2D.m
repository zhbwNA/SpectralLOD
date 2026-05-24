function A = assembleDiffusion2D(node, elem, degree, coef, opts)
% ASSEMBLEDIFFUSION2D  Assemble int (D grad u).grad v on triangles.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2 * degree, 2);
end
if isnumeric(coef) && isscalar(coef)
    A = coef * assembleStiffness2D(node, elem, degree);
    return;
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);

[lambda, weight] = quadtriangle(min(6, opts.quadOrder));
[~, DphiRef] = lagrange2D(degree, lambda);

[area, g1x, g1y, g2x, g2y, g3x, g3y, x1, y1, x2, y2, x3, y3] = elementGeometry2D(node, elem);
[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));

for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    [d11, d12, d21, d22] = evalDiffusionCoefficient(coef, xq, yq);

    Dq = squeeze(DphiRef(q, :, :));
    Dx = g1x * Dq(:,1).' + g2x * Dq(:,2).' + g3x * Dq(:,3).';
    Dy = g1y * Dq(:,1).' + g2y * Dq(:,2).' + g3y * Dq(:,3).';

    S = S + (2 * weight(q) * area) .* ...
        ((d11 .* Dx(:,bb) + d12 .* Dy(:,bb)) .* Dx(:,aa) + ...
         (d21 .* Dx(:,bb) + d22 .* Dy(:,bb)) .* Dy(:,aa));
end

A = sparse(reshape(elem(:, aa), [], 1), reshape(elem(:, bb), [], 1), reshape(S, [], 1), N, N);
end


function [d11, d12, d21, d22] = evalDiffusionCoefficient(coef, x, y)
if isnumeric(coef) || isa(coef, 'function_handle')
    c = evalPDECoefficient(coef, x, y, [], []);
    d11 = c; d22 = c;
    d12 = zeros(size(x)); d21 = zeros(size(x));
    return;
end

if ~isstruct(coef)
    error('assembleDiffusion2D:coef', 'Coefficient must be scalar, function handle, or struct.');
end

d11 = evalField(coef, 'd11', x, y, 1);
d22 = evalField(coef, 'd22', x, y, 1);
d12 = evalField(coef, 'd12', x, y, 0);
d21 = evalField(coef, 'd21', x, y, 0);
end


function val = evalField(coef, name, x, y, defaultValue)
if isfield(coef, name) && ~isempty(coef.(name))
    val = evalPDECoefficient(coef.(name), x, y, [], []);
else
    val = defaultValue * ones(size(x));
end
end


function [area, g1x, g1y, g2x, g2y, g3x, g3y, x1, y1, x2, y2, x3, y3] = elementGeometry2D(node, elem)
x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area = abs(area2) / 2;
g1x = (y2 - y3) ./ area2; g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2; g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2; g3y = (x2 - x1) ./ area2;
end
