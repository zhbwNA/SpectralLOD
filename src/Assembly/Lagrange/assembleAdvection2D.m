function B = assembleAdvection2D(node, elem, degree, beta, opts)
% ASSEMBLEADVECTION2D  Assemble int v beta.grad u on triangles.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(beta), beta = [0, 0]; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2 * degree, 2);
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);

[lambda, weight] = quadtriangle(min(6, opts.quadOrder));
[phi, DphiRef] = lagrange2D(degree, lambda);

[area, g1x, g1y, g2x, g2y, g3x, g3y, x1, y1, x2, y2, x3, y3] = elementGeometry2D(node, elem);
[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));

for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    [beta1, beta2] = evalAdvectionCoefficient(beta, xq, yq);

    Dq = squeeze(DphiRef(q, :, :));
    Dx = g1x * Dq(:,1).' + g2x * Dq(:,2).' + g3x * Dq(:,3).';
    Dy = g1y * Dq(:,1).' + g2y * Dq(:,2).' + g3y * Dq(:,3).';

    S = S + (2 * weight(q) * area) .* ...
        phi(q, aa) .* (beta1 .* Dx(:,bb) + beta2 .* Dy(:,bb));
end

B = sparse(reshape(elem(:, aa), [], 1), reshape(elem(:, bb), [], 1), reshape(S, [], 1), N, N);
end


function [beta1, beta2] = evalAdvectionCoefficient(beta, x, y)
if isnumeric(beta)
    if isscalar(beta)
        beta1 = beta * ones(size(x));
        beta2 = zeros(size(x));
    else
        beta1 = beta(1) * ones(size(x));
        beta2 = beta(2) * ones(size(x));
    end
    return;
end

if isstruct(beta)
    beta1 = evalField(beta, 'beta1', x, y, 0);
    beta2 = evalField(beta, 'beta2', x, y, 0);
    return;
end

if isa(beta, 'function_handle')
    try
        [beta1, beta2] = beta(x, y);
    catch
        val = beta(x, y);
        if isscalar(val)
            beta1 = val * ones(size(x));
            beta2 = zeros(size(x));
        elseif size(val, 2) == 2
            beta1 = reshape(val(:,1), size(x));
            beta2 = reshape(val(:,2), size(x));
        else
            beta1 = reshape(val(1,:), size(x));
            beta2 = reshape(val(2,:), size(x));
        end
    end
    if isscalar(beta1), beta1 = beta1 * ones(size(x)); end
    if isscalar(beta2), beta2 = beta2 * ones(size(x)); end
    beta1 = reshape(beta1, size(x));
    beta2 = reshape(beta2, size(x));
    return;
end

error('assembleAdvection2D:beta', 'Unsupported advection coefficient type.');
end


function val = evalField(beta, name, x, y, defaultValue)
if isfield(beta, name) && ~isempty(beta.(name))
    val = evalPDECoefficient(beta.(name), x, y, [], []);
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
