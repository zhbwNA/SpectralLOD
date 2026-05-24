function B = assembleNedAdvection2D(node, elem, beta, opts)
% ASSEMBLENEDADVECTION2D  Assemble int w_i . ((beta.grad) w_j) for NE_1 2D.

if nargin < 3 || isempty(beta), beta = [0, 0]; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = 2;
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
    [beta1, beta2] = evalAdvectionCoefficient(beta, xq, yq);

    phix = [lq(2)*g3(:,1) - lq(3)*g2(:,1), ...
            lq(3)*g1(:,1) - lq(1)*g3(:,1), ...
            lq(1)*g2(:,1) - lq(2)*g1(:,1)];
    phiy = [lq(2)*g3(:,2) - lq(3)*g2(:,2), ...
            lq(3)*g1(:,2) - lq(1)*g3(:,2), ...
            lq(1)*g2(:,2) - lq(2)*g1(:,2)];

    bg1 = beta1 .* g1(:,1) + beta2 .* g1(:,2);
    bg2 = beta1 .* g2(:,1) + beta2 .* g2(:,2);
    bg3 = beta1 .* g3(:,1) + beta2 .* g3(:,2);
    dphix = [bg2 .* g3(:,1) - bg3 .* g2(:,1), ...
             bg3 .* g1(:,1) - bg1 .* g3(:,1), ...
             bg1 .* g2(:,1) - bg2 .* g1(:,1)];
    dphiy = [bg2 .* g3(:,2) - bg3 .* g2(:,2), ...
             bg3 .* g1(:,2) - bg1 .* g3(:,2), ...
             bg1 .* g2(:,2) - bg2 .* g1(:,2)];

    S = S + (2 * weight(q) * area) .* sig(:,aa) .* sig(:,bb) .* ...
        (phix(:,aa) .* dphix(:,bb) + phiy(:,aa) .* dphiy(:,bb));
end

B = sparse(reshape(eid(:, aa), [], 1), reshape(eid(:, bb), [], 1), reshape(S, [], 1), NE, NE);
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

error('assembleNedAdvection2D:beta', 'Unsupported advection coefficient type.');
end


function val = evalField(beta, name, x, y, defaultValue)
if isfield(beta, name) && ~isempty(beta.(name))
    val = evalPDECoefficient(beta.(name), x, y, [], []);
else
    val = defaultValue * ones(size(x));
end
end
