function C = assembleCIP2D(node, elem, degree, gamma, opts)
% ASSEMBLECIP2D  Assemble sum_e sum_j coeff_j int_e [d_n^j u][d_n^j v] ds.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4, gamma = []; end
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'gammaIsCoefficient') || isempty(opts.gammaIsCoefficient)
    opts.gammaIsCoefficient = false;
end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2, degree + 1);
end

if degree < 1 || degree > 3
    error('assembleCIP2D:degree', 'degree must be 1, 2, or 3.');
end

jumpData = normalDerivativeJump2D(node, elem, degree, 1, opts);
N = jumpData.numDof;
nEdge = numel(jumpData.interiorEdge);
nQuad = numel(jumpData.quadWeight);
nTraceDof = size(jumpData.dof, 2);
if nEdge == 0
    C = sparse(N, N);
    return;
end

[aa, bb] = ndgrid(1:nTraceDof, 1:nTraceDof);
aa = aa(:).'; bb = bb(:).';
nPair = numel(aa);
nEntries = nEdge * nQuad * degree * nPair;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for j = 1:degree
    if j == 1
        jumpDataJ = jumpData;
    else
        jumpDataJ = normalDerivativeJump2D(node, elem, degree, j, opts);
    end
    coeff = cipCoefficientVector(gamma, j, jumpDataJ.midpoint, jumpDataJ.hEdge, opts) .* ...
        jumpDataJ.hEdge.^(2*j - 1);

    for q = 1:nQuad
        scale = coeff .* jumpDataJ.hEdge .* jumpDataJ.quadWeight(q);
        jumpQ = squeeze(jumpDataJ.jump(:, q, :));
        rows = idx + (1:(nEdge * nPair));
        ii(rows) = reshape(jumpDataJ.dof(:, aa), [], 1);
        jj(rows) = reshape(jumpDataJ.dof(:, bb), [], 1);
        ss(rows) = reshape(scale .* jumpQ(:, aa) .* jumpQ(:, bb), [], 1);
        idx = rows(end);
    end
end

C = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end


function coeff = cipCoefficientVector(gamma, j, midpoint, hEdge, opts)
if isempty(gamma)
    defaults = [0.1, 0.01, 0.001];
    g = defaults(j) * ones(size(hEdge));
elseif isa(gamma, 'function_handle')
    try
        g = gamma(j, midpoint, hEdge);
        if isscalar(g), g = g * ones(size(hEdge)); end
        g = reshape(g, size(hEdge));
    catch
        g = arrayfun(@(r) gamma(j, midpoint(r, :), hEdge(r)), (1:numel(hEdge)).');
    end
elseif isscalar(gamma)
    g = gamma * ones(size(hEdge));
else
    g = gamma(j) * ones(size(hEdge));
end

if opts.gammaIsCoefficient
    coeff = g;
else
    coeff = 1i * g;
end
end
