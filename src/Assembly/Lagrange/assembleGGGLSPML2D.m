function [A, b, freeDof, bdDof, info] = assembleGGGLSPML2D(node, elem, k, box, f, degree, opts)
% ASSEMBLEGGGLSPML2D  Assemble GGGLS non-divergence Cartesian PML form.

if nargin < 6 || isempty(degree), degree = 1; end
if nargin < 7 || isempty(opts), opts = struct(); end
if nargin < 5, f = []; end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(4, 2 * degree + 1);
end
if ~isfield(opts, 'pmlAlpha') || isempty(opts.pmlAlpha)
    opts.pmlAlpha = 5000;
end
if ~isfield(opts, 'cInvSq') || isempty(opts.cInvSq)
    opts.cInvSq = 1;
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);

[lambda, weight] = quadtriangle(min(6, opts.quadOrder));
[phi, DphiRef] = lagrange2D(degree, lambda);
nQuad = numel(weight);

x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area = abs(area2) / 2;

g1x = (y2 - y3) ./ area2; g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2; g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2; g3y = (x2 - x1) ./ area2;

[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:).'; bb = bb(:).';
S = zeros(NT, numel(aa));
rhs = zeros(NT, nLB);

for q = 1:nQuad
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    [d11, d22, beta1, beta2] = ggglsPMLCoefficients2D(xq, yq, box, opts.pmlAlpha);
    cinv = evalCInvSq(opts.cInvSq, xq, yq);
    Dq = squeeze(DphiRef(q, :, :));

    Dx = g1x * Dq(:,1).' + g2x * Dq(:,2).' + g3x * Dq(:,3).';
    Dy = g1y * Dq(:,1).' + g2y * Dq(:,2).' + g3y * Dq(:,3).';

    stiff = d11 .* Dx(:,aa) .* Dx(:,bb) + d22 .* Dy(:,aa) .* Dy(:,bb);
    firstOrder = phi(q, aa) .* (beta1 .* Dx(:,bb) + beta2 .* Dy(:,bb));
    mass = cinv .* phi(q, aa) .* phi(q, bb);
    wphys = 2 * weight(q) * area;
    S = S + wphys .* (k^(-2) * (stiff - firstOrder) - mass);

    if nargout > 1 && ~isempty(f)
        fq = evalLoad(f, xq, yq);
        rhs = rhs + wphys .* fq .* phi(q, :);
    end
end

ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(S, [], 1);
A = sparse(ii, jj, ss, N, N);

if nargout > 1
    b = accumarray(reshape(elem, [], 1), reshape(rhs, [], 1), [N, 1], @sum, 0);
end

bdDof = boxBoundaryNodes2D(node, box.outerBox);
freeDof = setdiff((1:N).', bdDof(:));
info = struct('box', box, 'degree', degree, 'bdDof', bdDof, 'freeDof', freeDof);
end


function [d11, d22, beta1, beta2] = ggglsPMLCoefficients2D(x, y, box, alpha)
[gp1, gpp1] = oneDimGGGLSProfile(x, box.physicalBox(1), box.physicalBox(2), alpha);
[gp2, gpp2] = oneDimGGGLSProfile(y, box.physicalBox(3), box.physicalBox(4), alpha);
gamma1 = 1 + 1i * gp1;
gamma2 = 1 + 1i * gp2;
d11 = 1 ./ (gamma1.^2);
d22 = 1 ./ (gamma2.^2);
beta1 = (1i * gpp1) ./ (gamma1.^3);
beta2 = (1i * gpp2) ./ (gamma2.^3);
end


function [gp, gpp] = oneDimGGGLSProfile(x, a, b, alpha)
gp = zeros(size(x));
gpp = zeros(size(x));
left = x <= a;
if any(left(:))
    r = a - x(left);
    gp(left) = alpha * r.^2;
    gpp(left) = -2 * alpha * r;
end
right = x >= b;
if any(right(:))
    r = x(right) - b;
    gp(right) = alpha * r.^2;
    gpp(right) = 2 * alpha * r;
end
end


function val = evalCInvSq(cInvSq, x, y)
if isnumeric(cInvSq)
    val = cInvSq * ones(size(x));
else
    val = cInvSq(x, y);
end
end


function val = evalLoad(f, x, y)
if isnumeric(f)
    val = f * ones(size(x));
else
    val = f(x, y);
end
end


function bd = boxBoundaryNodes2D(node, box)
tol = 100 * eps(max(1, max(abs(box))));
onX = abs(node(:,1) - box(1)) <= tol | abs(node(:,1) - box(2)) <= tol;
onY = abs(node(:,2) - box(3)) <= tol | abs(node(:,2) - box(4)) <= tol;
bd = find(onX | onY);
end
