function [A, b] = assembleHelmholtz2D(node, elem, bdFlag, k, f, g, degree)
% ASSEMBLEHELMHOLTZ2D  Assemble 2D Helmholtz matrix with impedance BC.

if nargin < 7 || isempty(degree), degree = 1; end
if nargin < 5, f = []; end
if nargin < 6, g = []; end

usePDE = isstruct(k);
if usePDE
    pde = k;
    if isempty(f) && isfield(pde, 'source'), f = pde.source; end
    if isempty(g) && isfield(pde, 'boundaryData'), g = pde.boundaryData; end
end

if degree == 1 || size(elem, 2) > 3
    nodeEval = node;
else
    [nodeEval, ~] = extendMesh2D(node, elem(:,1:3), degree);
end

K = assembleStiffness2D(node, elem, degree);

if usePDE
    qfun = @(x,y) helmholtzVolumeCoefficient(pde, x, y, []);
    etafun = @(x,y) helmholtzBoundaryCoefficient(pde, x, y, []);
    Mq = assembleWeightedMass2D(node, elem, degree, qfun);
    Meta = assembleWeightedBoundaryMass2D(node, elem, bdFlag, degree, etafun);
    A = K - Mq - 1i * Meta;
else
    M = assembleMass2D(node, elem, degree);
    Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree);
    A = K - (k^2) * M - 1i * k * Mb;
end

if nargout > 1
    N = size(A, 1);
    if isempty(f)
        b = zeros(N, 1);
    elseif isnumeric(f)
        if ~exist('M', 'var'), M = assembleMass2D(node, elem, degree); end
        b = M * (f * ones(N, 1));
    else
        if ~exist('M', 'var'), M = assembleMass2D(node, elem, degree); end
        b = M * f(nodeEval(:,1), nodeEval(:,2));
    end
    if ~isempty(g)
        if ~exist('Mb', 'var'), Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree); end
        if isnumeric(g)
            b = b + Mb * (g * ones(N, 1));
        else
            b = b + Mb * g(nodeEval(:,1), nodeEval(:,2));
        end
    end
end
end
