function [A, b] = assembleHelmholtz3D(node, elem, bdFlag, k, f, g, degree)
% ASSEMBLEHELMHOLTZ3D  Assemble 3D Helmholtz matrix with impedance BC.

if nargin < 7 || isempty(degree), degree = 1; end
if nargin < 5, f = []; end
if nargin < 6, g = []; end

usePDE = isstruct(k);
if usePDE
    pde = k;
    if isempty(f) && isfield(pde, 'source'), f = pde.source; end
    if isempty(g) && isfield(pde, 'boundaryData'), g = pde.boundaryData; end
end

K = assembleStiffness3D(node, elem, degree);

if usePDE
    qfun = @(x,y,z) helmholtzVolumeCoefficient(pde, x, y, z);
    etafun = @(x,y,z) helmholtzBoundaryCoefficient(pde, x, y, z);
    Mq = assembleWeightedMass3D(node, elem, degree, qfun);
    MbEta = assembleWeightedBoundaryMass3D(node, elem, bdFlag, degree, etafun);
    A = K - Mq - 1i * MbEta;
else
    M = assembleMass3D(node, elem, degree);
    Mb = assembleBoundaryMass3D(node, elem, bdFlag, degree);
    A = K - (k^2) * M - 1i * k * Mb;
end

if nargout > 1
    N = size(A, 1);
    if isempty(f)
        b = zeros(N, 1);
    elseif isnumeric(f)
        if ~exist('M', 'var'), M = assembleMass3D(node, elem, degree); end
        b = M * (f * ones(N, 1));
    else
        if ~exist('M', 'var'), M = assembleMass3D(node, elem, degree); end
        b = M * f(node(:,1), node(:,2), node(:,3));
    end
    if ~isempty(g)
        if ~exist('Mb', 'var'), Mb = assembleBoundaryMass3D(node, elem, bdFlag, degree); end
        if isnumeric(g)
            b = b + Mb * (g * ones(N, 1));
        else
            b = b + Mb * g(node(:,1), node(:,2), node(:,3));
        end
    end
end
end
