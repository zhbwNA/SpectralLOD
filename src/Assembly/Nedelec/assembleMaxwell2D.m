function [A, info] = assembleMaxwell2D(node, elem, bdFlag, k, varargin)
% ASSEMBLEMAXWELL2D  Assemble NE_1 curl-curl - k^2 mass - i eta boundary mass.
usePDE = isstruct(k);
if usePDE
    pde = k;
else
    pde = helmholtzPDE(k);
end

C = assembleCurlCurl2D(node, elem);
qfun = @(x,y) helmholtzVolumeCoefficient(pde, x, y, []);
etafun = @(x,y) helmholtzBoundaryCoefficient(pde, x, y, []);
Mq = assembleWeightedNedMass2D(node, elem, qfun);
Meta = assembleWeightedNedBndMass2D(node, elem, bdFlag, etafun);
A = C - Mq - 1i * Meta;

if nargout > 1
    info = struct('curlCurl', C, 'massShift', Mq, 'boundaryAbsorption', Meta, 'pde', pde);
end
end
