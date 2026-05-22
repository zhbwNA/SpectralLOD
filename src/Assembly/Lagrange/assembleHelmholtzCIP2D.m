function [A, b, C] = assembleHelmholtzCIP2D(node, elem, bdFlag, k, f, g, degree, gamma, opts)
% ASSEMBLEHELMHOLTZCIP2D  Helmholtz impedance matrix plus CIP stabilization.
%
%   [A,b,C] = ASSEMBLEHELMHOLTZCIP2D(node,elem,bdFlag,k,f,g,degree,gamma,opts)

if nargin < 7 || isempty(degree), degree = 1; end
if nargin < 8, gamma = []; end
if nargin < 9, opts = struct(); end

[A0, b] = assembleHelmholtz2D(node, elem, bdFlag, k, f, g, degree);
C = assembleCIP2D(node, elem, degree, gamma, opts);
A = A0 + C;
end
