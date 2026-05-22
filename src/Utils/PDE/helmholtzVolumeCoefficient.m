function q = helmholtzVolumeCoefficient(pde, x, y, z)
% HELMHOLTZVOLUMECOEFFICIENT  Evaluate k(x)^2+i epsilon(x) for shifted assembly.

if nargin < 4, z = []; end
kval = evalPDECoefficient(pde.k, x, y, z, []);
epsval = helmholtzShiftValue(pde.epsilon, x, y, z, kval);
q = kval.^2 + 1i * epsval;
end
