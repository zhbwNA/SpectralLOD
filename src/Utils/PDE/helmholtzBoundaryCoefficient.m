function eta = helmholtzBoundaryCoefficient(pde, x, y, z)
% HELMHOLTZBOUNDARYCOEFFICIENT  Evaluate impedance/absorbing boundary eta.

if nargin < 4, z = []; end
kval = evalPDECoefficient(pde.k, x, y, z, []);
epsval = helmholtzShiftValue(pde.epsilon, x, y, z, kval);

if ischar(pde.eta) || isstring(pde.eta)
    name = lower(char(pde.eta));
    switch name
        case 'k'
            eta = kval;
        case {'sqrtk2ieps', 'sqrt(k2+ieps)', 'sqrt'}
            eta = sqrt(kval.^2 + 1i * epsval);
        case {'none', 'zero'}
            eta = zeros(size(kval));
        otherwise
            error('helmholtzBoundaryCoefficient:eta', 'Unknown eta rule "%s".', name);
    end
else
    eta = evalPDECoefficient(pde.eta, x, y, z, kval);
end
end
