function epsval = helmholtzShiftValue(epsilon, x, y, z, kval)
% HELMHOLTZSHIFTVALUE  Evaluate shifted-Laplacian absorption epsilon.

if nargin < 4, z = []; end
if ischar(epsilon) || isstring(epsilon)
    name = lower(char(epsilon));
    switch name
        case {'none', 'zero'}
            epsval = zeros(size(kval));
        case {'linear', 'k'}
            epsval = abs(kval);
        case {'quadratic', 'k2', 'k^2'}
            epsval = abs(kval).^2;
        otherwise
            error('helmholtzShiftValue:epsilon', 'Unknown epsilon rule "%s".', name);
    end
else
    epsval = evalPDECoefficient(epsilon, x, y, z, kval);
end
end
