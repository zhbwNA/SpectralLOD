function pde = helmholtzPDE(k, varargin)
% HELMHOLTZPDE  Create coefficient data for Helmholtz/Maxwell wave assembly.
%
%   pde = helmholtzPDE(k, 'source', f, 'boundaryData', g, 'eta', eta)
%   accepts constant, numeric-array, or function-valued k. The scalar
%   Helmholtz assemblers interpret the structure as
%       K - int (k(x)^2+i*epsilon(x)) phi_j phi_i dx
%         - i int_Gamma eta(x) phi_j phi_i ds.
%   Named eta rules are 'k', 'sqrt', and 'zero'; epsilon rules are 'zero',
%   'linear', and 'quadratic'. Numeric values and function handles are also
%   accepted by the coefficient evaluators.

pde = struct();
pde.type = 'helmholtz';
pde.k = k;
pde.epsilon = 0;
pde.eta = 'k';
pde.source = [];
pde.boundaryData = [];
pde.description = 'unshifted Helmholtz operator';

for i = 1:2:numel(varargin)
    name = lower(varargin{i});
    value = varargin{i+1};
    switch name
        case {'epsilon', 'eps'}
            pde.epsilon = value;
        case 'eta'
            pde.eta = value;
        case {'source', 'f'}
            pde.source = value;
        case {'boundarydata', 'g'}
            pde.boundaryData = value;
        case 'description'
            pde.description = value;
        otherwise
            error('helmholtzPDE:option', 'Unknown option "%s".', varargin{i});
    end
end
end
