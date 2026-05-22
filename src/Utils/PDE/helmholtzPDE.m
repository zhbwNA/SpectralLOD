function pde = helmholtzPDE(k, varargin)
% HELMHOLTZPDE  Create coefficient data for Helmholtz/Maxwell wave assembly.

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
