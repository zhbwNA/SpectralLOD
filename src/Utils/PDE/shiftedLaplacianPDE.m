function pde = shiftedLaplacianPDE(k, varargin)
% SHIFTEDLAPLACIANPDE  Create shifted Helmholtz data Aeps=K-(k^2+i eps)M-i eta Mb.
%
%   pde = shiftedLaplacianPDE(k, 'epsilon', 'quadratic', 'eta', 'sqrt')
%   creates the shifted form used by the scalar Helmholtz and NE_1 Maxwell
%   wrappers. The default rules are epsilon=abs(k)^2 and eta=k. The eta
%   rule 'sqrt' means sqrt(k^2+i*epsilon).

epsilon = 'quadratic';
eta = 'k';
args = {};
for i = 1:2:numel(varargin)
    name = lower(varargin{i});
    value = varargin{i+1};
    switch name
        case {'epsilon', 'eps'}
            epsilon = value;
        case 'eta'
            eta = value;
        otherwise
            args = [args, varargin(i:i+1)]; %#ok<AGROW>
    end
end

pde = helmholtzPDE(k, 'epsilon', epsilon, 'eta', eta, args{:});
pde.type = 'shiftedLaplacian';
pde.description = 'shifted Laplacian Helmholtz preconditioner';
end
