function [applyPrecon, Aeps, pdeShift] = shiftedLaplacianPreconditionerMaxwell2D(node, elem, bdFlag, k, opts)
% SHIFTEDLAPLACIANPRECONDITIONERMAXWELL2D  Build shifted NE_1 Maxwell solver.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'epsilon'), opts.epsilon = 'quadratic'; end
if ~isfield(opts, 'eta'), opts.eta = 'k'; end
if ~isfield(opts, 'solverMode'), opts.solverMode = 'lu'; end

pdeShift = shiftedLaplacianPDE(k, 'epsilon', opts.epsilon, 'eta', opts.eta);
Aeps = assembleMaxwell2D(node, elem, bdFlag, pdeShift);

if strcmpi(opts.solverMode, 'direct')
    applyPrecon = @(r) Aeps \ r;
else
    [L, U, p, q] = lu(Aeps, 'vector');
    applyPrecon = @(r) applyLU(L, U, p, q, r);
end
end


function x = applyLU(L, U, p, q, r)
x = zeros(size(r));
x(q) = U \ (L \ r(p));
end
