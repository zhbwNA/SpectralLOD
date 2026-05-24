function B = assembleVectorAdvection2D(node, elem, degree, beta, nComp, opts)
% ASSEMBLEVECTORADVECTION2D  Componentwise vector Lagrange advection matrix.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(beta), beta = [0, 0]; end
if nargin < 5 || isempty(nComp), nComp = 2; end
if nargin < 6 || isempty(opts), opts = struct(); end

B = componentBlockAssembly(@assembleAdvection2D, node, elem, degree, beta, nComp, opts);
end


function B = componentBlockAssembly(assembler, node, elem, degree, coef, nComp, opts)
if iscell(coef)
    blocks = cell(1, nComp);
    for c = 1:nComp
        blocks{c} = assembler(node, elem, degree, coef{c}, opts);
    end
    B = blkdiag(blocks{:});
else
    S = assembler(node, elem, degree, coef, opts);
    B = kron(speye(nComp), S);
end
end
