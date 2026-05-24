function A = assembleVectorDiffusion2D(node, elem, degree, coef, nComp, opts)
% ASSEMBLEVECTORDIFFUSION2D  Componentwise vector Lagrange diffusion matrix.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end
if nargin < 5 || isempty(nComp), nComp = 2; end
if nargin < 6 || isempty(opts), opts = struct(); end

A = componentBlockAssembly(@assembleDiffusion2D, node, elem, degree, coef, nComp, opts);
end


function A = componentBlockAssembly(assembler, node, elem, degree, coef, nComp, opts)
if iscell(coef)
    blocks = cell(1, nComp);
    for c = 1:nComp
        blocks{c} = assembler(node, elem, degree, coef{c}, opts);
    end
    A = blkdiag(blocks{:});
else
    S = assembler(node, elem, degree, coef, opts);
    A = kron(speye(nComp), S);
end
end
