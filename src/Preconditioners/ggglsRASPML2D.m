function precon = ggglsRASPML2D(nodeH, elemH, parts, k, degree, opts)
% GGGLSRASPML2D  Paper-style RAS-PML/RMS-PML local solvers for GGGLS form.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'solverMode') || isempty(opts.solverMode), opts.solverMode = 'lu'; end
if ~isfield(opts, 'pmlAlpha') || isempty(opts.pmlAlpha), opts.pmlAlpha = 5000; end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder), opts.quadOrder = max(4, 2 * degree + 1); end

N = size(nodeH, 1);
nSub = numel(parts);
locSolvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
freeLocal = cell(nSub, 1);
wgt = cell(nSub, 1);

for j = 1:nSub
    eIdx = parts(j).elemIdx;
    gIdx{j} = unique(elemH(eIdx, :));
    g2l = zeros(N, 1);
    g2l(gIdx{j}) = (1:numel(gIdx{j})).';
    localNode = nodeH(gIdx{j}, :);
    localElem = g2l(elemH(eIdx, :));
    box = struct('physicalBox', parts(j).coreBox, 'outerBox', parts(j).pmlBox);
    localOpts = opts;
    [A_loc, ~, fd] = assembleGGGLSPML2D(localNode, localElem, k, box, [], degree, localOpts);
    freeLocal{j} = fd;
    A_free = A_loc(fd, fd);
    if strcmpi(opts.solverMode, 'direct')
        locSolvers{j} = A_free;
    else
        [L, U, p, q] = lu(A_free, 'vector');
        locSolvers{j} = {L, U, p(:), q(:)};
    end
    raw = parts(j).weightFun(nodeH(gIdx{j}, 1), nodeH(gIdx{j}, 2));
    wgt{j} = max(raw(:), 0);
end

precon = struct();
precon.applyRAS = @applyRAS;
precon.applySubdomain = @applySubdomain;
precon.solveLocalCorrection = @solveLocalCorrection;
precon.nSub = nSub;
precon.nGlobal = N;
precon.parts = parts;
precon.gIdx = gIdx;
precon.weight = wgt;

    function x = applyRAS(r)
        x = zeros(N, 1);
        for s = 1:nSub
            x = x + applySubdomain(r, s);
        end
    end

    function x = applySubdomain(r, s)
        x = zeros(N, 1);
        zloc = solveLocalCorrection(r, s);
        x(gIdx{s}) = zloc .* wgt{s};
    end

    function zloc = solveLocalCorrection(r, s)
        rloc = r(gIdx{s});
        zloc = zeros(size(rloc));
        fd = freeLocal{s};
        if strcmpi(opts.solverMode, 'direct')
            zloc(fd) = locSolvers{s} \ rloc(fd);
        else
            S = locSolvers{s};
            zfd = zeros(numel(fd), 1);
            zfd(S{4}) = S{2} \ (S{1} \ rloc(fd(S{3})));
            zloc(fd) = zfd;
        end
    end
end
