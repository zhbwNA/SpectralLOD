% VERIFY_GGGLS_RAS_PML_REPRODUCE  Section 8.3-style GGGLS RAS/RMS-PML pilot.

fprintf('========== GGGLS RAS-PML/RMS-PML Reproduction Pilot ==========\n\n');

outDir = fullfile(fileparts(mfilename('fullpath')), 'gggls_ras_pml');
if ~exist(outDir, 'dir'), mkdir(outDir); end

opts = ggglsOptions();
fprintf('Degree p=%d, delta=%.5g, kappa=%.5g, alpha=%.5g, tol=%.1e\n', ...
    opts.degree, opts.delta, opts.kappa, opts.pmlAlpha, opts.tol);
fprintf('Pilot k values: %s\n\n', mat2str(opts.kVals));

resultTemplate = struct('method', '', 'shape', '', 'grid', '', 'k', NaN, 'hInv', NaN, ...
    'dof', NaN, 'rasFixed', NaN, 'rasGMRES', NaN, 'rmsFixed', NaN, ...
    'rasRelres', NaN, 'gmresRelres', NaN, 'gmresPrecRelres', NaN, 'rmsRelres', NaN, ...
    'rasHist', NaN, 'rmsHist', NaN, 'status', '');
cases = ggglsCases(opts);
results = repmat(resultTemplate, numel(opts.kVals) * numel(cases), 1);
idxResult = 0;

for k = opts.kVals
    for cfg = cases
        idxResult = idxResult + 1;
        result = runGGGLSCase(k, cfg, opts);
        results(idxResult) = result;
        fprintf('%-12s k=%-4g grid=%-5s DOF=%-7d RAS=%-5s GMRES=%-5s RMS=%-5s\n', ...
            result.shape, result.k, result.grid, result.dof, iterString(result.rasFixed, opts.maxIter), ...
            iterString(result.rasGMRES, opts.maxIter), iterString(result.rmsFixed, opts.maxIter));
    end
end
results = results(1:idxResult);

plotGGGLSResults(outDir, results);
writeGGGLSReport(fullfile(outDir, 'gggls_ras_pml_results.md'), results, opts);
save(fullfile(outDir, 'gggls_ras_pml_results.mat'), 'results', 'opts');

okCases = strcmp({results.status}, 'ok');
assert(any([results(okCases).rasFixed] < opts.maxIter) || any([results(okCases).rmsFixed] < opts.maxIter), ...
    'No GGGLS RAS/RMS pilot case reached tolerance.');
fprintf('\nResults written to %s\n', outDir);
fprintf('========== GGGLS RAS-PML/RMS-PML pilot complete ==========\n');


function opts = ggglsOptions()
opts.degree = 2;
opts.delta = 1 / 40;
opts.kappa = 1 / 40;
opts.pmlAlpha = 5000;
opts.tol = numericEnvScalar('GGGLS_TOL', 1e-6);
opts.maxIter = numericEnvScalar('GGGLS_MAXITER', 80);
opts.kVals = numericEnvList('GGGLS_KVALS', 100);
opts.stripVals = numericEnvList('GGGLS_STRIPS', 2);
opts.checkerVals = numericEnvList('GGGLS_CHECKERS', 2);
opts.solverMode = strtrim(getenvDefault('GGGLS_SOLVER_MODE', 'lu'));
opts.memoryBudgetGB = numericEnvScalar('GGGLS_MEMORY_GB', 100);
end


function cases = ggglsCases(opts)
cases = struct('shape', {}, 'gridSize', {}, 'grid', {});
for n = opts.stripVals
    cases(end+1) = struct('shape', 'strip', 'gridSize', [n, 1], ...
        'grid', sprintf('%d', n)); %#ok<AGROW>
end
for n = opts.checkerVals
    cases(end+1) = struct('shape', 'checker', 'gridSize', [n, n], ...
        'grid', sprintf('%dx%d', n, n)); %#ok<AGROW>
end
end


function result = runGGGLSCase(k, cfg, opts)
n = alignedMeshDivisions(k, cfg.gridSize, opts);
h = 1 / n;
[estDof, estElem] = estimateP2MeshSize(n, opts.kappa);
estGB = estimateCaseMemoryGB(estDof, estElem, prod(cfg.gridSize));
if estGB > opts.memoryBudgetGB
    result = skippedResult(k, cfg, n, estDof, opts);
    result.status = sprintf('skip-memory-est-%.1fGB', estGB);
    return;
end

pbox = [-opts.kappa, 1 + opts.kappa, -opts.kappa, 1 + opts.kappa];
[node, elem] = squaremesh(pbox, h);
[nodeH, elemH] = extendMesh2D(node, elem, opts.degree);

box = struct('physicalBox', [0, 1, 0, 1], 'outerBox', pbox);
asmOpts = struct('pmlAlpha', opts.pmlAlpha, 'quadOrder', 6);
x0 = [0.5, 0.5];
f = @(x,y) besselj(0, k * sqrt((x - x0(1)).^2 + (y - x0(2)).^2));
[Afull, bfull, freeDof] = assembleGGGLSPML2D(nodeH, elemH, k, box, f, opts.degree, asmOpts);
A = Afull(freeDof, freeDof);
b = bfull(freeDof);

parts = ggglsSubdomains2D(node, elem, cfg.gridSize, opts.delta, opts.kappa);
parts = ggglsPartitionOfUnity2D(parts, cfg.gridSize, opts.delta);
precon = ggglsRASPML2D(nodeH, elemH, parts, k, opts.degree, ...
    struct('pmlAlpha', opts.pmlAlpha, 'quadOrder', 6, 'solverMode', opts.solverMode));
applyFree = @(r) restrictApply(precon.applyRAS, r, freeDof, size(nodeH,1));

[rasIts, rasHist] = rasFixedPoint(A, b, precon, freeDof, opts.tol, opts.maxIter);
[gmIts, gmRel, gmPrecRel] = rasGMRES(A, b, applyFree, opts.tol, opts.maxIter);
[rmsIts, rmsHist] = rmsFixedPoint(A, b, precon, cfg.gridSize, freeDof, opts.tol, opts.maxIter);

result = struct('method', 'GGGLS', 'shape', cfg.shape, 'grid', cfg.grid, ...
    'k', k, 'hInv', n, 'dof', size(A,1), 'rasFixed', rasIts, ...
    'rasGMRES', gmIts, 'rmsFixed', rmsIts, ...
    'rasRelres', rasHist(end), 'gmresRelres', gmRel, ...
    'gmresPrecRelres', gmPrecRel, 'rmsRelres', rmsHist(end), ...
    'rasHist', rasHist, 'rmsHist', rmsHist, 'status', 'ok');
end


function n = alignedMeshDivisions(k, gridSize, ~)
base = ceil(k^1.25);
align = 40;
align = lcm(align, gridSize(1));
align = lcm(align, gridSize(2));
n = ceil(base / align) * align;
end


function [nDof, nElem] = estimateP2MeshSize(nPhysical, kappa)
nx = round((1 + 2 * kappa) * nPhysical);
ny = nx;
nP1 = (nx + 1) * (ny + 1);
nEdge = nx * (ny + 1) + (nx + 1) * ny + nx * ny;
nDof = nP1 + nEdge;
nElem = 2 * nx * ny;
end


function gb = estimateCaseMemoryGB(nDof, nElem, nSub)
globalGB = 112 * nDof / 1e9;
elemGB = 128 * nElem / 1e9;
localGB = nSub * 350 * (max(nDof / nSub, 1)^1.35) / 1e9;
gb = globalGB + elemGB + localGB;
end


function result = skippedResult(k, cfg, n, dof, ~)
result = struct('method', 'GGGLS', 'shape', cfg.shape, 'grid', cfg.grid, ...
    'k', k, 'hInv', n, 'dof', dof, 'rasFixed', NaN, ...
    'rasGMRES', NaN, 'rmsFixed', NaN, ...
    'rasRelres', NaN, 'gmresRelres', NaN, 'gmresPrecRelres', NaN, 'rmsRelres', NaN, ...
    'rasHist', NaN, 'rmsHist', NaN, 'status', 'skip');
end


function [its, hist] = rasFixedPoint(A, b, precon, freeDof, tol, maxIter)
u = zeros(size(b));
r0 = norm(b);
hist = zeros(maxIter + 1, 1);
hist(1) = 1;
its = maxIter + 1;
for it = 1:maxIter
    r = b - A * u;
    rel = norm(r) / max(r0, eps);
    hist(it) = rel;
    if rel < tol
        its = it - 1;
        hist = hist(1:it);
        return;
    end
    z = restrictApply(precon.applyRAS, r, freeDof, precon.nGlobal);
    u = u + z;
end
hist(end) = norm(b - A * u) / max(r0, eps);
if hist(end) < tol, its = maxIter; end
end


function [its, trueRel, precRel] = rasGMRES(A, b, applyFree, tol, maxIter)
rhs = applyFree(b);
Afun = @(x) applyFree(A * x);
[x, ~, precRel, iter] = gmres(Afun, rhs, [], tol, maxIter);
trueRel = norm(b - A * x) / max(norm(b), eps);
its = maxIter + 1;
if ~isempty(iter)
    its = iter(2);
end
if precRel >= tol && its <= maxIter
    its = maxIter + 1;
end
end


function [its, hist] = rmsFixedPoint(A, b, precon, gridSize, freeDof, tol, maxIter)
orders = rmsOrders(gridSize);
uFull = zeros(precon.nGlobal, 1);
localU = cell(precon.nSub, 1);
for s = 1:precon.nSub
    localU{s} = zeros(numel(precon.gIdx{s}), 1);
end
r0 = norm(b);
hist = zeros(maxIter + 1, 1);
hist(1) = 1;
its = maxIter + 1;
for it = 1:maxIter
    u = uFull(freeDof);
    r = b - A * u;
    rel = norm(r) / max(r0, eps);
    hist(it) = rel;
    if rel < tol
        its = it - 1;
        hist = hist(1:it);
        return;
    end
    for oo = 1:numel(orders)
        for s = orders{oo}
            u = uFull(freeDof);
            r = b - A * u;
            rr = zeros(precon.nGlobal, 1);
            rr(freeDof) = r;
            cLoc = precon.solveLocalCorrection(rr, s);
            idx = precon.gIdx{s};
            newLoc = uFull(idx) + cLoc;
            uFull(idx) = uFull(idx) + precon.weight{s} .* (newLoc - localU{s});
            localU{s} = newLoc;
        end
    end
end
u = uFull(freeDof);
hist(end) = norm(b - A * u) / max(r0, eps);
if hist(end) < tol, its = maxIter; end
end


function orders = rmsOrders(gridSize)
nx = gridSize(1);
ny = gridSize(2);
if ny == 1
    orders = {1:nx, nx:-1:1};
    return;
end
orders = {
    lexOrder(nx, ny, 1:nx, 1:ny), ...
    lexOrder(nx, ny, nx:-1:1, ny:-1:1), ...
    lexOrder(nx, ny, 1:nx, ny:-1:1), ...
    lexOrder(nx, ny, nx:-1:1, 1:ny)};
end


function order = lexOrder(nx, ny, xOrder, yOrder)
order = zeros(1, nx * ny);
idx = 0;
for j = yOrder
    for i = xOrder
        idx = idx + 1;
        order(idx) = (j - 1) * nx + i;
    end
end
end


function z = restrictApply(applyFull, r, freeDof, nFull)
rr = zeros(nFull, 1);
rr(freeDof) = r;
zz = applyFull(rr);
z = zz(freeDof);
end


function plotGGGLSResults(outDir, results)
fig = figure('Name', 'GGGLS RAS RMS histories', 'Color', 'w');
hold on;
labels = {};
for i = 1:numel(results)
    if strcmp(results(i).status, 'ok')
        semilogy(0:numel(results(i).rasHist)-1, results(i).rasHist, '-', 'LineWidth', 1.2);
        labels{end+1} = sprintf('RAS %s k=%g', results(i).grid, results(i).k); %#ok<AGROW>
        semilogy(0:numel(results(i).rmsHist)-1, results(i).rmsHist, '--', 'LineWidth', 1.2);
        labels{end+1} = sprintf('RMS %s k=%g', results(i).grid, results(i).k); %#ok<AGROW>
    end
end
xlabel('iteration', 'Interpreter', 'latex');
ylabel('$\|b-Au_n\|_2/\|b\|_2$', 'Interpreter', 'latex');
title('GGGLS paper-style RAS-PML/RMS-PML pilot', 'Interpreter', 'latex');
legend(labels, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_gggls_histories.png'));

fig = figure('Name', 'GGGLS iterations', 'Color', 'w');
ok = strcmp({results.status}, 'ok');
r = results(ok);
hold on;
for shape = ["strip", "checker"]
    mask = strcmp({r.shape}, char(shape));
    rr = r(mask);
    plot([rr.k], [rr.rasFixed], 'o-', 'LineWidth', 1.3);
    plot([rr.k], [rr.rmsFixed], 's--', 'LineWidth', 1.3);
end
xlabel('$k$', 'Interpreter', 'latex');
ylabel('fixed-point iterations', 'Interpreter', 'latex');
title('GGGLS paper-style fixed-point counts', 'Interpreter', 'latex');
legend({'RAS strip','RMS strip','RAS checker','RMS checker'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_gggls_iterations.png'));
end


function writeGGGLSReport(fileName, results, opts)
fid = fopen(fileName, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# GGGLS RAS-PML/RMS-PML Reproduction Pilot\n\n');
fprintf(fid, 'This report uses a separate non-divergence assembly for the GGGLS form\n');
fprintf(fid, '`a(u,v)=int k^{-2}((D grad u).grad v - (beta.grad u)v)-u v`, with `f_PML(x)=5000*x^3/3`, `delta=kappa=1/40`, and P2 elements.\n\n');
fprintf(fid, 'The mesh divisor is rounded up to a multiple of `40`, and of the requested subdomain grid sizes, so that the global boundary, physical/PML interfaces, overlap interfaces, and local PML boundaries are resolved by mesh edges as assumed in the GGGLS discrete setup.\n\n');
fprintf(fid, '![Residual histories](fig_gggls_histories.png)\n\n');
fprintf(fid, '![Iteration counts](fig_gggls_iterations.png)\n\n');
fprintf(fid, '| shape | grid | k | 1/h | DOF | RAS fixed | RAS GMRES | RMS fixed | RAS relres | GMRES true relres | GMRES prec relres | RMS relres | status |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
for r = results
    fprintf(fid, '| %s | %s | %.0f | %d | %d | %s | %s | %s | %.3e | %.3e | %.3e | %.3e | %s |\n', ...
        r.shape, r.grid, r.k, r.hInv, r.dof, iterString(r.rasFixed, opts.maxIter), ...
        iterString(r.rasGMRES, opts.maxIter), iterString(r.rmsFixed, opts.maxIter), r.rasRelres, ...
        r.gmresRelres, r.gmresPrecRelres, r.rmsRelres, r.status);
end
fprintf(fid, '\nPaper Section 8.3 uses `k=100:50:350`, strips `N=2,4,8`, and checkerboards `2x2,4x4,8x8`. This script defaults to the smallest paper point and skips any requested case whose memory estimate exceeds `GGGLS_MEMORY_GB`.\n\n');
fprintf(fid, 'Target Case 1 counts from GGGLS Tables 1-4:\n\n');
fprintf(fid, '- RAS strips: `N=2`: 6,4,4,4,4,4 fixed-point counts for k=100..350; GMRES counts are the bracketed values.\n');
fprintf(fid, '- RMS strips: `N=2`: 4,3,2,2,2,2 fixed-point counts.\n');
fprintf(fid, '- RAS checkerboard `2x2`: 8,7,7,6,6,6 fixed-point counts.\n');
fprintf(fid, '- RMS checkerboard `2x2`: 2 fixed-point counts across k=100..350.\n');
fprintf(fid, '\n## Direct comparison for the requested cases\n\n');
fprintf(fid, '| shape | grid | k | paper RAS fixed | measured RAS fixed | paper GMRES | measured GMRES | paper RMS fixed | measured RMS fixed | note |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
for r = results
    [pRAS, pGMRES, pRMS] = paperCase1Target(r.shape, r.grid, r.k);
    note = 'ok';
    if ~strcmp(r.status, 'ok')
        note = r.status;
    elseif (~isnan(pRAS) && abs(r.rasFixed - pRAS) > 1) || ...
            (~isnan(pRMS) && abs(r.rmsFixed - pRMS) > 1)
        note = 'check';
    end
    fprintf(fid, '| %s | %s | %.0f | %s | %s | %s | %s | %s | %s | %s |\n', ...
        r.shape, r.grid, r.k, targetString(pRAS), iterString(r.rasFixed, opts.maxIter), ...
        targetString(pGMRES), iterString(r.rasGMRES, opts.maxIter), ...
        targetString(pRMS), iterString(r.rmsFixed, opts.maxIter), note);
end
fprintf(fid, '\nHere `measured GMRES` is the iteration count for the explicitly preconditioned GGGLS system `B_h^{-1}A_h x=B_h^{-1}b`, matching the GMRES system described in the paper. The adjacent measured table also reports the unpreconditioned physical residual `||b-Ax||/||b||`, which can be slightly larger at the same iteration.\n');
end


function saveFigure(fig, fileName)
set(fig, 'PaperPositionMode', 'auto');
print(fig, fileName, '-dpng', '-r160');
end


function s = iterString(its, maxIter)
if isnan(its)
    s = 'skip';
elseif its > maxIter
    s = sprintf('>%d', maxIter);
else
    s = sprintf('%d', its);
end


function s = targetString(v)
if isnan(v)
    s = 'n/a';
else
    s = sprintf('%d', v);
end
end


function [ras, gmres, rms] = paperCase1Target(shape, grid, k)
ras = NaN; gmres = NaN; rms = NaN;
kVals = [100, 150, 200, 250, 300, 350];
[tf, idx] = ismember(k, kVals);
if ~tf, return; end
if strcmp(shape, 'strip') && strcmp(grid, '2')
    rasVals = [6, 4, 4, 4, 4, 4];
    gmresVals = [6, 4, 4, 4, 4, 4];
    rmsVals = [4, 3, 2, 2, 2, 2];
elseif strcmp(shape, 'checker') && strcmp(grid, '2x2')
    rasVals = [8, 7, 7, 6, 6, 6];
    gmresVals = [7, 6, 6, 5, 5, 5];
    rmsVals = [2, 2, 2, 2, 2, 2];
else
    return;
end
ras = rasVals(idx);
gmres = gmresVals(idx);
rms = rmsVals(idx);
end
end


function value = getenvDefault(name, defaultValue)
value = getenv(name);
if isempty(value), value = defaultValue; end
end


function vals = numericEnvList(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    vals = defaultValue;
    return;
end
if any(strcmpi(txt, {'none', '[]', 'empty'}))
    vals = [];
    return;
end
txt = strrep(txt, ',', ' ');
vals = str2num(txt); %#ok<ST2NM>
if isempty(vals), error('Environment variable %s must be numeric.', name); end
end


function val = numericEnvScalar(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    val = defaultValue;
    return;
end
val = str2double(txt);
if isnan(val), error('Environment variable %s must be numeric.', name); end
end
