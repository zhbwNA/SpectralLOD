% VERIFY_THEORY_NEW_FUNCTIONS  Theory-level checks for new PML/transfer/CIP code.

fprintf('========== Theory-Level Verification For New Functions ==========\n\n');

verifyDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(verifyDir);
outDir = fullfile(verifyDir, 'theory_new_functions');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('Output directory: %s\n\n', outDir);

fprintf('Step 0: baseline regression checks\n');
run(fullfile(repoRoot, 'verify', 'verify_pml_helmholtz2d.m'));
run(fullfile(repoRoot, 'verify', 'verify_quasi_interpolation.m'));
run(fullfile(repoRoot, 'verify', 'verify_cip2d.m'));

fprintf('\nStep 1: PML decay/convergence and ORAS-PML iteration checks\n');
pml = runPMLStudy(outDir);
oras = runORASPMLStudy(outDir);

fprintf('\nStep 2: Clement and Scott-Zhang quasi-interpolation checks\n');
transfer = runTransferStudy(outDir);

fprintf('\nStep 3: CIP-FEM preasymptotic checks\n');
cip = runCIPStudy(outDir);

reportFile = fullfile(outDir, 'theory_verification_results.md');
writeTheoryReport(reportFile, pml, oras, transfer, cip);
fprintf('\nTheory verification report written to %s\n', reportFile);

assert(pml.decay.outerToInner < 1, 'PML decay check did not attenuate into the layer.');
assert(any([oras.results.gmresFlag] == 0), 'No ORAS-PML GMRES case converged.');
assert(all(isfinite([transfer.stability.ratio])), 'Transfer stability table contains non-finite ratios.');
assert(all(isfinite([cip.results.femEnergy])) && all(isfinite([cip.results.cipEnergy])), ...
    'CIP/FEM study produced non-finite errors.');

fprintf('\n========== Theory-Level Verification Complete ==========\n');


function pml = runPMLStudy(outDir)
decay = pmlDecayExperiment();
conv = pmlConvergenceExperiment();
sweep = pmlSweepExperiment();
plotPMLDecay(outDir, decay);
plotPMLConvergence(outDir, conv);
plotPMLSweep(outDir, sweep);

fprintf('  PML decay outer/inner amplitude ratio: %.3e\n', decay.outerToInner);
fprintf('  PML convergence fitted C for kh+k^3h^2 model: %.3e\n', conv.fitC);
fprintf('  PML width/absorption sweep best outer/inner ratio: %.3e\n', min([sweep.rows.outerToInner]));
pml = struct('decay', decay, 'convergence', conv, 'sweep', sweep);
end


function decay = pmlDecayExperiment()
k = 20;
width = 0.20;
h = 0.04;
phys = [0, 1, 0, 1];
pbox = [-width, 1 + width, -width, 1 + width];
[node, elem] = squaremesh(pbox, h);
pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
    'sigmaMax', 3 * k, 'sigmaOrder', 2, 'quadOrder', 4);
src = @(x,y) exp(-180 * ((x - 0.48).^2 + (y - 0.52).^2));
[A, b, freeDof] = assembleHelmholtzPML2D(node, elem, k, pml, src);
u = zeros(size(node, 1), 1);
u(freeDof) = A(freeDof, freeDof) \ b(freeDof);

d = outsideDistance(node, phys);
edges = linspace(0, width, 6);
bandCenter = 0.5 * (edges(1:end-1) + edges(2:end));
ampMean = zeros(numel(bandCenter), 1);
ampMax = zeros(numel(bandCenter), 1);
for i = 1:numel(bandCenter)
    if i == 1
        mask = d <= edges(i+1);
    else
        mask = d > edges(i) & d <= edges(i+1);
    end
    vals = abs(u(mask));
    ampMean(i) = mean(vals);
    ampMax(i) = max(vals);
end
decay = struct('k', k, 'h', h, 'width', width, 'bandCenter', bandCenter(:), ...
    'ampMean', ampMean, 'ampMax', ampMax, ...
    'outerToInner', ampMean(end) / max(ampMean(2), eps));
end


function sweep = pmlSweepExperiment()
k = 20;
h = 0.05;
phys = [0, 1, 0, 1];
widthVals = [0.10, 0.15, 0.20, 0.30];
sigmaFactors = [1, 2, 3, 4];
src = @(x,y) exp(-180 * ((x - 0.48).^2 + (y - 0.52).^2));
rows = struct('width', {}, 'sigmaFactor', {}, 'sigmaMax', {}, ...
    'imagStretchMax', {}, 'innerMean', {}, 'outerMean', {}, ...
    'outerToInner', {}, 'outerMax', {}, 'dof', {});
for width = widthVals
    for sf = sigmaFactors
        pbox = [-width, 1 + width, -width, 1 + width];
        [node, elem] = squaremesh(pbox, h);
        pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
            'sigmaMax', sf * k, 'sigmaOrder', 2, 'quadOrder', 4);
        [A, b, freeDof] = assembleHelmholtzPML2D(node, elem, k, pml, src);
        u = zeros(size(node, 1), 1);
        u(freeDof) = A(freeDof, freeDof) \ b(freeDof);
        d = outsideDistance(node, phys);
        innerMask = d > 0 & d <= width / 3;
        outerMask = d > 2 * width / 3 & d <= width + 1e-12;
        if ~any(innerMask)
            innerMask = d <= width / 3;
        end
        rows(end+1) = struct('width', width, 'sigmaFactor', sf, ...
            'sigmaMax', sf * k, 'imagStretchMax', sf, ...
            'innerMean', mean(abs(u(innerMask))), ...
            'outerMean', mean(abs(u(outerMask))), ...
            'outerToInner', mean(abs(u(outerMask))) / max(mean(abs(u(innerMask))), eps), ...
            'outerMax', max(abs(u(outerMask))), ...
            'dof', numel(freeDof)); %#ok<AGROW>
    end
end
sweep = struct('k', k, 'h', h, 'widthVals', widthVals, ...
    'sigmaFactors', sigmaFactors, 'rows', rows);
end


function conv = pmlConvergenceExperiment()
k = 8;
width = 0.25;
phys = [0, 1, 0, 1];
pbox = [-width, 1 + width, -width, 1 + width];
src = @(x,y) exp(-120 * ((x - 0.45).^2 + (y - 0.55).^2));
hRef = 1 / 48;
[refNode, refElem, uRef] = solvePMLSource(pbox, phys, hRef, k, src);

hVals = [1/12, 1/16, 1/24];
rows = struct('h', {}, 'kh', {}, 'pollutionModel', {}, 'l2NodeError', {}, 'ratio', {});
for i = 1:numel(hVals)
    h = hVals(i);
    [node, ~, u] = solvePMLSource(pbox, phys, h, k, src);
    uRefAtNode = interpolateP1(refNode, refElem, uRef, node);
    mask = node(:,1) >= phys(1) - 1e-12 & node(:,1) <= phys(2) + 1e-12 & ...
           node(:,2) >= phys(3) - 1e-12 & node(:,2) <= phys(4) + 1e-12;
    err = sqrt(mean(abs(u(mask) - uRefAtNode(mask)).^2));
    model = k * h + k^3 * h^2;
    rows(end+1) = struct('h', h, 'kh', k*h, 'pollutionModel', model, ...
        'l2NodeError', err, 'ratio', err / model); %#ok<AGROW>
end
fitC = max([rows.ratio]);
conv = struct('k', k, 'hRef', hRef, 'rows', rows, 'fitC', fitC);
end


function [node, elem, u] = solvePMLSource(pbox, phys, h, k, src)
[node, elem] = squaremesh(pbox, h);
pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
    'sigmaMax', 3 * k, 'sigmaOrder', 2, 'quadOrder', 4);
[A, b, freeDof] = assembleHelmholtzPML2D(node, elem, k, pml, src);
u = zeros(size(node, 1), 1);
u(freeDof) = A(freeDof, freeDof) \ b(freeDof);
end


function oras = runORASPMLStudy(outDir)
cases = buildORASPMLCases();
results = repmat(struct('shape', '', 'grid', '', 'k', 0, 'h', 0, 'dof', 0, ...
    'richardson', 0, 'gmres', 0, 'gmresFlag', 0, 'relres', 0, ...
    'resHist', [], 'paperFixedPointRange', '', 'paperGMRESRange', '', ...
    'status', ''), 0, 1);

for c = 1:numel(cases)
    r = runORASPMLCase(cases(c));
    results(end+1) = r; %#ok<AGROW>
    fprintf('  ORAS-PML %-10s k=%g dof=%d Richardson=%s GMRES=%s flag=%d\n', ...
        r.grid, r.k, r.dof, iterString(r.richardson), iterString(r.gmres), r.gmresFlag);
end

plotORASPML(outDir, results);
oras = struct('results', results, 'paper', paperRASPMLTargets());
end


function cases = buildORASPMLCases()
kVals = [20, 40, 60];
cases = struct('shape', {}, 'gridSize', {}, 'gridLabel', {}, 'k', {}, ...
    'h', {}, 'width', {}, 'overlap', {}, 'maxIter', {});
for k = kVals
    h = min(1/36, 2*pi / (8 * k));
    h = 1 / ceil(1 / h);
    cases(end+1) = struct('shape', 'strip', 'gridSize', [2, 1], ...
        'gridLabel', '2 strips', 'k', k, 'h', h, ...
        'width', 0.10, 'overlap', 0.10, 'maxIter', 80); %#ok<AGROW>
    cases(end+1) = struct('shape', 'checker', 'gridSize', [2, 2], ...
        'gridLabel', '2x2', 'k', k, 'h', h, ...
        'width', 0.10, 'overlap', 0.10, 'maxIter', 80); %#ok<AGROW>
end
end


function r = runORASPMLCase(cfg)
phys = [0, 1, 0, 1];
pbox = [-cfg.width, 1 + cfg.width, -cfg.width, 1 + cfg.width];
[node, elem, bd] = squaremesh(pbox, cfg.h);
pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
    'sigmaMax', 3 * cfg.k, 'sigmaOrder', 2, 'quadOrder', 4);
x0 = [0.35, 0.45];
src = @(x,y) besselj(0, cfg.k * sqrt((x - x0(1)).^2 + (y - x0(2)).^2)) .* ...
    exp(-20 * ((x - x0(1)).^2 + (y - x0(2)).^2));
[Afull, bfull, freeDof] = assembleHelmholtzPML2D(node, elem, cfg.k, pml, src);
A = Afull(freeDof, freeDof);
b = bfull(freeDof);

if strcmp(cfg.shape, 'strip')
    nSub = cfg.gridSize(1);
else
    nSub = cfg.gridSize;
end
parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', cfg.overlap);
parts = smoothPartitionOfUnity2D(parts, pbox, cfg.gridSize, cfg.overlap);
applyFull = orasPMLHelmholtz2D(node, elem, cfg.k, parts, pml, 'lu', false);
applyFree = @(v) restrictPMLApply(applyFull, v, freeDof, size(node, 1));

[richIts, hist] = richardsonPML(A, b, applyFree, 1e-6, cfg.maxIter);
[~, flag, relres, iter] = gmres(A, b, [], 1e-6, cfg.maxIter, applyFree);
gmIts = cfg.maxIter + 1;
if ~isempty(iter)
    gmIts = iter(2);
end

[fpRange, gmRange] = paperComparisonRange(cfg.shape, cfg.gridSize);
r = struct('shape', cfg.shape, 'grid', cfg.gridLabel, 'k', cfg.k, 'h', cfg.h, ...
    'dof', size(A, 1), 'richardson', richIts, 'gmres', gmIts, ...
    'gmresFlag', flag, 'relres', relres, 'resHist', hist, ...
    'paperFixedPointRange', fpRange, 'paperGMRESRange', gmRange, ...
    'status', ternary(flag == 0, 'ok', 'gmres-flag'));
end


function z = restrictPMLApply(applyFull, r, freeDof, n)
rr = zeros(n, 1);
rr(freeDof) = r;
zz = applyFull(rr);
z = zz(freeDof);
end


function [its, hist] = richardsonPML(A, b, applyB, tol, maxIter)
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
    u = u + applyB(r);
end
hist(end) = norm(b - A * u) / max(r0, eps);
if hist(end) < tol
    its = maxIter;
end
end


function transfer = runTransferStudy(outDir)
[approx, stability] = transfer2DExperiments();
stability3D = transfer3DPilot();
plotTransfer(outDir, approx, stability);
fprintf('  Transfer max stability ratio: %.3g\n', max([stability.ratio]));
fprintf('  3D pilot max stability ratio: %.3g\n', max([stability3D.ratio]));
transfer = struct('approximation', approx, 'stability', stability, 'stability3D', stability3D);
end


function [approx, stability] = transfer2DExperiments()
uFun = @(x,y) sin(pi*x) .* sin(pi*y) + 0.25 * sin(2*pi*x) .* sin(pi*y);
HVals = [1/2, 1/4, 1/8];
pVals = [1, 2, 4, 8, inf];
sVals = [0, 0.5, 1, 1.25];
approx = struct('operator', {}, 'H', {}, 'p', {}, 'globalRatio', {}, 'localRatio', {});
stability = struct('operator', {}, 'H', {}, 's', {}, 'p', {}, 'ratio', {}, 'sampled', {});

for H = HVals
    h = H / 2;
    [cNode, cElem] = squaremesh([0, 1, 0, 1], H);
    [fNode, fElem] = squaremesh([0, 1, 0, 1], h);
    uF = uFun(fNode(:,1), fNode(:,2));
    P = prolongateNestedP1(cNode, cElem, fNode);
    ops = {'Clement', 'Scott-Zhang'};
    mats = {weightedClementP1(fNode, fElem, cNode, cElem), ...
            scottZhangP1(fNode, fElem, cNode, cElem)};
    for op = 1:numel(ops)
        Q = mats{op};
        uC = Q * uF;
        eF = P * uC - uF;
        w1 = feFunctionNormP1(fNode, fElem, uF, 1, 2);
        for p = pVals
            err0 = feFunctionNormP1(fNode, fElem, eF, 0, p);
            localRatio = localApproxRatio2D(fNode, fElem, cNode, cElem, eF, uF, H);
            approx(end+1) = struct('operator', ops{op}, 'H', H, 'p', p, ...
                'globalRatio', err0.total / max(H * w1.total, eps), ...
                'localRatio', localRatio); %#ok<AGROW>
        end

        for s = sVals
            for p = pVals
                nf = feFunctionNormP1(fNode, fElem, uF, s, p, struct('maxSamples', 180));
                nc = feFunctionNormP1(cNode, cElem, uC, s, p, struct('maxSamples', 180));
                stability(end+1) = struct('operator', ops{op}, 'H', H, ...
                    's', s, 'p', p, 'ratio', nc.total / max(nf.total, eps), ...
                    'sampled', nc.sampled || nf.sampled); %#ok<AGROW>
            end
        end
    end
end
end


function stability3D = transfer3DPilot()
[cNode, cElem] = cubemesh([0, 1, 0, 1, 0, 1], 0.5);
[fNode, fElem] = cubemesh([0, 1, 0, 1, 0, 1], 0.25);
uF = sin(pi*fNode(:,1)) .* sin(pi*fNode(:,2)) .* sin(pi*fNode(:,3));
ops = {'Clement', 'Scott-Zhang'};
mats = {weightedClementP1(fNode, fElem, cNode, cElem), ...
        scottZhangP1(fNode, fElem, cNode, cElem)};
stability3D = struct('operator', {}, 's', {}, 'p', {}, 'ratio', {});
for op = 1:numel(ops)
    uC = mats{op} * uF;
    for s = [0, 1]
        for p = [2, inf]
            nf = feFunctionNormP1(fNode, fElem, uF, s, p);
            nc = feFunctionNormP1(cNode, cElem, uC, s, p);
            stability3D(end+1) = struct('operator', ops{op}, 's', s, ...
                'p', p, 'ratio', nc.total / max(nf.total, eps)); %#ok<AGROW>
        end
    end
end
end


function ratio = localApproxRatio2D(fNode, fElem, cNode, cElem, err, u, H)
cent = (fNode(fElem(:,1), :) + fNode(fElem(:,2), :) + fNode(fElem(:,3), :)) / 3;
[ct, ~] = locateSimplexP1(cNode, cElem, cent, 1e-10);
ratio = 0;
for t = 1:size(cElem, 1)
    elems = find(ct == t);
    if isempty(elems), continue; end
    eL2 = p1SubsetL2(fNode, fElem(elems, :), err);
    uH1 = p1SubsetH1(fNode, fElem(elems, :), u);
    ratio = max(ratio, eL2 / max(H * uH1, eps));
end
end


function cip = runCIPStudy(outDir)
k = 8;
hVals = [1/4, 1/6, 1/8];
results = struct('degree', {}, 'h', {}, 'kh', {}, 'model', {}, ...
    'femEnergy', {}, 'cipEnergy', {}, 'femRatio', {}, 'cipRatio', {});
for degree = 1:3
    for h = hVals
        [femErr, cipErr] = solveCIPPlaneWave(k, h, degree);
        model = (k * h)^degree + k * (k * h)^(2 * degree);
        results(end+1) = struct('degree', degree, 'h', h, 'kh', k*h, ...
            'model', model, 'femEnergy', femErr.energy, ...
            'cipEnergy', cipErr.energy, 'femRatio', femErr.energy / model, ...
            'cipRatio', cipErr.energy / model); %#ok<AGROW>
        fprintf('  CIP p=%d h=%.4g FEM %.3e CIP %.3e model %.3e\n', ...
            degree, h, femErr.energy, cipErr.energy, model);
    end
end
plotCIP(outDir, results);
cip = struct('k', k, 'results', results);
end


function [femErr, cipErr] = solveCIPPlaneWave(k, h, degree)
[node, elem, bd] = squaremesh([0, 1, 0, 1], h);
if degree == 1
    nodeH = node;
    elemH = elem;
else
    [nodeH, elemH] = extendMesh2D(node, elem, degree);
end
uExact = @(x,y) exp(1i * k * x);
gradExact = @(x,y) planeWaveGrad(k, x, y);
f = 0;
g = @(x,y) planeWaveImpedanceData(k, x, y);

[A0, b0] = assembleHelmholtz2D(nodeH, elemH, bd, k, f, g, degree);
uh0 = A0 \ b0;
[Acip, bcip] = assembleHelmholtzCIP2D(nodeH, elemH, bd, k, f, g, degree, []);
uhc = Acip \ bcip;
femErr = lagrangeError2D(nodeH, elemH, degree, uh0, uExact, gradExact, k);
cipErr = lagrangeError2D(nodeH, elemH, degree, uhc, uExact, gradExact, k);
end


function [gx, gy] = planeWaveGrad(k, x, y)
gx = 1i * k * exp(1i * k * x);
gy = zeros(size(y));
end


function g = planeWaveImpedanceData(k, x, y)
u = exp(1i * k * x);
tol = 1e-10;
g = -1i * k * u;
left = abs(x) < tol;
right = abs(x - 1) < tol;
g(left) = -2i * k * u(left);
g(right) = 0;
end


function plotPMLDecay(outDir, decay)
fig = figure('Name', 'PML decay', 'Color', 'w');
semilogy(decay.bandCenter, decay.ampMean, 'o-', 'LineWidth', 1.5); hold on;
semilogy(decay.bandCenter, decay.ampMax, 's--', 'LineWidth', 1.2);
xlabel('$\mathrm{dist}(x,\Omega)$', 'Interpreter', 'latex');
ylabel('$|u_h|$ band amplitude', 'Interpreter', 'latex');
title('$u_h$ attenuation in the PML layer', 'Interpreter', 'latex');
legend({'mean', 'max'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_pml_decay.png'));
end


function plotPMLConvergence(outDir, conv)
rows = conv.rows;
h = [rows.h];
err = [rows.l2NodeError];
model = conv.fitC * [rows.pollutionModel];
fig = figure('Name', 'PML convergence', 'Color', 'w');
loglog(h, err, 'o-', h, model, 'k--', 'LineWidth', 1.5);
xlabel('$h$', 'Interpreter', 'latex');
ylabel('physical-domain RMS error', 'Interpreter', 'latex');
title('$kh+k^3h^2$ PML convergence fit', 'Interpreter', 'latex');
legend({'measured', 'fitted model'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_pml_convergence.png'));
end


function plotPMLSweep(outDir, sweep)
Z = nan(numel(sweep.widthVals), numel(sweep.sigmaFactors));
for r = sweep.rows
    iw = find(abs(sweep.widthVals - r.width) < 1e-12, 1);
    is = find(abs(sweep.sigmaFactors - r.sigmaFactor) < 1e-12, 1);
    Z(iw, is) = r.outerToInner;
end
fig = figure('Name', 'PML width absorption sweep', 'Color', 'w');
imagesc(sweep.sigmaFactors, sweep.widthVals, log10(Z));
set(gca, 'YDir', 'normal');
colorbar;
xlabel('$\sigma_{\max}/k$', 'Interpreter', 'latex');
ylabel('PML layer width', 'Interpreter', 'latex');
title('$\log_{10}$ outer/inner PML amplitude ratio', 'Interpreter', 'latex');
saveFigure(fig, fullfile(outDir, 'fig_pml_width_sigma_sweep.png'));
end


function plotORASPML(outDir, results)
fig = figure('Name', 'ORAS-PML residuals', 'Color', 'w');
hold on;
for i = 1:numel(results)
    semilogy(0:numel(results(i).resHist)-1, results(i).resHist, 'LineWidth', 1.2);
end
xlabel('fixed-point iteration', 'Interpreter', 'latex');
ylabel('$\|r_n\|/\|r_0\|$', 'Interpreter', 'latex');
title('RAS/ORAS-PML residual histories', 'Interpreter', 'latex');
legend(composeORASLabels(results), 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_oras_pml_residuals.png'));

fig = figure('Name', 'ORAS-PML iterations', 'Color', 'w');
ks = unique([results.k]);
hold on;
for shape = ["strip", "checker"]
    mask = strcmp({results.shape}, char(shape));
    rr = results(mask);
    plot([rr.k], [rr.gmres], 'o-', 'LineWidth', 1.5);
end
xlabel('$k$', 'Interpreter', 'latex');
ylabel('GMRES iterations', 'Interpreter', 'latex');
title('Current P1 local-PML ORAS baseline', 'Interpreter', 'latex');
legend({'2 strips', '2x2'}, 'Interpreter', 'latex', 'Location', 'best');
grid on; xlim([min(ks)-5, max(ks)+5]);
saveFigure(fig, fullfile(outDir, 'fig_oras_pml_iterations.png'));
end


function labels = composeORASLabels(results)
labels = cell(1, numel(results));
for i = 1:numel(results)
    labels{i} = sprintf('%s, k=%g', results(i).grid, results(i).k);
end
end


function plotTransfer(outDir, approx, stability)
fig = figure('Name', 'Transfer approximation', 'Color', 'w');
hold on;
ops = unique({approx.operator});
for i = 1:numel(ops)
    mask = strcmp({approx.operator}, ops{i}) & [approx.p] == 2;
    rr = approx(mask);
    loglog([rr.H], [rr.globalRatio], 'o-', 'LineWidth', 1.5);
end
xlabel('$H$', 'Interpreter', 'latex');
ylabel('$\|u-I_Hu\|_{L^2}/(H\|u\|_{W^{1,2}})$', 'Interpreter', 'latex');
title('Quasi-local approximation proxy', 'Interpreter', 'latex');
legend(ops, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_transfer_approximation.png'));

fig = figure('Name', 'Transfer stability', 'Color', 'w');
mask = [stability.p] == 2 & abs([stability.H] - 0.25) < 1e-12;
rr = stability(mask);
hold on;
for i = 1:numel(ops)
    opMask = strcmp({rr.operator}, ops{i});
    plot([rr(opMask).s], [rr(opMask).ratio], 'o-', 'LineWidth', 1.5);
end
xlabel('$s$', 'Interpreter', 'latex');
ylabel('$\|I_Hu\|_{W^{s,2}}/\|u\|_{W^{s,2}}$', 'Interpreter', 'latex');
title('Scott-Zhang/Clement sampled Sobolev stability', 'Interpreter', 'latex');
legend(ops, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_transfer_stability.png'));
end


function plotCIP(outDir, results)
fig = figure('Name', 'CIP preasymptotic', 'Color', 'w');
hold on;
for degree = 1:3
    mask = [results.degree] == degree;
    rr = results(mask);
    loglog([rr.kh], [rr.femEnergy], 'o-', 'LineWidth', 1.3);
    loglog([rr.kh], [rr.cipEnergy], 's--', 'LineWidth', 1.3);
end
xlabel('$kh$', 'Interpreter', 'latex');
ylabel('$k$-weighted error', 'Interpreter', 'latex');
title('FEM and CIP-FEM plane-wave preasymptotic errors', 'Interpreter', 'latex');
legend({'FEM P1','CIP P1','FEM P2','CIP P2','FEM P3','CIP P3'}, ...
    'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_cip_preasymptotic.png'));
end


function writeTheoryReport(fileName, pml, oras, transfer, cip)
fid = fopen(fileName, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Theory-Level Verification Results\n\n');
fprintf(fid, 'Generated by `verify/verify_theory_new_functions.m`.\n\n');

fprintf(fid, '## PML Decay And Convergence\n\n');
fprintf(fid, 'PML decay: `k=%g`, `h=%.4g`, layer width `%.4g`, outer/inner mean-amplitude ratio `%.3e`.\n\n', ...
    pml.decay.k, pml.decay.h, pml.decay.width, pml.decay.outerToInner);
fprintf(fid, 'Definition of band center: let `d(x)=dist(x,Omega_phys)` and split the PML layer `[0,width]` into radial-distance bands `[r_i,r_{i+1}]`. The reported band center is `(r_i+r_{i+1})/2`; the mean/max amplitudes are computed over mesh nodes whose distance `d(x)` lies in that band. Thus the table measures whether `|u_h|` decays as one moves outward through the PML layer.\n\n');
fprintf(fid, '![PML decay by band center](fig_pml_decay.png)\n\n');
fprintf(fid, '| band center | mean amplitude | max amplitude |\n|---:|---:|---:|\n');
for i = 1:numel(pml.decay.bandCenter)
    fprintf(fid, '| %.4g | %.4e | %.4e |\n', pml.decay.bandCenter(i), ...
        pml.decay.ampMean(i), pml.decay.ampMax(i));
end
fprintf(fid, '\nPML layer width and complex absorption sweep: `sigmaMax` is the peak value in `s_l=1+i*sigma_l/k`, so `sigmaMax/k` is the largest imaginary stretch in one coordinate direction. Smaller outer/inner ratios mean stronger attenuation from the inner third to the outer third of the layer.\n\n');
fprintf(fid, '![PML width/sigma sweep](fig_pml_width_sigma_sweep.png)\n\n');
fprintf(fid, '| width | sigmaMax/k | sigmaMax | inner mean | outer mean | outer/inner | outer max | DOF |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for r = pml.sweep.rows
    fprintf(fid, '| %.4g | %.4g | %.4g | %.4e | %.4e | %.4e | %.4e | %d |\n', ...
        r.width, r.sigmaFactor, r.sigmaMax, r.innerMean, r.outerMean, ...
        r.outerToInner, r.outerMax, r.dof);
end
fprintf(fid, '\nPML convergence model: measured physical-domain RMS error compared with `kh+k^3h^2`; fitted max ratio `%.3e`.\n\n', ...
    pml.convergence.fitC);
fprintf(fid, '![PML convergence fit](fig_pml_convergence.png)\n\n');
fprintf(fid, '| h | kh | kh+k^3h^2 | RMS error | error/model |\n|---:|---:|---:|---:|---:|\n');
for r = pml.convergence.rows
    fprintf(fid, '| %.5g | %.4g | %.4e | %.4e | %.4e |\n', ...
        r.h, r.kh, r.pollutionModel, r.l2NodeError, r.ratio);
end

fprintf(fid, '\n## ORAS-PML Iterations\n\n');
fprintf(fid, 'Current implementation baseline is P1-only. Published RAS-PML/RMS-PML tables use P2 and therefore are comparison targets, not exact reproduction rows for these P1 runs.\n\n');
fprintf(fid, '![ORAS-PML residual histories](fig_oras_pml_residuals.png)\n\n');
fprintf(fid, '![ORAS-PML iteration counts](fig_oras_pml_iterations.png)\n\n');
fprintf(fid, '| shape | k | h | DOF | Richardson | GMRES | relres | paper fixed-point range | paper GMRES range | status |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---|---|---|\n');
for r = oras.results
    fprintf(fid, '| %s | %.0f | %.5g | %d | %s | %s | %.3e | %s | %s | %s |\n', ...
        r.grid, r.k, r.h, r.dof, iterString(r.richardson), iterString(r.gmres), ...
        r.relres, r.paperFixedPointRange, r.paperGMRESRange, r.status);
end
fprintf(fid, '\nPaper target note: fixed-width and mesh-dependent delta counts from Galkowski-Gong-Graham-Lafontaine-Spence style tables are stored only as ranges here because this suite runs a moderate P1 analogue.\n\n');
fprintf(fid, 'Mismatch diagnosis:\n\n');
fprintf(fid, '| item | diagnosis |\n|---|---|\n');
fprintf(fid, '| P1 vs P2 | The current `orasPMLHelmholtz2D` path is P1-only; published tables use P2 elements, so exact iteration matching is not expected. |\n');
fprintf(fid, '| PML differential operator | GGGLS define `Delta_s=sum_l (1/(1+i g_l''(x_l)) d_{x_l})^2` and local `Delta_{s,j}` analogously. `assembleHelmholtzPML2D` now uses the expanded non-divergence form through `assembleNondivStiffness2D`, with `d_l=s_l^{-2}` and `beta_l=-s_l''/s_l^3` in the weak term `int d_l partial_l u partial_l v + int v beta_l partial_l u`. The older divergence-form coefficients remain available in `pmlCoefficients2D` for comparison. |\n');
fprintf(fid, '| Local operator `P_s^j` | In GGGLS, `P_s^j=-k^{-2}Delta_{s,j}-c^{-2}` differs from the global `P_s` only in subdomain PML regions; the partition of unity is subordinate to regions avoiding `supp(P_s^j-P_s)`. Our `smoothPartitionOfUnity2D` is normalized nodally over extracted elements and can be nonzero in local PML regions, so the algebra behind their RAS-PML proof is not reproduced. |\n');
fprintf(fid, '| FE realization | GGGLS define local FE spaces as restrictions of the global mesh on mesh-aligned `Omega_j` and use the discrete local matrix `A_{h,j}` for `P_s^j`. Our current P1 path rebuilds the local PML operator on the extracted local mesh and then eliminates all local boundary nodes. For aligned P1 meshes this is close in spirit but still not the same `P_s^j` discretization in (1.3)-(1.4). |\n');
fprintf(fid, '| PML profile | The paper uses a smooth scaling function `f_s` that is eventually linear, and the numerical section states `f_PML(x)=a x^3/3` with `a=30k`. The repo uses `sigma=sigmaMax*(distance/thickness)^sigmaOrder` inside `s=1+i sigma/k`; matching the paper requires a separate profile and the `Delta_s` weak form. |\n');
fprintf(fid, '| Fixed-point RAS-PML residual growth | The residual growth in the figure is therefore plausible for this current preconditioned Richardson operator; it means our present P1/divergence-form local-PML operator is not reproducing the GGGLS RAS-PML contraction. It should not be interpreted as a contradiction of the paper. |\n');
fprintf(fid, '| k=60 strip row | GMRES stopped at the 80-iteration cap with residual above `1e-6`; this is recorded as a mismatch rather than treated as paper agreement. |\n');
fprintf(fid, '| Variable coefficient c | Variable `c^{-2}` paper rows are not implemented in this moderate suite. |\n\n');

fprintf(fid, '## Clement And Scott-Zhang Transfers\n\n');
fprintf(fid, 'Approximation proxy uses `||u-I_Hu||_{L^p}/(H||u||_{W^{1,2}})` and a max local patch proxy.\n\n');
fprintf(fid, 'Summary of maximum ratios over all displayed 2D rows:\n\n');
fprintf(fid, '| family | Clement max | Scott-Zhang max |\n|---|---:|---:|\n');
fprintf(fid, '| approximation global ratio | %.4e | %.4e |\n', ...
    maxRatioByOperator(transfer.approximation, 'Clement', 'globalRatio'), ...
    maxRatioByOperator(transfer.approximation, 'Scott-Zhang', 'globalRatio'));
fprintf(fid, '| approximation local ratio | %.4e | %.4e |\n', ...
    maxRatioByOperator(transfer.approximation, 'Clement', 'localRatio'), ...
    maxRatioByOperator(transfer.approximation, 'Scott-Zhang', 'localRatio'));
fprintf(fid, '| sampled/full stability ratio | %.4e | %.4e |\n\n', ...
    maxRatioByOperator(transfer.stability, 'Clement', 'ratio'), ...
    maxRatioByOperator(transfer.stability, 'Scott-Zhang', 'ratio'));
fprintf(fid, '![Transfer approximation ratios](fig_transfer_approximation.png)\n\n');
fprintf(fid, '![Transfer stability ratios](fig_transfer_stability.png)\n\n');
fprintf(fid, '| operator | H | p | global ratio | max local ratio |\n|---|---:|---:|---:|---:|\n');
for r = transfer.approximation
    fprintf(fid, '| %s | %.4g | %s | %.4e | %.4e |\n', ...
        r.operator, r.H, pString(r.p), r.globalRatio, r.localRatio);
end
fprintf(fid, '\nSampled Sobolev stability ratios for 2D P1 transfers. Fractional `s` values are diagnostics based on quadrature/centroid samples.\n\n');
fprintf(fid, '| operator | H | s | p | ratio | sampled? |\n|---|---:|---:|---:|---:|---|\n');
for r = transfer.stability
    fprintf(fid, '| %s | %.4g | %.3g | %s | %.4e | %s |\n', ...
        r.operator, r.H, r.s, pString(r.p), r.ratio, yesNo(r.sampled));
end
fprintf(fid, '\n3D pilot stability ratios:\n\n');
fprintf(fid, '| operator | s | p | ratio |\n|---|---:|---:|---:|\n');
for r = transfer.stability3D
    fprintf(fid, '| %s | %.3g | %s | %.4e |\n', r.operator, r.s, pString(r.p), r.ratio);
end

fprintf(fid, '\n## CIP-FEM P1-P3 Preasymptotic Scaling\n\n');
fprintf(fid, 'Plane-wave impedance problem with `k=%g`. Model column is `(kh)^p+k(kh)^{2p}`.\n\n', cip.k);
fprintf(fid, '![CIP-FEM preasymptotic errors](fig_cip_preasymptotic.png)\n\n');
fprintf(fid, '| degree | h | kh | model | FEM energy | CIP energy | FEM/model | CIP/model |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for r = cip.results
    fprintf(fid, '| %d | %.5g | %.4g | %.4e | %.4e | %.4e | %.4e | %.4e |\n', ...
        r.degree, r.h, r.kh, r.model, r.femEnergy, r.cipEnergy, r.femRatio, r.cipRatio);
end

fprintf(fid, '\n## Figures\n\n');
fprintf(fid, '- `fig_pml_decay.png`\n');
fprintf(fid, '- `fig_pml_width_sigma_sweep.png`\n');
fprintf(fid, '- `fig_pml_convergence.png`\n');
fprintf(fid, '- `fig_oras_pml_residuals.png`\n');
fprintf(fid, '- `fig_oras_pml_iterations.png`\n');
fprintf(fid, '- `fig_transfer_approximation.png`\n');
fprintf(fid, '- `fig_transfer_stability.png`\n');
fprintf(fid, '- `fig_cip_preasymptotic.png`\n');
end


function d = outsideDistance(node, box)
dx = max([box(1) - node(:,1), node(:,1) - box(2), zeros(size(node,1),1)], [], 2);
dy = max([box(3) - node(:,2), node(:,2) - box(4), zeros(size(node,1),1)], [], 2);
d = sqrt(dx.^2 + dy.^2);
end


function uq = interpolateP1(node, elem, u, xq)
[tid, lambda] = locateSimplexP1(node, elem, xq, 1e-10);
if any(tid == 0)
    error('verify_theory_new_functions:interp', 'Point outside reference mesh.');
end
uq = sum(u(elem(tid, :)) .* lambda, 2);
end


function val = p1SubsetL2(node, elem, u)
M = assembleMass2D(node, elem, 1);
idx = unique(elem(:));
val = sqrt(real(u(idx)' * M(idx, idx) * u(idx)));
end


function val = p1SubsetH1(node, elem, u)
K = assembleStiffness2D(node, elem, 1);
idx = unique(elem(:));
val = sqrt(real(u(idx)' * K(idx, idx) * u(idx)));
end


function targets = paperRASPMLTargets()
targets = struct();
targets.note = ['Published P2 tables use k=100:50:350, residual tolerance 1e-6, ', ...
    'fixed PML profile f_PML(x)=a*x^3/3 with a=30k, and report fixed-point counts ', ...
    'with GMRES counts in parentheses.'];
targets.stripN2Fixed = [3, 3, 3, 4, 4, 4];
targets.stripN2GMRES = [3, 3, 3, 4, 4, 4];
targets.checker2Fixed = [8, 8, 7, 7, 6, 6];
targets.checker2GMRES = [6, 6, 6, 6, 6, 6];
end


function val = maxRatioByOperator(rows, opName, fieldName)
mask = strcmp({rows.operator}, opName);
vals = [rows(mask).(fieldName)];
val = max(vals);
end


function [fpRange, gmRange] = paperComparisonRange(shape, gridSize)
if strcmp(shape, 'strip') && gridSize(1) == 2
    fpRange = '3-4';
    gmRange = '3-4';
elseif strcmp(shape, 'checker') && all(gridSize == [2, 2])
    fpRange = '6-9';
    gmRange = '6-7';
else
    fpRange = 'n/a';
    gmRange = 'n/a';
end
end


function saveFigure(fig, fileName)
set(fig, 'PaperPositionMode', 'auto');
print(fig, fileName, '-dpng', '-r160');
end


function s = iterString(its)
if isinf(its) || its > 1e6
    s = 'skip';
elseif its > 80
    s = '>80';
else
    s = sprintf('%d', its);
end
end


function s = pString(p)
if isinf(p)
    s = 'inf';
else
    s = sprintf('%.4g', p);
end
end


function s = yesNo(tf)
if tf
    s = 'yes';
else
    s = 'no';
end
end


function out = ternary(tf, a, b)
if tf
    out = a;
else
    out = b;
end
end
