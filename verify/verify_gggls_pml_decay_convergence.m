% VERIFY_GGGLS_PML_DECAY_CONVERGENCE  PML decay and preasymptotic checks for GGGLS assembly.

fprintf('========== GGGLS PML Decay And Preasymptotic Verification ==========\n\n');

verifyDir = fileparts(mfilename('fullpath'));
outDir = fullfile(verifyDir, 'gggls_pml_theory');
if ~exist(outDir, 'dir'), mkdir(outDir); end

decay = ggglsPMLDecayExperiment();
sweep = ggglsPMLSweepExperiment();
conv = ggglsPMLConvergenceExperiment();

plotGGGLSPMLDecay(outDir, decay);
plotGGGLSPMLSweep(outDir, sweep);
plotGGGLSPMLConvergence(outDir, conv);
writeGGGLSPMLReport(fullfile(outDir, 'gggls_pml_decay_convergence_results.md'), decay, sweep, conv);
save(fullfile(outDir, 'gggls_pml_decay_convergence_results.mat'), 'decay', 'sweep', 'conv');

fprintf('GGGLS PML decay outer/inner amplitude ratio: %.3e\n', decay.outerToInner);
fprintf('GGGLS PML sweep best outer/inner amplitude ratio: %.3e\n', min([sweep.rows.outerToInner]));
fprintf('GGGLS PML convergence max error/(kh+k^3h^2): %.3e\n', conv.fitC);
fprintf('Results written to %s\n', outDir);
fprintf('========== GGGLS PML verification complete ==========\n');

assert(decay.outerToInner < 1, 'GGGLS PML decay check did not attenuate into the layer.');
assert(all(isfinite([conv.rows.l2NodeError])), 'GGGLS PML convergence produced non-finite errors.');


function decay = ggglsPMLDecayExperiment()
k = 20;
width = 0.20;
h = 0.04;
alpha = 100;
phys = [0, 1, 0, 1];
pbox = [-width, 1 + width, -width, 1 + width];
src = @(x,y) exp(-180 * ((x - 0.48).^2 + (y - 0.52).^2));
[node, ~, u] = solveGGGLSPMLSource(pbox, phys, h, k, alpha, src, 1);

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
decay = struct('k', k, 'h', h, 'width', width, 'alpha', alpha, ...
    'bandCenter', bandCenter(:), 'ampMean', ampMean, 'ampMax', ampMax, ...
    'outerToInner', ampMean(end) / max(ampMean(2), eps));
end


function sweep = ggglsPMLSweepExperiment()
k = 20;
h = 0.05;
phys = [0, 1, 0, 1];
widthVals = [0.10, 0.15, 0.20, 0.30];
alphaVals = [50, 100, 250, 500, 1000, 2500, 5000];
src = @(x,y) exp(-180 * ((x - 0.48).^2 + (y - 0.52).^2));
rows = struct('width', {}, 'alpha', {}, 'maxImagStretch', {}, ...
    'innerMean', {}, 'outerMean', {}, 'outerToInner', {}, 'outerMax', {}, 'dof', {});

for width = widthVals
    for alpha = alphaVals
        pbox = [-width, 1 + width, -width, 1 + width];
        [node, ~, u, freeDof] = solveGGGLSPMLSource(pbox, phys, h, k, alpha, src, 1);
        d = outsideDistance(node, phys);
        innerMask = d > 0 & d <= width / 3;
        outerMask = d > 2 * width / 3 & d <= width + 1e-12;
        if ~any(innerMask), innerMask = d <= width / 3; end
        rows(end+1) = struct('width', width, 'alpha', alpha, ...
            'maxImagStretch', alpha * width^2, ...
            'innerMean', mean(abs(u(innerMask))), ...
            'outerMean', mean(abs(u(outerMask))), ...
            'outerToInner', mean(abs(u(outerMask))) / max(mean(abs(u(innerMask))), eps), ...
            'outerMax', max(abs(u(outerMask))), ...
            'dof', numel(freeDof)); %#ok<AGROW>
    end
end
sweep = struct('k', k, 'h', h, 'widthVals', widthVals, 'alphaVals', alphaVals, 'rows', rows);
end


function conv = ggglsPMLConvergenceExperiment()
k = 8;
width = 0.25;
alpha = 100;
degree = 1;
phys = [0, 1, 0, 1];
pbox = [-width, 1 + width, -width, 1 + width];
src = @(x,y) exp(-120 * ((x - 0.45).^2 + (y - 0.55).^2));
hRef = 1 / 96;
[refNode, refElem, uRef] = solveGGGLSPMLSource(pbox, phys, hRef, k, alpha, src, degree);

hVals = [1/12, 1/16, 1/24, 1/32];
rows = struct('h', {}, 'kh', {}, 'pollutionModel', {}, 'l2NodeError', {}, 'ratio', {});
for h = hVals
    [node, ~, u] = solveGGGLSPMLSource(pbox, phys, h, k, alpha, src, degree);
    uRefAtNode = interpolateP1(refNode, refElem, uRef, node);
    mask = node(:,1) >= phys(1) - 1e-12 & node(:,1) <= phys(2) + 1e-12 & ...
           node(:,2) >= phys(3) - 1e-12 & node(:,2) <= phys(4) + 1e-12;
    err = sqrt(mean(abs(u(mask) - uRefAtNode(mask)).^2));
    model = k * h + k^3 * h^2;
    rows(end+1) = struct('h', h, 'kh', k*h, 'pollutionModel', model, ...
        'l2NodeError', err, 'ratio', err / model); %#ok<AGROW>
end
fitC = max([rows.ratio]);
conv = struct('k', k, 'width', width, 'alpha', alpha, 'degree', degree, ...
    'hRef', hRef, 'rows', rows, 'fitC', fitC);
end


function [node, elem, u, freeDof] = solveGGGLSPMLSource(pbox, phys, h, k, alpha, src, degree)
[node, elem] = squaremesh(pbox, h);
box = struct('physicalBox', phys, 'outerBox', pbox);
opts = struct('pmlAlpha', alpha, 'quadOrder', max(4, 2 * degree + 1));
[A, b, freeDof] = assembleGGGLSPML2D(node, elem, k, box, src, degree, opts);
u = zeros(size(A, 1), 1);
u(freeDof) = A(freeDof, freeDof) \ b(freeDof);
end


function d = outsideDistance(node, box)
dx = max([box(1) - node(:,1), node(:,1) - box(2), zeros(size(node,1),1)], [], 2);
dy = max([box(3) - node(:,2), node(:,2) - box(4), zeros(size(node,1),1)], [], 2);
d = sqrt(dx.^2 + dy.^2);
end


function uq = interpolateP1(node, elem, u, xq)
[tid, lambda] = locateSimplexP1(node, elem, xq, 1e-10);
if any(tid == 0)
    error('verify_gggls_pml_decay_convergence:interp', 'Point outside reference mesh.');
end
uq = sum(u(elem(tid, :)) .* lambda, 2);
end


function plotGGGLSPMLDecay(outDir, decay)
fig = figure('Name', 'GGGLS PML decay', 'Color', 'w');
semilogy(decay.bandCenter, decay.ampMean, 'o-', 'LineWidth', 1.5); hold on;
semilogy(decay.bandCenter, decay.ampMax, 's--', 'LineWidth', 1.2);
xlabel('$\mathrm{dist}(x,\Omega_{\rm int})$', 'Interpreter', 'latex');
ylabel('$|u_h|$ band amplitude', 'Interpreter', 'latex');
title('GGGLS non-divergence PML attenuation', 'Interpreter', 'latex');
legend({'mean', 'max'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_gggls_pml_decay.png'));
end


function plotGGGLSPMLSweep(outDir, sweep)
Z = nan(numel(sweep.widthVals), numel(sweep.alphaVals));
for r = sweep.rows
    iw = find(abs(sweep.widthVals - r.width) < 1e-12, 1);
    ia = find(abs(sweep.alphaVals - r.alpha) < 1e-12, 1);
    Z(iw, ia) = r.outerToInner;
end
fig = figure('Name', 'GGGLS PML width alpha sweep', 'Color', 'w');
imagesc(sweep.alphaVals, sweep.widthVals, log10(Z));
set(gca, 'YDir', 'normal');
colorbar;
xlabel('$\alpha$ in $f_{\rm PML}(t)=\alpha t^3/3$', 'Interpreter', 'latex');
ylabel('PML layer width', 'Interpreter', 'latex');
title('$\log_{10}$ outer/inner amplitude ratio', 'Interpreter', 'latex');
saveFigure(fig, fullfile(outDir, 'fig_gggls_pml_width_alpha_sweep.png'));
end


function plotGGGLSPMLConvergence(outDir, conv)
rows = conv.rows;
h = [rows.h];
err = [rows.l2NodeError];
model = conv.fitC * [rows.pollutionModel];
fig = figure('Name', 'GGGLS PML convergence', 'Color', 'w');
loglog(h, err, 'o-', h, model, 'k--', 'LineWidth', 1.5);
xlabel('$h$', 'Interpreter', 'latex');
ylabel('physical-domain RMS error', 'Interpreter', 'latex');
title('GGGLS PML check against $kh+k^3h^2$', 'Interpreter', 'latex');
legend({'measured', 'fitted $kh+k^3h^2$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on;
saveFigure(fig, fullfile(outDir, 'fig_gggls_pml_convergence.png'));
end


function writeGGGLSPMLReport(fileName, decay, sweep, conv)
fid = fopen(fileName, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# GGGLS PML Decay And Preasymptotic Verification\n\n');
fprintf(fid, 'This run repeats the earlier PML decay and P1 pre-asymptotic convergence checks using `assembleGGGLSPML2D`, i.e. the expanded non-divergence GGGLS operator with `f_PML(t)=alpha*t^3/3`.\n\n');
fprintf(fid, '## Decay By Band Center\n\n');
fprintf(fid, 'Parameters: `k=%g`, `h=%.4g`, PML width `%.4g`, `alpha=%g`. Outer/inner mean-amplitude ratio: `%.3e`.\n\n', ...
    decay.k, decay.h, decay.width, decay.alpha, decay.outerToInner);
fprintf(fid, 'Band center is `(r_i+r_{i+1})/2` after splitting the PML distance interval `[0,width]` into equal bands, where `r=dist(x,Omega_int)`.\n\n');
fprintf(fid, '![GGGLS PML decay](fig_gggls_pml_decay.png)\n\n');
fprintf(fid, '| band center | mean amplitude | max amplitude |\n|---:|---:|---:|\n');
for i = 1:numel(decay.bandCenter)
    fprintf(fid, '| %.4g | %.4e | %.4e |\n', decay.bandCenter(i), decay.ampMean(i), decay.ampMax(i));
end

fprintf(fid, '\n## Width And Absorption Sweep\n\n');
fprintf(fid, 'The complex absorbing strength is controlled by `alpha`; the maximum imaginary stretch over a layer of width `w` is `alpha*w^2`.\n\n');
fprintf(fid, '![GGGLS PML width/alpha sweep](fig_gggls_pml_width_alpha_sweep.png)\n\n');
fprintf(fid, '| width | alpha | max imag stretch | inner mean | outer mean | outer/inner | outer max | DOF |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for r = sweep.rows
    fprintf(fid, '| %.3g | %.4g | %.4g | %.4e | %.4e | %.4e | %.4e | %d |\n', ...
        r.width, r.alpha, r.maxImagStretch, r.innerMean, r.outerMean, ...
        r.outerToInner, r.outerMax, r.dof);
end

fprintf(fid, '\n## P1 Preasymptotic Convergence\n\n');
fprintf(fid, 'Parameters: `k=%g`, reference `h=%.4g`, PML width `%.4g`, `alpha=%g`. Model column is `kh+k^3h^2`; fitted max ratio is `%.3e`.\n\n', ...
    conv.k, conv.hRef, conv.width, conv.alpha, conv.fitC);
fprintf(fid, '![GGGLS PML convergence](fig_gggls_pml_convergence.png)\n\n');
fprintf(fid, '| h | kh | kh+k^3h^2 | RMS error | error/model |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|\n');
for r = conv.rows
    fprintf(fid, '| %.5g | %.5g | %.5g | %.4e | %.4e |\n', ...
        r.h, r.kh, r.pollutionModel, r.l2NodeError, r.ratio);
end
end


function saveFigure(fig, fileName)
set(fig, 'PaperPositionMode', 'auto');
print(fig, fileName, '-dpng', '-r160');
end
