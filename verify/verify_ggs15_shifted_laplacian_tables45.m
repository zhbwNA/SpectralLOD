% VERIFY_GGS15_SHIFTED_LAPLACIAN_TABLES45  Reproduce GGS15 Tables 4-5 GMRES trends.

fprintf('========== GGS15 shifted-Laplacian Tables 4-5 reproduction ==========\n');

cfg = ggs15Config();
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'verify', 'ggs15_shifted_laplacian');
docDir = fullfile(repoRoot, 'docs', 'GGS15_shifted_laplacian');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(docDir, 'dir'), mkdir(docDir); end

paper = ggs15PaperTables();
allResults = struct([]);

for tableId = 4:5
    if tableId == 4
        kvals = cfg.table4K;
    else
        kvals = cfg.table5K;
    end

    fprintf('Table %d configured k values: %s\n', tableId, mat2str(kvals));
    for ik = 1:numel(kvals)
        k = kvals(ik);
        n = ggs15MeshN(tableId, k);
        ndofEstimate = (n + 1)^2;
        if ndofEstimate > cfg.maxDof
            fprintf('  skip Table %d, k=%g: estimated dofs %d exceed GGS15_MAXDOF=%d\n', ...
                tableId, k, ndofEstimate, cfg.maxDof);
            allResults = ggs15AppendResult(allResults, ...
                ggs15SkippedResult(tableId, k, n, ndofEstimate, ...
                sprintf('skipped: dofs %d exceed maxDof %d', ndofEstimate, cfg.maxDof)));
            continue;
        end

        fprintf('  assembling Table %d, k=%g, n=%d, h=%.4g ...\n', tableId, k, n, 1/n);
        [node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/n);
        [A, ~] = assembleHelmholtz2D(node, elem, bdFlag, k, 0, 0, 1);
        b = ones(size(A, 1), 1);

        for is = 1:numel(paper.shiftLabels)
            label = paper.shiftLabels{is};
            epsilon = ggs15ShiftValue(label, k);
            target = ggs15Target(paper, tableId, k, is);
            fprintf('    shift %-7s epsilon=%.6g ... ', label, epsilon);
            try
                opts = struct('epsilon', epsilon, 'eta', 'sqrt', 'solverMode', 'lu');
                [applyPrecon, Aeps] = shiftedLaplacianPreconditioner2D(node, elem, bdFlag, k, 1, opts);
                rhs = applyPrecon(b);
                Afun = @(x) applyPrecon(A * x);
                maxIt = min(cfg.maxIt, size(A, 1));
                [~, flag, relres, iter, resvec] = gmres(Afun, rhs, [], cfg.tol, maxIt);
                iterCount = ggs15GmresIters(iter);
                dApprox = NaN;
                if cfg.computeD && size(A, 1) <= cfg.dMaxDof
                    B = ggs15DensePreconditionedMatrix(Aeps, A);
                    dApprox = ggs15NumericalRangeDistance(B, cfg.dAngles);
                end
                fprintf('it=%d flag=%d relres=%.2e\n', iterCount, flag, relres);
                allResults = ggs15AppendResult(allResults, ...
                    ggs15Result(tableId, k, n, size(A,1), label, epsilon, ...
                    target.d, target.it, dApprox, iterCount, flag, relres, numel(resvec)-1, 'ran'));
            catch ME
                fprintf('failed: %s\n', ME.message);
                allResults = ggs15AppendResult(allResults, ...
                    ggs15Result(tableId, k, n, size(A,1), label, epsilon, ...
                    target.d, target.it, NaN, NaN, NaN, NaN, NaN, ME.message));
            end
        end
    end
end

save(fullfile(outDir, 'ggs15_shifted_laplacian_tables45_results.mat'), ...
    'cfg', 'paper', 'allResults');
ggs15WriteReport(fullfile(docDir, 'tables_4_5_results.md'), cfg, paper, allResults);

fprintf('Report written to %s\n', fullfile(docDir, 'tables_4_5_results.md'));
fprintf('========== GGS15 shifted-Laplacian reproduction complete ==========\n');


function cfg = ggs15Config()
cfg.table4K = ggs15EnvNumeric('GGS15_TABLE4_KVALS', [10 20]);
cfg.table5K = ggs15EnvNumeric('GGS15_TABLE5_KVALS', 10);
cfg.maxDof = ggs15EnvScalar('GGS15_MAXDOF', 20000);
cfg.tol = ggs15EnvScalar('GGS15_GMRES_TOL', 1e-6);
cfg.maxIt = ggs15EnvScalar('GGS15_MAXIT', 400);
cfg.computeD = ggs15EnvLogical('GGS15_COMPUTE_D', true);
cfg.dMaxDof = ggs15EnvScalar('GGS15_D_MAXDOF', 700);
cfg.dAngles = ggs15EnvScalar('GGS15_D_ANGLES', 48);
cfg.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
end


function values = ggs15EnvNumeric(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValue;
else
    values = str2num(raw); %#ok<ST2NM>
    if isempty(values)
        error('Environment variable %s must be a numeric vector.', name);
    end
end
end


function value = ggs15EnvScalar(name, defaultValue)
values = ggs15EnvNumeric(name, defaultValue);
value = values(1);
end


function value = ggs15EnvLogical(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    value = defaultValue;
else
    value = any(strcmp(raw, {'1', 'true', 'yes', 'on'}));
end
end


function paper = ggs15PaperTables()
paper.source = 'Gander, Graham, and Spence, Numerische Mathematik 131 (2015), Tables 4-5';
paper.url = 'https://purehost.bath.ac.uk/ws/portalfiles/portal/109280714/GaGrSp14_revised.pdf';
paper.shiftLabels = {'k/4', 'k/2', 'k', '2k', '4k', 'k^(3/2)', 'k^2'};
paper.table4.k = [10 20 40 80 160];
paper.table4.d = [
    .9328 .8714 .7641 .5971 .3861 .4594 .1466;
    .9272 .8618 .7493 .5797 .3729 .3413 .0538;
    .9246 .8569 .7411 .5675 .3590 .2311 .0156;
    .9230 .8540 .7360 .5610 .3525 .1477 .0039;
    .9223 .8525 .7336 .5547 .3439 .0870 .0030];
paper.table4.it = [
    4 5 6 7 9 8 13;
    4 5 6 8 10 11 25;
    4 5 6 8 11 13 47;
    4 5 6 7 10 16 84;
    4 5 6 7 10 19 148];
paper.table5.k = [10 20 40 80];
paper.table5.d = [
    .9323 .8706 .7627 .5943 .3812 .4550 .1432;
    .9260 .8595 .7458 .5749 .3704 .3367 .0525;
    .9226 .8535 .7358 .5609 .3529 .2275 .0150;
    .9201 .8490 .7283 .5504 .3417 .1443 .0056];
paper.table5.it = [
    4 5 6 7 9 8 13;
    4 5 6 8 11 11 24;
    4 5 6 8 11 14 48;
    4 5 6 8 10 16 86];
end


function n = ggs15MeshN(tableId, k)
if tableId == 4
    n = 2 * k;
else
    n = ceil(k^(3/2));
end
end


function epsilon = ggs15ShiftValue(label, k)
switch label
    case 'k/4'
        epsilon = k / 4;
    case 'k/2'
        epsilon = k / 2;
    case 'k'
        epsilon = k;
    case '2k'
        epsilon = 2 * k;
    case '4k'
        epsilon = 4 * k;
    case 'k^(3/2)'
        epsilon = k^(3/2);
    case 'k^2'
        epsilon = k^2;
    otherwise
        error('Unknown shift label %s.', label);
end
end


function target = ggs15Target(paper, tableId, k, shiftIdx)
if tableId == 4
    rows = paper.table4.k;
    dTable = paper.table4.d;
    itTable = paper.table4.it;
else
    rows = paper.table5.k;
    dTable = paper.table5.d;
    itTable = paper.table5.it;
end
row = find(rows == k, 1);
target = struct('d', NaN, 'it', NaN);
if ~isempty(row)
    target.d = dTable(row, shiftIdx);
    target.it = itTable(row, shiftIdx);
end
end


function r = ggs15Result(tableId, k, n, ndof, label, epsilon, paperD, paperIt, repoD, repoIt, flag, relres, resSteps, note)
r = struct('table', tableId, 'k', k, 'n', n, 'ndof', ndof, 'shift', label, ...
    'epsilon', epsilon, 'paperD', paperD, 'paperIt', paperIt, 'repoD', repoD, ...
    'repoIt', repoIt, 'flag', flag, 'relres', relres, 'resSteps', resSteps, 'note', note);
end


function r = ggs15SkippedResult(tableId, k, n, ndof, note)
r = ggs15Result(tableId, k, n, ndof, '', NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, note);
end


function results = ggs15AppendResult(results, result)
if isempty(results)
    results = result;
else
    results(end+1) = result;
end
end


function it = ggs15GmresIters(iter)
if numel(iter) == 2
    it = iter(2);
else
    it = iter(end);
end
end


function B = ggs15DensePreconditionedMatrix(Aeps, A)
B = Aeps \ full(A);
end


function d = ggs15NumericalRangeDistance(B, nAngles)
theta = linspace(0, 2*pi, nAngles + 1);
theta(end) = [];
z = zeros(2 * nAngles, 1);
idx = 0;
for itheta = 1:nAngles
    w = exp(-1i * theta(itheta));
    H = (w * B + conj(w) * B') / 2;
    [V, D] = eig(full(H), 'vector');
    [~, imax] = max(real(D));
    [~, imin] = min(real(D));
    idx = idx + 1;
    v = V(:, imax);
    z(idx) = v' * B * v;
    idx = idx + 1;
    v = V(:, imin);
    z(idx) = v' * B * v;
end
d = min(abs(z));
end


function ggs15WriteReport(path, cfg, paper, results)
fid = fopen(path, 'w');
if fid < 0
    error('Could not open report path %s.', path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Reproduction target: Gander-Graham-Spence shifted-Laplacian GMRES paper, Tables 4-5.\n');
fprintf(fid, 'Created: 2026-05-25\n');
fprintf(fid, 'Updated: 2026-05-25\n');
fprintf(fid, 'Verification entry point: `verify/verify_ggs15_shifted_laplacian_tables45.m`\n');
fprintf(fid, 'Main utilities: `squaremesh`, `assembleHelmholtz2D`, `shiftedLaplacianPDE`, `shiftedLaplacianPreconditioner2D`, MATLAB `gmres`\n\n');

fprintf(fid, '# GGS15 Shifted-Laplacian GMRES Tables 4-5\n\n');
fprintf(fid, 'Source: %s. Author manuscript: <%s>.\n\n', paper.source, paper.url);
fprintf(fid, 'This report records a repo-local reproduction path for Tables 4 and 5 of the shifted-Laplacian GMRES paper. ');
fprintf(fid, 'The paper reports the distance `d` from the origin to the numerical range of the left-preconditioned matrix and unrestarted GMRES iterations to relative residual `1e-6`.\n\n');

fprintf(fid, '## Mathematical Form\n\n');
fprintf(fid, 'The model problem is the 2D interior impedance Helmholtz problem on `Omega=(0,1)^2`. ');
fprintf(fid, 'The repo discretization uses continuous P1 finite elements on the uniform triangulation returned by `squaremesh([0,1,0,1],1/n)`. ');
fprintf(fid, 'The unshifted matrix is\n\n');
fprintf(fid, '```text\nA = K - k^2 M - i k M_boundary.\n```\n\n');
fprintf(fid, 'For a shift `epsilon > 0`, the shifted-Laplacian preconditioner matrix is\n\n');
fprintf(fid, '```text\nA_epsilon = K - (k^2 + i epsilon) M - i sqrt(k^2 + i epsilon) M_boundary.\n```\n\n');
fprintf(fid, 'The verified linear system is the left-preconditioned system\n\n');
fprintf(fid, '```text\nA_epsilon^{-1} A u = A_epsilon^{-1} b,\n```\n\n');
fprintf(fid, 'with `b = ones(N,1)`, zero initial guess, and unrestarted MATLAB `gmres` tolerance `%.1e`. ', cfg.tol);
fprintf(fid, 'This follows the paper table setup except for the numerical-range distance, which is computed here only as a dense angular approximation when `N <= GGS15_D_MAXDOF`.\n\n');

fprintf(fid, '## Paper Parameters\n\n');
fprintf(fid, '- Table 4 mesh rule: `n = 2*k`.\n');
fprintf(fid, '- Table 5 mesh rule: `n = ceil(k^(3/2))`.\n');
fprintf(fid, '- Shift columns: `epsilon = k/4, k/2, k, 2k, 4k, k^(3/2), k^2`.\n');
fprintf(fid, '- Default configured Table 4 rows: `%s`.\n', mat2str(cfg.table4K));
fprintf(fid, '- Default configured Table 5 rows: `%s`.\n', mat2str(cfg.table5K));
fprintf(fid, '- Default maximum DOFs: `%d`; dense `d` approximation maximum DOFs: `%d`.\n\n', cfg.maxDof, cfg.dMaxDof);

ggs15WritePaperTable(fid, paper, 4);
ggs15WritePaperTable(fid, paper, 5);

fprintf(fid, '## Repo Run Results\n\n');
if isempty(results)
    fprintf(fid, 'No rows were run.\n\n');
else
    for tableId = 4:5
        rows = results([results.table] == tableId);
        fprintf(fid, '### Table %d Configured Run\n\n', tableId);
        if isempty(rows)
            fprintf(fid, 'No rows were configured for this table.\n\n');
        else
            fprintf(fid, '| k | n | DOFs | epsilon | paper d | repo d | paper it | repo it | flag | relres | note |\n');
            fprintf(fid, '|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|:---|\n');
            for ir = 1:numel(rows)
                r = rows(ir);
                fprintf(fid, '| %g | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
                    r.k, r.n, r.ndof, ggs15Md(r.shift), ggs15Fmt(r.paperD, '%.3f'), ...
                    ggs15Fmt(r.repoD, '%.3f'), ggs15Fmt(r.paperIt, '%.0f'), ...
                    ggs15Fmt(r.repoIt, '%.0f'), ggs15Fmt(r.flag, '%.0f'), ...
                    ggs15Fmt(r.relres, '%.2e'), ggs15Md(r.note));
            end
            fprintf(fid, '\n');
        end
    end
end

fprintf(fid, '## Re-run Controls\n\n');
fprintf(fid, 'PowerShell examples:\n\n');
fprintf(fid, '```powershell\n');
fprintf(fid, '$env:GGS15_TABLE4_KVALS = ''10 20 40''\n');
fprintf(fid, '$env:GGS15_TABLE5_KVALS = ''10 20''\n');
fprintf(fid, '$env:GGS15_MAXDOF = ''50000''\n');
fprintf(fid, 'matlab -nosplash -nodesktop -batch "addpath(genpath(''.'')); run(''verify/verify_ggs15_shifted_laplacian_tables45.m'');"\n');
fprintf(fid, '```\n\n');
fprintf(fid, 'Set `GGS15_COMPUTE_D=0` to skip dense numerical-range estimates. ');
fprintf(fid, 'Set `GGS15_D_MAXDOF` and `GGS15_D_ANGLES` higher only for small matrices, because the approximation forms the dense preconditioned matrix and diagonalizes angular Hermitian parts.\n\n');

fprintf(fid, '## Status\n\n');
fprintf(fid, 'The verifier now gives a self-contained reproduction harness for Tables 4-5 and stores the paper target values in code and in this report. ');
fprintf(fid, 'The default run is a small executable subset of the paper tables; larger rows should be run explicitly after choosing the resource envelope.\n');
end


function ggs15WritePaperTable(fid, paper, tableId)
if tableId == 4
    kRows = paper.table4.k;
    dVals = paper.table4.d;
    itVals = paper.table4.it;
    label = 'Table 4 Paper Targets: fixed points per wavelength';
else
    kRows = paper.table5.k;
    dVals = paper.table5.d;
    itVals = paper.table5.it;
    label = 'Table 5 Paper Targets: fixed scaled points per wavelength';
end
fprintf(fid, '### %s\n\n', label);
fprintf(fid, '| k | metric |');
for is = 1:numel(paper.shiftLabels)
    fprintf(fid, ' %s |', paper.shiftLabels{is});
end
fprintf(fid, '\n|---:|:---|');
for is = 1:numel(paper.shiftLabels)
    fprintf(fid, '---:|');
end
fprintf(fid, '\n');
for ik = 1:numel(kRows)
    fprintf(fid, '| %g | d |', kRows(ik));
    for is = 1:numel(paper.shiftLabels)
        fprintf(fid, ' %.3f |', dVals(ik, is));
    end
    fprintf(fid, '\n| %g | it |', kRows(ik));
    for is = 1:numel(paper.shiftLabels)
        fprintf(fid, ' %d |', itVals(ik, is));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n');
end


function s = ggs15Fmt(x, fmt)
if isempty(x) || (isnumeric(x) && isnan(x))
    s = '--';
else
    s = sprintf(fmt, x);
end
end


function s = ggs15Md(x)
if isempty(x)
    s = '--';
else
    s = char(string(x));
    s = strrep(s, '|', '\|');
end
end
