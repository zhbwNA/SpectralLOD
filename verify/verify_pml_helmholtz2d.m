% VERIFY_PML_HELMHOLTZ2D  Smoke checks for P1 Helmholtz PML assembly/preconditioner.

fprintf('========== PML Helmholtz 2D Verification ==========\n\n');

phys = [0, 1, 0, 1];
pbox = [-0.25, 1.25, -0.25, 1.25];
[node, elem, bdFlag] = squaremesh(pbox, 0.25); %#ok<ASGLU>
k = 8;
pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
    'sigmaMax', 2*k, 'sigmaOrder', 2, 'quadOrder', 4);

fprintf('Test 1: PML coefficients ... ');
[a11, a22, bcoef] = pmlCoefficients2D(0.5, 0.5, k, pml);
assert(abs(a11 - 1) < 1e-14 && abs(a22 - 1) < 1e-14 && abs(bcoef - 1) < 1e-14, ...
    'PML coefficients must reduce to physical coefficients inside Omega.');
[a11p, a22p, bcoefp] = pmlCoefficients2D(-0.2, 0.5, k, pml);
assert(abs(imag(a11p)) > 0 || abs(imag(a22p)) > 0 || abs(imag(bcoefp)) > 0, ...
    'PML coefficients must be complex in the layer.');
fprintf('PASSED\n');

fprintf('Test 2: zero-PML assembly equals K-k^2M ... ');
pml0 = struct('physicalBox', pbox, 'pmlBox', pbox, 'sigmaMax', 0, 'quadOrder', 2);
[A0, ~, freeDof, bdDof] = assembleHelmholtzPML2D(node, elem, k, pml0, 0);
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
rel = norm(A0 - (K - k^2*M), 'fro') / max(1, norm(K - k^2*M, 'fro'));
assert(rel < 1e-12, 'Zero-PML matrix mismatch: %.3e', rel);
assert(~isempty(freeDof) && ~isempty(bdDof), 'Dirichlet dof split must be nonempty.');
fprintf('PASSED  (rel %.2e)\n', rel);

fprintf('Test 3: PML local-solver preconditioner apply ... ');
parts = partitionMesh2D(node, elem, bdFlag, 2, 'overlap', 0.25);
parts = smoothPartitionOfUnity2D(parts, pbox, [2, 1], 0.25);
applyB = orasPMLHelmholtz2D(node, elem, k, parts, pml, 'lu', false);
rng(1);
r = randn(size(node,1), 1) + 1i*randn(size(node,1), 1);
z = applyB(r);
assert(all(isfinite(real(z))) && all(isfinite(imag(z))), ...
    'PML preconditioner returned non-finite values.');
assert(norm(z) > 0, 'PML preconditioner returned the zero vector for nonzero input.');
fprintf('PASSED\n');

fprintf('\n========== PML Helmholtz 2D tests PASSED ==========\n');
