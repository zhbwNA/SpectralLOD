% VERIFY_PML_HELMHOLTZ2D  Smoke checks for P1 Helmholtz PML assembly/preconditioner.

fprintf('========== PML Helmholtz 2D Verification ==========\n\n');

phys = [0, 1, 0, 1];
pbox = [-0.25, 1.25, -0.25, 1.25];
[node, elem, bdFlag] = squaremesh(pbox, 0.25); %#ok<ASGLU>
k = 8;
pml = struct('physicalBox', phys, 'pmlBox', pbox, ...
    'sigmaMax', 2*k, 'sigmaOrder', 2, 'quadOrder', 4);

fprintf('Test 1: PML coefficients ... ');
[d11, d22, beta1, beta2] = pmlNondivCoefficients2D(0.5, 0.5, k, pml);
assert(abs(d11 - 1) < 1e-14 && abs(d22 - 1) < 1e-14 && abs(beta1) < 1e-14 && abs(beta2) < 1e-14, ...
    'PML coefficients must reduce to physical coefficients inside Omega.');
[d11p, d22p, beta1p, beta2p] = pmlNondivCoefficients2D(-0.2, 0.5, k, pml);
assert(abs(imag(d11p)) > 0 || abs(imag(d22p)) > 0 || abs(imag(beta1p)) > 0 || abs(imag(beta2p)) > 0, ...
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

fprintf('Test 3: PML assembly uses non-divergence stiffness ... ');
[Ap, ~] = assembleHelmholtzPML2D(node, elem, k, pml, 0);
coef = struct();
coef.d11 = @(x,y) pmlCoefField(x, y, k, pml, 'd11');
coef.d22 = @(x,y) pmlCoefField(x, y, k, pml, 'd22');
coef.beta1 = @(x,y) pmlCoefField(x, y, k, pml, 'beta1');
coef.beta2 = @(x,y) pmlCoefField(x, y, k, pml, 'beta2');
Knd = assembleNondivStiffness2D(node, elem, 1, coef, struct('quadOrder', pml.quadOrder));
relp = norm(Ap - (Knd - k^2*M), 'fro') / max(1, norm(Ap, 'fro'));
assert(relp < 1e-14, 'PML is not assembled through the non-divergence form: %.3e', relp);
fprintf('PASSED  (rel %.2e)\n', relp);

fprintf('Test 4: PML local-solver preconditioner apply ... ');
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


function val = pmlCoefField(x, y, k, pml, name)
[d11, d22, beta1, beta2] = pmlNondivCoefficients2D(x, y, k, pml);
switch name
    case 'd11'
        val = d11;
    case 'd22'
        val = d22;
    case 'beta1'
        val = beta1;
    case 'beta2'
        val = beta2;
end
end
