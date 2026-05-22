% VERIFY_VARIABLE_K_SHIFTED_LAPLACIAN  Check variable k and shifted wrappers.

fprintf('========== Variable k and shifted-Laplacian verification ==========\n');

[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/4);
k0 = 7;

fprintf('Test 1: Helmholtz PDE struct preserves constant-k assembly ... ');
[Aold, bold] = assembleHelmholtz2D(node, elem, bdFlag, k0, 1, 0, 1);
pde0 = helmholtzPDE(k0, 'source', 1, 'boundaryData', 0);
[Apde, bpde] = assembleHelmholtz2D(node, elem, bdFlag, pde0, [], [], 1);
relA = norm(Aold - Apde, 'fro') / max(1, norm(Aold, 'fro'));
relb = norm(bold - bpde) / max(1, norm(bold));
assert(relA < 1e-12 && relb < 1e-12, 'constant-k PDE struct changed Helmholtz assembly.');
fprintf('passed (relA %.2e)\n', relA);

fprintf('Test 2: shifted Laplacian constant formula ... ');
eps0 = 0.5 * k0^2;
eta0 = sqrt(k0^2 + 1i * eps0);
pdeShift = shiftedLaplacianPDE(k0, 'epsilon', eps0, 'eta', 'sqrt');
[Ashift, ~] = assembleHelmholtz2D(node, elem, bdFlag, pdeShift, 0, 0, 1);
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
Mb = assembleBoundaryMass2D(node, elem, bdFlag, 1);
Aexpected = K - (k0^2 + 1i * eps0) * M - 1i * eta0 * Mb;
relShift = norm(Ashift - Aexpected, 'fro') / max(1, norm(Aexpected, 'fro'));
assert(relShift < 1e-12, 'shifted Laplacian formula mismatch.');
fprintf('passed (relA %.2e)\n', relShift);

fprintf('Test 3: variable k(x,y), epsilon(k), eta(k) assembly ... ');
kfun = @(x,y) 6 + x + 0.5 * y;
epsfun = @(kv) 0.25 * kv.^2;
etafun = @(kv) sqrt(kv.^2 + 1i * 0.25 * kv.^2);
pdeVar = shiftedLaplacianPDE(kfun, 'epsilon', epsfun, 'eta', etafun);
[Avar, bvar] = assembleHelmholtz2D(node, elem, bdFlag, pdeVar, @(x,y) x + y, 0, 1);
assert(all(isfinite(nonzeros(Avar))) && all(isfinite(bvar)), 'variable-k Helmholtz assembly produced non-finite entries.');
assert(norm(Avar - Aold, 'fro') > 1e-8, 'variable-k matrix unexpectedly equals constant-k matrix.');
fprintf('passed (nnz %d)\n', nnz(Avar));

fprintf('Test 4: shifted-Laplacian preconditioner wrapper ... ');
[applyPrecon, Aeps] = shiftedLaplacianPreconditioner2D(node, elem, bdFlag, k0, 1, ...
    struct('epsilon', eps0, 'eta', 'sqrt', 'solverMode', 'lu'));
r = ones(size(Aeps, 1), 1);
z = applyPrecon(r);
relSolve = norm(Aeps * z - r) / norm(r);
assert(relSolve < 1e-10, 'shifted-Laplacian LU wrapper is inaccurate.');
fprintf('passed (relSolve %.2e)\n', relSolve);

fprintf('Test 5: shifted-Laplacian preconditioner accepts variable k(x,y) ... ');
[applyVarPrecon, AepsVar, pdeShiftVar] = shiftedLaplacianPreconditioner2D(node, elem, bdFlag, kfun, 1, ...
    struct('epsilon', epsfun, 'eta', etafun, 'solverMode', 'lu'));
rVar = (1:size(AepsVar, 1)).';
zVar = applyVarPrecon(rVar);
relSolveVar = norm(AepsVar * zVar - rVar) / norm(rVar);
assert(relSolveVar < 1e-10, 'variable-k shifted-Laplacian LU wrapper is inaccurate.');
assert(isa(pdeShiftVar.k, 'function_handle'), 'preconditioner did not preserve function-valued k.');
assert(norm(AepsVar - Aeps, 'fro') > 1e-8, 'variable-k shifted preconditioner unexpectedly equals constant-k preconditioner.');
fprintf('passed (relSolve %.2e)\n', relSolveVar);

fprintf('Test 6: Maxwell NE_1 wrapper preserves constant-k assembly ... ');
[Amax, info] = assembleMaxwell2D(node, elem, bdFlag, k0);
C = assembleCurlCurl2D(node, elem);
Mn = assembleNedMass2D(node, elem);
Mbn = assembleNedBndMass2D(node, elem, bdFlag);
AmaxExpected = C - k0^2 * Mn - 1i * k0 * Mbn;
relMax = norm(Amax - AmaxExpected, 'fro') / max(1, norm(AmaxExpected, 'fro'));
assert(relMax < 1e-12 && isfield(info, 'pde'), 'constant-k Maxwell wrapper mismatch.');
fprintf('passed (relA %.2e)\n', relMax);

fprintf('Test 7: 3D Helmholtz PDE struct preserves constant-k assembly ... ');
[node3, elem3, bd3] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[Aold3, ~] = assembleHelmholtz3D(node3, elem3, bd3, k0, 0, 0, 1);
pde3 = helmholtzPDE(k0, 'source', 0, 'boundaryData', 0);
[Apde3, ~] = assembleHelmholtz3D(node3, elem3, bd3, pde3, [], [], 1);
rel3 = norm(Aold3 - Apde3, 'fro') / max(1, norm(Aold3, 'fro'));
assert(rel3 < 1e-12, 'constant-k 3D PDE struct changed Helmholtz assembly.');
fprintf('passed (relA %.2e)\n', rel3);

fprintf('Test 8: 3D shifted-Laplacian preconditioner wrapper ... ');
[applyPrecon3, Aeps3] = shiftedLaplacianPreconditioner3D(node3, elem3, bd3, k0, 1, ...
    struct('epsilon', eps0, 'eta', 'sqrt', 'solverMode', 'lu'));
r3 = ones(size(Aeps3, 1), 1);
z3 = applyPrecon3(r3);
relSolve3 = norm(Aeps3 * z3 - r3) / norm(r3);
assert(relSolve3 < 1e-10, '3D shifted-Laplacian LU wrapper is inaccurate.');
fprintf('passed (relSolve %.2e)\n', relSolve3);

fprintf('========== Variable k and shifted-Laplacian verification complete ==========\n');
