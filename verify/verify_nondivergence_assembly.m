% VERIFY_NONDIVERGENCE_ASSEMBLY  Checks for reusable non-divergence assemblers.

fprintf('========== Non-Divergence Assembly Verification ==========\n\n');

[node, elem] = squaremesh([0, 1, 0, 1], 0.25);

fprintf('Test 1: scalar gradient form reduces to stiffness ... ');
K = assembleStiffness2D(node, elem, 1);
Knd = assembleDiffusion2D(node, elem, 1, 1);
rel = norm(Knd - K, 'fro') / max(1, norm(K, 'fro'));
assert(rel < 1e-13, 'P1 non-divergence stiffness mismatch: %.3e', rel);
fprintf('PASSED  (rel %.2e)\n', rel);

fprintf('Test 2: scalar split reproduces non-divergence wrapper ... ');
coef = struct('d11', 2, 'd22', 3, 'beta1', 0.25, 'beta2', -0.5);
Diff = assembleDiffusion2D(node, elem, 1, coef);
Adv = assembleAdvection2D(node, elem, 1, coef);
B = assembleNondivStiffness2D(node, elem, 1, coef);
relSplit = norm(B - (Diff + Adv), 'fro') / max(1, norm(B, 'fro'));
assert(relSplit < 1e-14, 'Scalar non-divergence split mismatch: %.3e', relSplit);
res = norm(B * ones(size(B, 1), 1), inf);
assert(res < 1e-12, 'Gradient plus drift operator should annihilate constants: %.3e', res);
assert(norm(Adv - Adv.', 'fro') > 1e-12, 'Drift contribution should make the matrix nonsymmetric.');
fprintf('PASSED  (split %.2e, res %.2e)\n', relSplit, res);

fprintf('Test 3: vector Lagrange split is componentwise ... ');
Avec = assembleVectorNondivStiffness2D(node, elem, 1, coef, 2);
AvecSplit = assembleVectorDiffusion2D(node, elem, 1, coef, 2) + ...
    assembleVectorAdvection2D(node, elem, 1, coef, 2);
relVec = norm(Avec - AvecSplit, 'fro') / max(1, norm(Avec, 'fro'));
assert(relVec < 1e-14, 'Vector Lagrange non-divergence split mismatch: %.3e', relVec);
assert(isequal(size(Avec), 2 * [size(B,1), size(B,2)]), 'Unexpected vector block matrix size.');
fprintf('PASSED  (rel %.2e)\n', relVec);

fprintf('Test 4: curl form reduces to NE1 curl-curl ... ');
C = assembleCurlCurl2D(node, elem);
Cnd = assembleNondivCurlCurl2D(node, elem, 1);
relCurl = norm(Cnd - C, 'fro') / max(1, norm(C, 'fro'));
assert(relCurl < 1e-13, 'NE1 non-divergence curl-curl mismatch: %.3e', relCurl);
fprintf('PASSED  (rel %.2e)\n', relCurl);

fprintf('Test 5: variable curl-curl dispatch matches weighted form ... ');
Cv = assembleWeightedCurlCurl2D(node, elem, @(x,y) 1 + x + 0*y);
Cdispatch = assembleCurlCurl2D(node, elem, @(x,y) 1 + x + 0*y);
relDispatch = norm(Cdispatch - Cv, 'fro') / max(1, norm(Cv, 'fro'));
assert(relDispatch < 1e-14, 'Variable curl-curl dispatch mismatch: %.3e', relDispatch);
assert(all(isfinite(nonzeros(Cv))), 'Variable curl coefficient matrix has non-finite entries.');
assert(norm(Cv - Cv.', 'fro') < 1e-12, 'Scalar weighted curl-curl matrix should remain symmetric.');
fprintf('PASSED  (rel %.2e)\n', relDispatch);

fprintf('Test 6: vector FE first order piece is reusable ... ');
Nadv0 = assembleNedAdvection2D(node, elem, [0, 0]);
assert(nnz(Nadv0) == 0, 'Zero Nedelec advection field should assemble the zero matrix.');
Nadv = assembleNedAdvection2D(node, elem, [0.25, -0.5]);
assert(all(isfinite(nonzeros(Nadv))), 'Nedelec advection matrix has non-finite entries.');
assert(norm(Nadv - Nadv.', 'fro') > 1e-12, 'Nedelec advection matrix should be nonsymmetric.');
fprintf('PASSED\n');

fprintf('\n========== Non-Divergence Assembly tests PASSED ==========\n');
