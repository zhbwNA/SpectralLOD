% VERIFY_NONDIVERGENCE_ASSEMBLY  Checks for reusable non-divergence assemblers.

fprintf('========== Non-Divergence Assembly Verification ==========\n\n');

[node, elem] = squaremesh([0, 1, 0, 1], 0.25);

fprintf('Test 1: scalar gradient form reduces to stiffness ... ');
K = assembleStiffness2D(node, elem, 1);
Knd = assembleNondivStiffness2D(node, elem, 1, 1);
rel = norm(Knd - K, 'fro') / max(1, norm(K, 'fro'));
assert(rel < 1e-13, 'P1 non-divergence stiffness mismatch: %.3e', rel);
fprintf('PASSED  (rel %.2e)\n', rel);

fprintf('Test 2: first-order drift annihilates constants ... ');
coef = struct('d11', 2, 'd22', 3, 'beta1', 0.25, 'beta2', -0.5);
B = assembleNondivStiffness2D(node, elem, 1, coef);
res = norm(B * ones(size(B, 1), 1), inf);
assert(res < 1e-12, 'Gradient plus drift operator should annihilate constants: %.3e', res);
assert(norm(B - B.', 'fro') > 1e-12, 'Drift contribution should make the matrix nonsymmetric.');
fprintf('PASSED  (res %.2e)\n', res);

fprintf('Test 3: curl form reduces to NE1 curl-curl ... ');
C = assembleCurlCurl2D(node, elem);
Cnd = assembleNondivCurlCurl2D(node, elem, 1);
relCurl = norm(Cnd - C, 'fro') / max(1, norm(C, 'fro'));
assert(relCurl < 1e-13, 'NE1 non-divergence curl-curl mismatch: %.3e', relCurl);
fprintf('PASSED  (rel %.2e)\n', relCurl);

fprintf('Test 4: variable curl coefficient remains finite ... ');
Cv = assembleNondivCurlCurl2D(node, elem, @(x,y) 1 + x + 0*y);
assert(all(isfinite(nonzeros(Cv))), 'Variable curl coefficient matrix has non-finite entries.');
assert(norm(Cv - Cv.', 'fro') < 1e-12, 'Scalar weighted curl-curl matrix should remain symmetric.');
fprintf('PASSED\n');

fprintf('\n========== Non-Divergence Assembly tests PASSED ==========\n');
