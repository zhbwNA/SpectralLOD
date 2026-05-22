% VERIFY_CIP2D  Verify 2D Lagrange CIP assembly for P1-P3.

fprintf('========== CIP 2D Verification ==========\n\n');

[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 0.5);

for degree = 1:3
    fprintf('Test degree %d: polynomial consistency ... ', degree);
    if degree == 1
        nodeH = node;
    else
        [nodeH, ~] = extendMesh2D(node, elem, degree);
    end
    C = assembleCIP2D(node, elem, degree, []);
    u = polynomialData(nodeH, degree);
    res = norm(C * u, inf);
    assert(res < 1e-10, 'CIP must vanish on global P%d polynomial, got %.3e.', degree, res);
    symErr = norm(C - C.', 'fro') / max(1, norm(C, 'fro'));
    assert(symErr < 1e-12, 'CIP matrix must be complex symmetric, got %.3e.', symErr);
    fprintf('PASSED  (res %.2e)\n', res);
end

fprintf('Test Helmholtz-CIP wrapper ... ');
[A, b, C] = assembleHelmholtzCIP2D(node, elem, bdFlag, 5, 1, 0, 2, [0.1, 0.01]);
[node2, ~] = extendMesh2D(node, elem, 2);
assert(isequal(size(A), [size(node2,1), size(node2,1)]), 'Unexpected wrapper matrix size.');
assert(isequal(size(b), [size(node2,1), 1]), 'Unexpected wrapper rhs size.');
assert(nnz(C) > 0, 'CIP matrix should be nonzero on this mesh.');
fprintf('PASSED\n');

fprintf('\n========== CIP 2D tests PASSED ==========\n');


function u = polynomialData(x, degree)
switch degree
    case 1
        u = 1 + 2*x(:,1) - x(:,2);
    case 2
        u = 1 + x(:,1) + 2*x(:,2) + 3*x(:,1).^2 - 2*x(:,1).*x(:,2) + x(:,2).^2;
    case 3
        u = 1 + x(:,1) - x(:,2) + x(:,1).^2 + 2*x(:,1).*x(:,2) - x(:,2).^2 + ...
            0.5*x(:,1).^3 - x(:,1).^2.*x(:,2) + 0.25*x(:,1).*x(:,2).^2 + x(:,2).^3;
end
end
