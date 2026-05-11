function fem = main_task8_global_lod_correction()
%% TASK 8: Construct global LOD-corrected basis
%
% We construct two global corrected bases:
%
%   1. SPD auxiliary correction using
%
%        B_k = K + k Mbd.
%
%      Phi_B = P - Z (Z' B_k Z)^{-1} Z' B_k P.
%
%   2. Helmholtz correction using
%
%        A_k = K - k^2 M - 1i k Mbd.
%
%      Phi_A = P - Z (Z' A_k Z)^{-1} Z' A_k P.
%
% Here:
%
%   W_h = ker(C) = range(Z),
%
% and P is a right inverse of an effective independent constraint matrix:
%
%   Ceff P = I.
%
% The corrected bases satisfy
%
%   Ceff Phi_B = I,        Z' B_k Phi_B = 0,
%   Ceff Phi_A = I,        Z' A_k Phi_A = 0.
%
% Author: ChatGPT
% Style: close to Chen Long's iFEM convention

clear; clc;

%% Build constrained space from previous tasks

fem = main_task6_build_constraint_space();

K   = fem.K;
M   = fem.M;
Mbd = fem.Mbd;

Z = fem.constraints.Z;
C = fem.constraints.C;

k = fem.stek.k;

%% Rebuild Helmholtz and auxiliary matrices with the current k

Bk = K + k*Mbd;
Ak = K - k^2*M - 1i*k*Mbd;

%% Build effective independent constraints and a right inverse

coarse = buildCoarseRightInverse(C);

Ceff = coarse.Ceff;
P    = coarse.P;

fprintf('\n===== Task 8: global LOD correction =====\n');
fprintf('Global FEM DOFs N                      : %d\n', size(K,1));
fprintf('Original number of constraints          : %d\n', size(C,1));
fprintf('Effective number of constraints         : %d\n', size(Ceff,1));
fprintf('dim(W_h)                               : %d\n', size(Z,2));
fprintf('dim(multiscale space)                  : %d\n', size(P,2));
fprintf('Wave number k                          : %.6g\n', k);

%% Construct global B_k-LOD basis

lodB = buildGlobalLODBasis(P,Z,Bk);

%% Construct global A_k-LOD basis

lodA = buildGlobalLODBasis(P,Z,Ak);

%% Store

fem.lod.Ceff = Ceff;
fem.lod.P = P;
fem.lod.coarse = coarse;

fem.lod.Bk = Bk;
fem.lod.Ak = Ak;

fem.lod.PhiB = lodB.Phi;
fem.lod.QP_B = lodB.QP;
fem.lod.reducedB = lodB.reducedMatrix;

fem.lod.PhiA = lodA.Phi;
fem.lod.QP_A = lodA.QP;
fem.lod.reducedA = lodA.reducedMatrix;

%% Sanity checks

checkGlobalLODBasis(fem,lodB,lodA);

%% Optional visualization

plotSomeLODBasisFunctions(fem);

end


%% ------------------------------------------------------------------------
function coarse = buildCoarseRightInverse(C)
%%BUILDCOARSERIGHTINVERSE Build an independent constraint matrix Ceff and P.
%
% If C has full row rank, we keep Ceff = C and construct
%
%     P = C' (C C')^{-1},
%
% so that C P = I.
%
% If C has dependent rows, we replace C by an orthonormal row-space basis
% Ceff, obtained from the SVD:
%
%     C = U S V',
%     Ceff = V_r',
%
% so that
%
%     ker(Ceff) = ker(C),
%     Ceff Ceff' = I,
%     P = Ceff',
%     Ceff P = I.
%
% In the second case, the coarse coordinates are orthonormal combinations
% of the original constraints.

[m,N] = size(C);

if m == 0
    Ceff = sparse(0,N);
    P = sparse(N,0);

    coarse.Ceff = Ceff;
    coarse.P = P;
    coarse.rankC = 0;
    coarse.mode = 'no constraints';
    coarse.svals = [];
    return;
end

Cfull = full(C);
svals = svd(Cfull);
tol = max(size(Cfull))*eps(max(svals));
rankC = nnz(svals > tol);

if rankC == m
    % Full row rank: preserve the original constraints.
    Ceff = C;

    G = Ceff*Ceff';
    G = 0.5*(G+G');

    P = Ceff' * (G \ speye(m));

    mode = 'original full-row-rank constraints';

else
    % Dependent constraints: use an orthonormal row-space basis.
    [~,~,V] = svd(Cfull,'econ');

    Ceff = V(:,1:rankC)';
    P = Ceff';

    mode = 'orthonormalized row-space constraints';

    warning(['Constraint matrix C is rank deficient. ', ...
             'Using an orthonormal row-space basis Ceff instead of C.']);
end

coarse.Ceff = sparse(Ceff);
coarse.P = sparse(P);
coarse.rankC = rankC;
coarse.mode = mode;
coarse.svals = svals;

end


%% ------------------------------------------------------------------------
function lod = buildGlobalLODBasis(P,Z,A)
%%BUILDGLOBALLODBASIS Build corrected basis using matrix A.
%
% Given:
%
%     W_h = range(Z),
%
% and an initial right inverse P, compute
%
%     Phi = P - Z (Z' A Z)^{-1} Z' A P.
%
% This ensures
%
%     Z' A Phi = 0.
%
% If A is Hermitian SPD, this is the usual orthogonal projection.
% If A is the Helmholtz matrix, this is the Petrov/Galerkin correction
% induced by the sesquilinear form.

if isempty(P)
    lod.Phi = P;
    lod.QP = sparse(size(Z,1),0);
    lod.reducedMatrix = sparse(0,0);
    return;
end

if isempty(Z)
    Phi = P;
    QP = sparse(size(P,1),size(P,2));
    reducedMatrix = sparse(0,0);

    lod.Phi = Phi;
    lod.QP = QP;
    lod.reducedMatrix = reducedMatrix;
    return;
end

reducedMatrix = Z'*A*Z;
rhs = Z'*A*P;

Y = reducedMatrix \ rhs;

QP = Z*Y;
Phi = P - QP;

lod.Phi = Phi;
lod.QP = QP;
lod.reducedMatrix = reducedMatrix;

end


%% ------------------------------------------------------------------------
function checkGlobalLODBasis(fem,lodB,lodA)
%%CHECKGLOBALLODBASIS Verify constraints and orthogonality.

Ceff = fem.lod.Ceff;
P    = fem.lod.P;
Z    = fem.constraints.Z;

Bk = fem.lod.Bk;
Ak = fem.lod.Ak;

PhiB = lodB.Phi;
PhiA = lodA.Phi;

fprintf('\n===== Task 8 sanity checks =====\n');

%% Right inverse check

if ~isempty(P)
    rightInvErr = norm(Ceff*P - speye(size(Ceff,1)),'fro');
else
    rightInvErr = 0;
end

fprintf('right inverse error ||Ceff P - I||       : %.3e\n', rightInvErr);

%% Constraint preservation

if ~isempty(PhiB)
    conB = norm(Ceff*PhiB - speye(size(Ceff,1)),'fro');
else
    conB = 0;
end

if ~isempty(PhiA)
    conA = norm(Ceff*PhiA - speye(size(Ceff,1)),'fro');
else
    conA = 0;
end

fprintf('constraint error Ceff PhiB - I           : %.3e\n', conB);
fprintf('constraint error Ceff PhiA - I           : %.3e\n', conA);

%% Original constraint matrix check

C = fem.constraints.C;

if ~isempty(C)
    rowConsistencyB = norm(C*PhiB - C*P,'fro')/max(1,norm(C*P,'fro'));
    rowConsistencyA = norm(C*PhiA - C*P,'fro')/max(1,norm(C*P,'fro'));

    fprintf('original C row-coordinate error, PhiB    : %.3e\n', rowConsistencyB);
    fprintf('original C row-coordinate error, PhiA    : %.3e\n', rowConsistencyA);
end

%% Orthogonality checks

if ~isempty(Z)
    orthoB = norm(Z'*Bk*PhiB,'fro')/max(1,norm(Bk*PhiB,'fro'));
    orthoA = norm(Z'*Ak*PhiA,'fro')/max(1,norm(Ak*PhiA,'fro'));

    fprintf('B_k-orthogonality ||Z'' Bk PhiB||        : %.3e\n', orthoB);
    fprintf('A_k-orthogonality ||Z'' Ak PhiA||        : %.3e\n', orthoA);
else
    fprintf('Orthogonality checks skipped because W_h is empty.\n');
end

%% Reduced multiscale matrices

Bms = PhiB'*Bk*PhiB;
Ams = PhiA'*Ak*PhiA;

symBms = norm(Bms-Bms','fro')/max(1,norm(Bms,'fro'));
fprintf('symmetry error of PhiB'' Bk PhiB         : %.3e\n', symBms);

fprintf('size of multiscale B matrix              : %d x %d\n', size(Bms,1), size(Bms,2));
fprintf('size of multiscale A matrix              : %d x %d\n', size(Ams,1), size(Ams,2));

%% Condition numbers of reduced fine corrector matrices

if ~isempty(lodB.reducedMatrix)
    eigB = eig(full(0.5*(lodB.reducedMatrix+lodB.reducedMatrix')));
    eigB = real(eigB);
    eigB = eigB(eigB > 1e-14*max(eigB));

    if ~isempty(eigB)
        condB = max(eigB)/min(eigB);
        fprintf('cond(Z'' Bk Z)                           : %.3e\n', condB);
    end
end

if ~isempty(lodA.reducedMatrix)
    sA = svd(full(lodA.reducedMatrix));
    condA = max(sA)/min(sA);
    fprintf('cond(Z'' Ak Z)                           : %.3e\n', condA);
end

fprintf('\nTask 8 sanity check completed.\n');

end


%% ------------------------------------------------------------------------
function plotSomeLODBasisFunctions(fem)
%%PLOTSOMELODBASISFUNCTIONS Plot a few corrected basis functions.

node = fem.node;
elem = fem.elem;

PhiB = fem.lod.PhiB;
PhiA = fem.lod.PhiA;

nBasis = size(PhiB,2);

if nBasis == 0
    return;
end

nPlot = min(3,nBasis);

for j = 1:nPlot
    figure;
    trisurf(elem,node(:,1),node(:,2),real(PhiB(:,j)));
    shading interp;
    view(2);
    axis equal tight;
    colorbar;
    title(sprintf('Real part of B_k-LOD basis function %d',j));

    figure;
    trisurf(elem,node(:,1),node(:,2),real(PhiA(:,j)));
    shading interp;
    view(2);
    axis equal tight;
    colorbar;
    title(sprintf('Real part of A_k-LOD basis function %d',j));
end

end