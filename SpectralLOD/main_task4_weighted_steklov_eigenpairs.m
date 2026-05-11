function fem = main_task4_weighted_steklov_eigenpairs()
%% TASK 4: Impedance-weighted Steklov/harmonic eigenpairs
%
% We solve the generalized eigenproblem
%
%     (S_sigma + k Mbd_sigma) x = rho (H' M H) x,
%
% where
%
%     S_sigma     = assembled full-skeleton Schur complement,
%     H           = discrete harmonic extension matrix,
%     Mbd_sigma   = Mbd restricted to Sigma DOFs,
%     H' M H      = volume mass of harmonic extensions.
%
% The harmonic eigenfunction in the volume is
%
%     zeta = H*x.
%
% The eigenvectors are normalized so that
%
%     zeta_i' M zeta_j = delta_ij,
%
% equivalently
%
%     x_i' (H' M H) x_j = delta_ij.
%
% Author: ChatGPT
% Style: close to Chen Long's iFEM convention

clear; clc;

%% Parameters

k = 20;              % wave number for the weighted eigenproblem
C0 = 4;              % threshold constant: rho < C0*k^2
nev = 70;            % number of eigenpairs to compute
useFullEig = true;   % true for small tests; false uses eigs

plotFlag = true;

%% Build data from previous tasks

fem = main_task3_schur_harmonic_extension();

% Override the wave number if needed.
fem.k = k;

%% Solve impedance-weighted Steklov eigenproblem

stek = solveWeightedSteklovEigenproblem(fem,k,nev,useFullEig);

fem.stek = stek;

%% Count modes below threshold C0*k^2

threshold = C0*k^2;

idxLow = find(stek.rho < threshold);
mPick = numel(idxLow);

fprintf('\n===== Impedance-weighted Steklov eigenproblem =====\n');
fprintf('Wave number k                         : %.6g\n', k);
fprintf('Threshold C0*k^2                      : %.6e\n', threshold);
fprintf('Number of computed eigenpairs          : %d\n', numel(stek.rho));
fprintf('Number of modes below threshold        : %d\n', mPick);
fprintf('Smallest eigenvalue rho_1              : %.6e\n', stek.rho(1));
fprintf('Largest computed eigenvalue            : %.6e\n', stek.rho(end));

if mPick < numel(stek.rho)
    fprintf('rho_{m+1} after threshold              : %.6e\n', stek.rho(mPick+1));
else
    fprintf('rho_{m+1} after threshold              : not computed\n');
end

%% Sanity checks

checkWeightedSteklovEigenpairs(fem,stek);

%% Plot eigenvalue distribution

if plotFlag
    plotWeightedSteklovEigenvalues(stek.rho,k,C0);
end

end


%% ------------------------------------------------------------------------
function stek = solveWeightedSteklovEigenproblem(fem,k,nev,useFullEig)
%%SOLVEWEIGHTEDSTEKLOVEIGENPROBLEM Solve
%
%     A_sigma x = rho B_sigma x,
%
% with
%
%     A_sigma = S_sigma + k Mbd_sigma,
%     B_sigma = H' M H.
%
% Input:
%     fem       : data structure from Task 3
%     k         : wave number
%     nev       : number of requested eigenpairs
%     useFullEig: if true, use eig; otherwise use eigs
%
% Output:
%     stek.rho        : eigenvalues
%     stek.X          : skeleton eigenvectors
%     stek.Zeta       : harmonic extensions H*X
%     stek.Asigma     : left-hand matrix
%     stek.Bsigma     : right-hand mass matrix
%     stek.MbdSigma   : boundary mass restricted to Sigma
%     stek.threshold  : default threshold 4*k^2

H = fem.schur.H;
Ssigma = fem.schur.Ssigma;
M = fem.M;
Mbd = fem.Mbd;
sigmaNode = fem.schur.sigmaNode;

nSigma = numel(sigmaNode);

% Boundary mass restricted to Sigma.
MbdSigma = Mbd(sigmaNode,sigmaNode);

% Harmonic volume mass.
Bsigma = H'*M*H;

% Weighted Steklov matrix.
Asigma = Ssigma + k*MbdSigma;

% Symmetrize to remove tiny roundoff asymmetry.
Asigma = 0.5*(Asigma + Asigma');
Bsigma = 0.5*(Bsigma + Bsigma');

% Safety check: Bsigma should be SPD.
if nSigma <= 500
    eigB = eig(full(Bsigma));
    if min(eigB) <= 1e-12*max(eigB)
        warning('Bsigma may be nearly singular. min/max eig(Bsigma) = %.3e', ...
            min(eigB)/max(eigB));
    end
end

if useFullEig || nSigma <= nev + 5
    %% Full generalized eigenvalue solve

    [X,D] = eig(full(Asigma),full(Bsigma));
    rho = real(diag(D));

    [rho,idx] = sort(rho,'ascend');
    X = X(:,idx);

    % Remove tiny negative eigenvalues caused by roundoff.
    rho(abs(rho) < 1e-12*max(1,abs(rho(end)))) = 0;

    % Keep at most nev eigenpairs.
    nevEff = min(nev,numel(rho));
    rho = rho(1:nevEff);
    X = X(:,1:nevEff);

else
    %% Sparse generalized eigenvalue solve

    opts.isreal = true;
    opts.issym = true;
    opts.tol = 1e-10;
    opts.maxit = 1000;
    opts.disp = 0;

    nevEff = min(nev,nSigma-2);

    try
        [X,D] = eigs(Asigma,Bsigma,nevEff,'smallestreal',opts);
    catch
        % Older MATLAB versions may not support 'smallestreal'.
        [X,D] = eigs(Asigma,Bsigma,nevEff,'sm',opts);
    end

    rho = real(diag(D));
    [rho,idx] = sort(rho,'ascend');
    X = X(:,idx);
end

%% Normalize eigenvectors in Bsigma inner product

for j = 1:size(X,2)
    nj = sqrt(abs(X(:,j)'*Bsigma*X(:,j)));
    X(:,j) = X(:,j)/nj;
end

% Re-orthogonalize mildly using Cholesky of Gram matrix if necessary.
G = X'*Bsigma*X;
orthErr = norm(G-eye(size(G)),'fro');

if orthErr > 1e-8
    [R,flag] = chol(G);
    if flag == 0
        X = X/R;
    else
        warning('B-orthogonalization failed because Gram matrix is not SPD.');
    end
end

%% Harmonic extensions

Zeta = H*X;

%% Store

stek.rho = rho;
stek.X = X;
stek.Zeta = Zeta;

stek.Asigma = sparse(Asigma);
stek.Bsigma = sparse(Bsigma);
stek.MbdSigma = sparse(MbdSigma);

stek.k = k;
stek.threshold = 4*k^2;

end


%% ------------------------------------------------------------------------
function checkWeightedSteklovEigenpairs(fem,stek)
%%CHECKWEIGHTEDSTEKLOVEIGENPAIRS Verify eigenproblem and normalization.

Asigma = stek.Asigma;
Bsigma = stek.Bsigma;
X = stek.X;
rho = stek.rho;

M = fem.M;
K = fem.K;
Mbd = fem.Mbd;
H = fem.schur.H;

Zeta = stek.Zeta;

fprintf('\n===== Task 4 sanity checks =====\n');

%% Symmetry checks

symA = norm(Asigma-Asigma','fro')/max(1,norm(Asigma,'fro'));
symB = norm(Bsigma-Bsigma','fro')/max(1,norm(Bsigma,'fro'));

fprintf('symmetry error Asigma                 : %.3e\n', symA);
fprintf('symmetry error Bsigma                 : %.3e\n', symB);

%% Generalized eigenproblem residuals

nEig = numel(rho);
res = zeros(nEig,1);

for j = 1:nEig
    rj = Asigma*X(:,j) - rho(j)*Bsigma*X(:,j);
    denom = max(1,norm(Asigma*X(:,j)) + abs(rho(j))*norm(Bsigma*X(:,j)));
    res(j) = norm(rj)/denom;
end

fprintf('max generalized eigen residual        : %.3e\n', max(res));
fprintf('mean generalized eigen residual       : %.3e\n', mean(res));

%% B-orthonormality

Gx = X'*Bsigma*X;
orthX = norm(Gx-eye(size(Gx)),'fro');
fprintf('B_sigma orthonormality error          : %.3e\n', orthX);

Gz = Zeta'*M*Zeta;
orthZ = norm(Gz-eye(size(Gz)),'fro');
fprintf('volume M orthonormality of H*X        : %.3e\n', orthZ);

%% Check Rayleigh quotient in volume form

nCheck = min(nEig,10);
rqErr = zeros(nCheck,1);

for j = 1:nCheck
    z = Zeta(:,j);

    numerator = z'*K*z + stek.k*(z'*Mbd*z);
    denominator = z'*M*z;
    rq = real(numerator/denominator);

    rqErr(j) = abs(rq-rho(j))/max(1,abs(rho(j)));
end

fprintf('max Rayleigh quotient error first 10  : %.3e\n', max(rqErr));

%% Check eigenvalues are nonnegative

minRho = min(rho);
fprintf('minimum computed rho                  : %.6e\n', minRho);

if minRho < -1e-8
    warning('Computed negative eigenvalue detected.');
end

fprintf('\nTask 4 sanity check completed.\n');

end


%% ------------------------------------------------------------------------
function plotWeightedSteklovEigenvalues(rho,k,C0)
%%PLOTWEIGHTEDSTEKLOVEIGENVALUES Plot eigenvalues and threshold.

threshold = C0*k^2;

figure;
semilogy(1:numel(rho),rho,'o-','LineWidth',1.2);
hold on;
yline(threshold,'--','LineWidth',1.5);
grid on;

xlabel('mode index \ell');
ylabel('\rho_\ell^{(k)}');
title(sprintf('Impedance-weighted Steklov eigenvalues, k = %.3g',k));
legend('\rho_\ell^{(k)}',sprintf('threshold %.1f k^2',C0), ...
       'Location','best');

hold off;

end