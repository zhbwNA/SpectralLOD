function fem = main_task7_verify_constrained_poincare()
%% TASK 7: Verify constrained Poincare-type inequality on W_h
%
% We check numerically whether
%
%     k^2 ||w||_{L2(Omega)}^2
%       <= C_W ( ||grad w||_{L2(Omega)}^2
%                + k ||w||_{L2(partial Omega)}^2 )
%
% for w in W_h = ker(C).
%
% Matrix form:
%
%     eta(w) =
%       k^2 w' M w / ( w' (K + k Mbd) w ).
%
% Since W_h = range(Z), w = Z y, the worst-case discrete constant is
%
%     eta_max =
%       lambda_max( k^2 Z' M Z, Z' (K + k Mbd) Z ).
%
% We compute:
%   1. random-sampling estimates of eta(w);
%   2. the exact reduced generalized eigenvalue eta_max;
%   3. a related Helmholtz coercivity quotient
%
%        |w' Ahelm w| / ( w' (K + k Mbd) w ).
%
% Author: ChatGPT
% Style: close to Chen Long's iFEM convention

clear; clc;

%% Parameters

nSample = 500;       % number of random samples
rngSeed = 1;         % reproducibility
doExactEig = true;   % compute exact worst-case constant on W_h

%% Build constrained space from Task 6

fem = main_task6_build_constraint_space();

K = fem.K;
M = fem.M;
Mbd = fem.Mbd;
Z = fem.constraints.Z;

k = fem.stek.k;

N = size(K,1);
nW = size(Z,2);

fprintf('\n===== Task 7: constrained Poincare verification =====\n');
fprintf('Global FEM DOFs N                    : %d\n', N);
fprintf('dim(W_h)                             : %d\n', nW);
fprintf('Wave number k                        : %.6g\n', k);

if nW == 0
    warning('The constrained space W_h is empty. Nothing to test.');
    return;
end

%% Define matrices for the inequality

Gk = K + k*Mbd;       % denominator matrix
Ak = K - k^2*M - 1i*k*Mbd;

%% Random sampling test

stats = randomTestConstrainedPoincare(Z,M,K,Mbd,Ak,k,nSample,rngSeed);

fprintf('\n===== Random sampling results =====\n');
fprintf('Number of random samples              : %d\n', nSample);
fprintf('max eta(w)                            : %.6e\n', stats.maxEta);
fprintf('mean eta(w)                           : %.6e\n', stats.meanEta);
fprintf('min eta(w)                            : %.6e\n', stats.minEta);
fprintf('std eta(w)                            : %.6e\n', stats.stdEta);

fprintf('\nHelmholtz quotient on random samples:\n');
fprintf('min |a(w,w)| / ||w||_G^2              : %.6e\n', stats.minHelmQuot);
fprintf('mean |a(w,w)| / ||w||_G^2             : %.6e\n', stats.meanHelmQuot);

%% Exact reduced generalized eigenvalue test

if doExactEig
    exact = exactConstrainedPoincareConstant(Z,M,Gk,k);

    fprintf('\n===== Exact reduced generalized eigenvalue result =====\n');
    fprintf('eta_max = max k^2 w^T M w / w^T Gk w : %.6e\n', exact.etaMax);
    fprintf('eta_min                              : %.6e\n', exact.etaMin);

    if exact.etaMax < 1
        fprintf('\nResult: eta_max < 1. The discrete constrained estimate is strong enough.\n');
    else
        fprintf('\nResult: eta_max >= 1. The chosen constraints are not strong enough for this k.\n');
        fprintf('Suggestion: increase C0, compute more Steklov modes, or refine the selection threshold.\n');
    end

    fem.poincareExact = exact;
end

%% Store

fem.poincareRandom = stats;

%% Optional plot

plotPoincareSamples(stats);

end


%% ------------------------------------------------------------------------
function stats = randomTestConstrainedPoincare(Z,M,K,Mbd,Ak,k,nSample,rngSeed)
%%RANDOMTESTCONSTRAINEDPOINCARE Randomly sample w in W_h and test quotients.

rng(rngSeed);

nW = size(Z,2);

eta = zeros(nSample,1);
helmQuot = zeros(nSample,1);

for j = 1:nSample
    y = randn(nW,1);
    w = Z*y;

    L2sq = real(w'*M*w);
    H1semiSq = real(w'*K*w);
    Bdsq = real(w'*Mbd*w);

    denom = H1semiSq + k*Bdsq;

    eta(j) = k^2*L2sq/denom;

    aval = w'*Ak*w;
    helmQuot(j) = abs(aval)/denom;
end

stats.eta = eta;
stats.helmQuot = helmQuot;

stats.maxEta = max(eta);
stats.meanEta = mean(eta);
stats.minEta = min(eta);
stats.stdEta = std(eta);

stats.minHelmQuot = min(helmQuot);
stats.meanHelmQuot = mean(helmQuot);
stats.maxHelmQuot = max(helmQuot);
stats.stdHelmQuot = std(helmQuot);

end


%% ------------------------------------------------------------------------
function exact = exactConstrainedPoincareConstant(Z,M,Gk,k)
%%EXACTCONSTRAINEDPOINCARECONSTANT Compute exact worst-case discrete constant.
%
% We compute
%
%     eta_max =
%       lambda_max( k^2 Z' M Z, Z' Gk Z ).
%
% where
%
%     Gk = K + k Mbd.
%
% This gives the sharp discrete constant on W_h.

Mr = Z'*M*Z;
Gr = Z'*Gk*Z;

Mr = 0.5*(Mr+Mr');
Gr = 0.5*(Gr+Gr');

A = k^2*Mr;

% For moderate-size experiments, use dense eig.
[V,D] = eig(full(A),full(Gr));

etaVals = real(diag(D));
etaVals = sort(etaVals,'ascend');

% Remove tiny negative roundoff values.
scale = max(1,max(abs(etaVals)));
etaVals(abs(etaVals) < 1e-12*scale) = 0;

exact.etaVals = etaVals;
exact.etaMin = etaVals(1);
exact.etaMax = etaVals(end);

% Worst-case vector in W_h.
[~,idxMax] = max(real(diag(D)));
yMax = V(:,idxMax);
wMax = Z*yMax;

% Normalize by Gk-energy.
normG = sqrt(real(wMax'*Gk*wMax));
wMax = wMax/normG;

exact.wMax = wMax;

% Check quotient directly.
exact.etaMaxCheck = real(k^2*(wMax'*M*wMax)/(wMax'*Gk*wMax));

end


%% ------------------------------------------------------------------------
function plotPoincareSamples(stats)
%%PLOTPOINCARESAMPLES Plot sampled quotients.

figure;
plot(1:numel(stats.eta),stats.eta,'o','MarkerSize',4);
grid on;
xlabel('sample index');
ylabel('\eta(w)');
title('Random samples of constrained Poincare quotient');

figure;
histogram(stats.eta,30);
grid on;
xlabel('\eta(w)');
ylabel('frequency');
title('Histogram of constrained Poincare quotient');

figure;
plot(1:numel(stats.helmQuot),stats.helmQuot,'o','MarkerSize',4);
grid on;
xlabel('sample index');
ylabel('|a(w,w)| / (||grad w||^2 + k||w||_{bd}^2)');
title('Random samples of Helmholtz coercivity quotient');

end