function fem = main_task3_schur_harmonic_extension()
%% TASK 3: Build local Schur complements and harmonic extension matrix
%
% We construct:
%
%   1. local subdomain Schur complements
%
%        S_i = K_BB^i - K_BI^i (K_II^i)^{-1} K_IB^i
%
%   2. assembled full-skeleton Schur complement
%
%        S_sigma = sum_i R_i^T S_i R_i
%
%   3. harmonic extension matrix H
%
%        H : skeleton values on Sigma --> globally conforming FEM vector
%
%      such that for each subdomain Omega_i, the interior DOFs solve
%
%        K_II^i u_I = - K_IB^i u_B.
%
%   4. consistency check
%
%        S_sigma = H^T K H
%
%      up to roundoff.
%
% Author: ChatGPT
% Style: close to Chen Long's iFEM convention

clear; clc;

%% Load mesh and global FEM matrices from Task 2

fem = main_task2_assemble_matrices();

node = fem.node;
elem = fem.elem;
mesh = fem.mesh;
dd   = fem.dd;
K    = fem.K;

%% Build local Schur complements and harmonic extension

schur = buildSubdomainSchurAndHarmonicExtension(node,elem,dd,K);

fem.schur = schur;

%% Print information

fprintf('\n===== Schur complement and harmonic extension information =====\n');
fprintf('Number of full skeleton DOFs N_sigma     : %d\n', numel(dd.sigmaNode));
fprintf('Size of H                                : %d x %d\n', size(schur.H,1), size(schur.H,2));
fprintf('Size of Ssigma                           : %d x %d\n', size(schur.Ssigma,1), size(schur.Ssigma,2));
fprintf('nnz(H)                                   : %d\n', nnz(schur.H));
fprintf('nnz(Ssigma)                              : %d\n', nnz(schur.Ssigma));

%% Sanity checks

checkSchurAndHarmonicExtension(fem,schur);

end


%% ------------------------------------------------------------------------
function schur = buildSubdomainSchurAndHarmonicExtension(node,elem,dd,Kglobal)
%%BUILDSUBDOMAINSCHURANDHARMONICEXTENSION
%
% Build local Schur complements and global harmonic extension matrix.
%
% Input:
%     node    : global node coordinates
%     elem    : global element connectivity
%     dd      : domain decomposition structure
%     Kglobal : global stiffness matrix
%
% Output:
%     schur.Slocal{s}     : local Schur complement on subdomain boundary
%     schur.localK{s}     : local stiffness matrix of subdomain
%     schur.localNode{s}  : global node indices used by subdomain
%     schur.localI{s}     : local interior indices
%     schur.localB{s}     : local boundary indices
%     schur.globalI{s}    : global interior node indices
%     schur.globalB{s}    : global boundary node indices
%     schur.Sglobal       : global-size assembled Schur complement
%     schur.Ssigma        : Schur complement restricted to Sigma DOFs
%     schur.H             : harmonic extension matrix from Sigma DOFs to all DOFs

N = size(node,1);
nSub = dd.nSub;

sigmaNode = dd.sigmaNode(:);
nSigma = numel(sigmaNode);

% Map global node index --> Sigma column index.
sigmaMap = zeros(N,1);
sigmaMap(sigmaNode) = 1:nSigma;

Sglobal = sparse(N,N);

H = sparse(N,nSigma);

% Values on Sigma are copied directly.
H(sigmaNode,1:nSigma) = speye(nSigma);

Slocal = cell(nSub,1);
localK = cell(nSub,1);
localNode = cell(nSub,1);
localI = cell(nSub,1);
localB = cell(nSub,1);
globalI = cell(nSub,1);
globalB = cell(nSub,1);

for s = 1:nSub
    %% Assemble local subdomain stiffness matrix

    elemId = dd.subElem{s};

    [Ki,locNode,locElem] = assembleLocalP1Stiffness(node,elem,elemId);

    localK{s} = Ki;
    localNode{s} = locNode;

    %% Local interior and boundary indices

    gB = dd.subBdNode{s}(:);
    gI = dd.subIntNode{s}(:);

    % Map global node index to local node index.
    locMap = zeros(N,1);
    locMap(locNode) = 1:numel(locNode);

    lB = locMap(gB);
    lI = locMap(gI);

    if any(lB == 0)
        error('Subdomain %d: some boundary nodes are not in localNode.',s);
    end
    if any(lI == 0)
        error('Subdomain %d: some interior nodes are not in localNode.',s);
    end

    localB{s} = lB;
    localI{s} = lI;
    globalB{s} = gB;
    globalI{s} = gI;

    %% Local Schur complement

    KBB = Ki(lB,lB);

    if isempty(lI)
        Si = KBB;
    else
        KBI = Ki(lB,lI);
        KIB = Ki(lI,lB);
        KII = Ki(lI,lI);

        Si = KBB - KBI*(KII\KIB);
    end

    Slocal{s} = sparse(Si);

    %% Assemble into global Schur complement

    Sglobal(gB,gB) = Sglobal(gB,gB) + Si;

    %% Harmonic extension matrix on this subdomain

    % Boundary values are already identity in H.
    % For interior nodes:
    %
    %     u_I = - K_II^{-1} K_IB u_B.
    %
    % We express u_B in terms of Sigma DOFs.
    %
    % Since every subdomain boundary node is in Sigma, sigmaMap(gB) is valid.

    if ~isempty(lI)
        KII = Ki(lI,lI);
        KIB = Ki(lI,lB);

        Bcols = sigmaMap(gB);

        if any(Bcols == 0)
            error('Subdomain %d: subdomain boundary node not found in sigmaNode.',s);
        end

        Tlocal = - KII\KIB;  % maps local boundary values to local interior values

        H(gI,Bcols) = H(gI,Bcols) + Tlocal;
    end
end

%% Restrict global Schur complement to Sigma

Ssigma = Sglobal(sigmaNode,sigmaNode);

%% Store

schur.Slocal = Slocal;
schur.localK = localK;
schur.localNode = localNode;
schur.localI = localI;
schur.localB = localB;
schur.globalI = globalI;
schur.globalB = globalB;

schur.Sglobal = Sglobal;
schur.Ssigma = sparse(Ssigma);
schur.H = sparse(H);

schur.sigmaNode = sigmaNode;
schur.sigmaMap = sigmaMap;

end


%% ------------------------------------------------------------------------
function [Kloc,locNode,locElem] = assembleLocalP1Stiffness(node,elem,elemId)
%%ASSEMBLELOCALP1STIFFNESS Assemble stiffness matrix on one subdomain.
%
% Input:
%     node   : global node coordinates
%     elem   : global element connectivity
%     elemId : elements belonging to one subdomain
%
% Output:
%     Kloc    : local stiffness matrix
%     locNode : global node indices used by this subdomain
%     locElem : element connectivity in local numbering

elemSub = elem(elemId,:);
locNode = unique(elemSub(:));

nloc = numel(locNode);

% global node index --> local node index
maxNode = size(node,1);
g2l = zeros(maxNode,1);
g2l(locNode) = 1:nloc;

locElem = g2l(elemSub);

Kloc = sparse(nloc,nloc);

for t = 1:size(locElem,1)
    vidLocal = locElem(t,:);
    vidGlobal = locNode(vidLocal);

    p = node(vidGlobal,:);

    x1 = p(1,1); y1 = p(1,2);
    x2 = p(2,1); y2 = p(2,2);
    x3 = p(3,1); y3 = p(3,2);

    detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
    area = abs(detJ)/2;

    if area <= 0
        error('Local element has nonpositive area.');
    end

    b = [y2-y3; y3-y1; y1-y2];
    c = [x3-x2; x1-x3; x2-x1];

    localK = (b*b' + c*c')/(4*area);

    Kloc(vidLocal,vidLocal) = Kloc(vidLocal,vidLocal) + localK;
end

end


%% ------------------------------------------------------------------------
function checkSchurAndHarmonicExtension(fem,schur)
%%CHECKSCHURANDHARMONICEXTENSION Consistency tests.

K = fem.K;
dd = fem.dd;

H = schur.H;
Ssigma = schur.Ssigma;
sigmaNode = schur.sigmaNode;

fprintf('\n===== Task 3 sanity checks =====\n');

%% 1. Symmetry of Ssigma

symS = norm(Ssigma-Ssigma','fro')/max(1,norm(Ssigma,'fro'));
fprintf('symmetry error Ssigma                 : %.3e\n', symS);

%% 2. Check Ssigma = H^T K H

SH = H'*K*H;

relErr = norm(Ssigma-SH,'fro')/max(1,norm(Ssigma,'fro'));
fprintf('relative error ||Ssigma - H^T K H||   : %.3e\n', relErr);

%% 3. Check harmonicity of H columns inside each subdomain

maxHarmonicResidual = 0;

for s = 1:dd.nSub
    Ki = schur.localK{s};
    lI = schur.localI{s};
    lB = schur.localB{s};
    gI = schur.globalI{s};
    gB = schur.globalB{s};

    if isempty(lI)
        continue;
    end

    KII = Ki(lI,lI);
    KIB = Ki(lI,lB);

    % Extract H restricted to local interior and local boundary.
    HI = H(gI,:);
    HB = H(gB,:);

    residual = KII*HI + KIB*HB;

    resNorm = norm(residual,'fro')/max(1,norm(KIB*HB,'fro'));
    maxHarmonicResidual = max(maxHarmonicResidual,resNorm);
end

fprintf('max local harmonic residual           : %.3e\n', maxHarmonicResidual);

%% 4. Check constants: harmonic extension of constant trace is constant

nSigma = numel(sigmaNode);
oneSigma = ones(nSigma,1);
u = H*oneSigma;

constErr = norm(u-ones(size(u)),inf);
fprintf('constant trace extension error        : %.3e\n', constErr);

%% 5. Energy identity for a random skeleton vector

rng(1);
x = randn(nSigma,1);
u = H*x;

energyVolume = u'*K*u;
energySigma = x'*Ssigma*x;

energyErr = abs(energyVolume-energySigma)/max(1,abs(energySigma));
fprintf('random energy identity error          : %.3e\n', energyErr);

%% Optional warnings

tol = 1e-9;

if symS > tol
    warning('Ssigma symmetry error is larger than tolerance.');
end

if relErr > 1e-8
    warning('Ssigma and H^T K H differ more than expected.');
end

if maxHarmonicResidual > 1e-8
    warning('Harmonic extension residual is larger than expected.');
end

if constErr > 1e-10
    warning('Constant trace is not extended to a constant function.');
end

if energyErr > 1e-8
    warning('Energy identity check failed.');
end

fprintf('\nTask 3 sanity check completed.\n');

end