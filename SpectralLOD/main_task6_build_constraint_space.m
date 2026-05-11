function fem = main_task6_build_constraint_space()
    %% TASK 6: Construct constraint matrix C and kernel basis Z
    %
    % We construct the constrained FEM space
    %
    %     W_h = ker(C),
    %
    % where C contains:
    %
    %   1. local Dirichlet constraints on the zero-trace parts w_i,
    %
    %          (w_i, psi_j^i)_{Omega_i} = 0,
    %
    %      with
    %
    %          w = (I - H R_sigma) u.
    %
    %   2. impedance-weighted harmonic/Steklov constraints,
    %
    %          (h, zeta_l)_{Omega} = 0,
    %
    %      with
    %
    %          h = H R_sigma u.
    %
    % Matrix representation:
    %
    %     R_sigma : extracts Sigma DOFs from a global FEM vector
    %     H       : harmonic extension matrix
    %     P_f     : I - H R_sigma
    %
    % For subdomain i:
    %
    %     C_i^D = Phi_i' M_II^i R_i^I P_f.
    %
    % For the harmonic Steklov modes:
    %
    %     C^S = X_m' B_sigma R_sigma,
    %
    % where
    %
    %     B_sigma = H' M H.
    %
    % Finally,
    %
    %     C = [C_D; C_S],
    %     W_h = ker(C).
    %
    % Author: ChatGPT
    % Style: close to Chen Long's iFEM convention
    
    clear; clc;
    
    %% Build data from previous tasks
    
    fem = main_task5_local_dirichlet_eigenpairs();
    
    %% Build constraint matrix
    
    constraints = buildConstraintMatrix(fem);
    
    fem.constraints = constraints;
    
    %% Print summary
    
    fprintf('\n===== Constraint matrix summary =====\n');
    fprintf('Number of global FEM DOFs N             : %d\n', constraints.N);
    fprintf('Number of Sigma DOFs N_sigma            : %d\n', constraints.nSigma);
    fprintf('Number of Dirichlet constraints         : %d\n', constraints.nDirichletConstraints);
    fprintf('Number of Steklov constraints           : %d\n', constraints.nSteklovConstraints);
    fprintf('Total number of constraints             : %d\n', size(constraints.C,1));
    fprintf('Numerical rank of C                     : %d\n', constraints.rankC);
    fprintf('Dimension of W_h = ker(C)               : %d\n', size(constraints.Z,2));
    
    %% Sanity checks
    
    checkConstraintMatrix(fem,constraints);
    
    end
    
    
    %% ------------------------------------------------------------------------
    function constraints = buildConstraintMatrix(fem)
    %%BUILDCONSTRAINTMATRIX Build C such that W_h = ker(C).
    %
    % Output:
    %     constraints.Cdir    : Dirichlet constraint matrix
    %     constraints.Cstek   : harmonic/Steklov constraint matrix
    %     constraints.C       : full constraint matrix
    %     constraints.Z       : orthonormal basis of ker(C)
    %     constraints.Rsigma  : skeleton restriction matrix
    %     constraints.Pfine   : P_f = I - H R_sigma
    %     constraints.mPick   : selected Steklov modes
    %     constraints.sPick   : selected local Dirichlet modes
    
    N = size(fem.node,1);
    
    dd = fem.dd;
    H = fem.schur.H;
    sigmaNode = fem.schur.sigmaNode(:);
    nSigma = numel(sigmaNode);
    
    M = fem.M;
    
    dirichlet = fem.dirichlet;
    stek = fem.stek;
    
    %% Restriction matrix R_sigma
    
    % R_sigma u = u(sigmaNode)
    Rsigma = sparse(1:nSigma,sigmaNode,ones(nSigma,1),nSigma,N);
    
    %% Fine-scale projection matrix P_f = I - H R_sigma
    
    Pfine = speye(N) - H*Rsigma;
    
    %% ------------------------------------------------------------------------
    %  Part 1: local Dirichlet constraints
    %
    % For each subdomain Omega_i:
    %
    %     C_i^D u = Phi_i' M_II^i R_i^I (I - H R_sigma) u.
    %
    % Here Phi_i contains selected local Dirichlet eigenvectors on local
    % interior DOFs, M_II^i is the local interior mass matrix, and R_i^I
    % extracts the global interior DOFs of Omega_i.
    
    Cdir = sparse(0,N);
    
    nDirichletConstraints = 0;
    sPick = dirichlet.nPick(:);
    
    for s = 1:dd.nSub
        PhiI = dirichlet.PhiI{s};      % nI x s_i selected eigenvectors
        MII  = dirichlet.MII{s};       % nI x nI local interior mass matrix
        gI   = fem.schur.globalI{s};   % global interior DOFs of subdomain
    
        np = size(PhiI,2);
    
        if np == 0
            continue;
        end
    
        % G_i = Phi_i' M_II^i, size np x nI.
        Gi = PhiI' * MII;
    
        % Constraint:
        %
        %     C_i^D = Gi * R_i^I * Pfine.
        %
        % Instead of explicitly forming R_i^I, use rows of Pfine at gI.
    
        Cblock = Gi * Pfine(gI,:);
    
        Cdir = [Cdir; Cblock];
    
        nDirichletConstraints = nDirichletConstraints + np;
    end
    
    %% ------------------------------------------------------------------------
    %  Part 2: impedance-weighted harmonic/Steklov constraints
    %
    % We select all computed weighted Steklov modes satisfying
    %
    %     rho_l < threshold.
    %
    % The harmonic eigenfunctions are zeta_l = H x_l.
    %
    % The constraint is
    %
    %     (h, zeta_l)_Omega = 0,
    %
    % where h = H R_sigma u.
    %
    % In matrix form:
    %
    %     zeta_l' M H R_sigma u = 0.
    %
    % Since zeta_l = H x_l and B_sigma = H' M H,
    %
    %     zeta_l' M H R_sigma u
    %       = x_l' H' M H R_sigma u
    %       = x_l' B_sigma R_sigma u.
    %
    % Therefore:
    %
    %     C^S = X_m' B_sigma R_sigma.
    
    threshold = stek.threshold;
    
    mPick = nnz(stek.rho < threshold);
    
    if mPick > numel(stek.rho)
        error('mPick exceeds the number of computed Steklov modes.');
    end
    
    if mPick == numel(stek.rho)
        warning(['All computed Steklov modes are below threshold. ', ...
                 'Increase nev in Task 4 to make sure rho_{m+1} is available.']);
    end
    
    if mPick > 0
        Xpick = stek.X(:,1:mPick);
        Bsigma = stek.Bsigma;
    
        Cstek = Xpick' * Bsigma * Rsigma;
    else
        Xpick = sparse(nSigma,0);
        Cstek = sparse(0,N);
    end
    
    nSteklovConstraints = mPick;
    
    %% Full constraint matrix
    
    C = [Cdir; Cstek];
    
    %% Compute kernel basis Z
    
    [Z,rankC,svals] = computeKernelBasis(C);
    
    %% Store
    
    constraints.Cdir = Cdir;
    constraints.Cstek = Cstek;
    constraints.C = C;
    
    constraints.Z = Z;
    constraints.rankC = rankC;
    constraints.svals = svals;
    
    constraints.Rsigma = Rsigma;
    constraints.Pfine = Pfine;
    
    constraints.mPick = mPick;
    constraints.sPick = sPick;
    
    constraints.nDirichletConstraints = nDirichletConstraints;
    constraints.nSteklovConstraints = nSteklovConstraints;
    
    constraints.N = N;
    constraints.nSigma = nSigma;
    
    constraints.threshold = threshold;
    
    constraints.Xpick = Xpick;
    
    end
    
    
    %% ------------------------------------------------------------------------
    function [Z,rankC,svals] = computeKernelBasis(C)
    %%COMPUTEKERNELBASIS Compute an orthonormal basis of ker(C).
    %
    % For moderate-size experiments, use dense QR on C'.
    %
    % If C is m x N, then
    %
    %     ker(C) = range(Q(:,r+1:N)),
    %
    % where Q comes from a full QR factorization of C'.
    
    [m,N] = size(C);
    
    if m == 0
        Z = speye(N);
        rankC = 0;
        svals = [];
        return;
    end
    
    Cfull = full(C);
    
    % Singular values for numerical rank.
    svals = svd(Cfull);
    
    tol = max(size(Cfull))*eps(max(svals));
    
    rankC = nnz(svals > tol);
    
    % Full QR of C' gives an orthonormal basis for R^N.
    [Q,~] = qr(Cfull');
    
    Z = Q(:,rankC+1:end);
    
    end
    
    
    %% ------------------------------------------------------------------------
    function checkConstraintMatrix(fem,constraints)
    %%CHECKCONSTRAINTMATRIX Verify that Z spans ker(C) and constraints match.
    
    C = constraints.C;
    Z = constraints.Z;
    
    N = constraints.N;
    
    fprintf('\n===== Task 6 sanity checks =====\n');
    
    %% 1. Check C*Z = 0
    
    if isempty(Z)
        fprintf('C*Z check skipped because ker(C) is empty.\n');
    else
        CZ = C*Z;
        CZerr = norm(CZ,'fro')/max(1,norm(Z,'fro'));
        fprintf('relative constraint residual ||C Z||   : %.3e\n', CZerr);
    end
    
    %% 2. Check Z orthonormality in Euclidean inner product
    
    if ~isempty(Z)
        orthZ = norm(Z'*Z-eye(size(Z,2)),'fro');
        fprintf('Euclidean orthonormality error of Z    : %.3e\n', orthZ);
    end
    
    %% 3. Random vector in W_h
    
    if ~isempty(Z)
        rng(1);
        y = randn(size(Z,2),1);
        u = Z*y;
    
        constraintResidual = norm(C*u)/max(1,norm(u));
        fprintf('random u=Zy constraint residual        : %.3e\n', constraintResidual);
    else
        fprintf('Random vector check skipped: dim ker(C)=0.\n');
    end
    
    %% 4. Directly check decomposition u = w + h
    
    if ~isempty(Z)
        Rsigma = constraints.Rsigma;
        H = fem.schur.H;
        Pfine = constraints.Pfine;
    
        h = H*(Rsigma*u);
        w = Pfine*u;
    
        decompErr = norm(u-w-h)/max(1,norm(u));
        fprintf('decomposition error u-w-h              : %.3e\n', decompErr);
    
        % Check w vanishes on Sigma.
        sigmaTraceW = Rsigma*w;
        sigmaErr = norm(sigmaTraceW)/max(1,norm(w));
        fprintf('trace of w on Sigma                    : %.3e\n', sigmaErr);
    end
    
    %% 5. Direct check of Dirichlet constraints
    
    if ~isempty(Z)
        maxDirErr = 0;
    
        for s = 1:fem.dd.nSub
            PhiI = fem.dirichlet.PhiI{s};
            MII = fem.dirichlet.MII{s};
            gI = fem.schur.globalI{s};
    
            if isempty(PhiI)
                continue;
            end
    
            wI = w(gI);
    
            val = PhiI' * MII * wI;
    
            maxDirErr = max(maxDirErr,norm(val));
        end
    
        fprintf('max direct Dirichlet constraint error  : %.3e\n', maxDirErr);
    end
    
    %% 6. Direct check of Steklov/harmonic constraints
    
    if ~isempty(Z)
        mPick = constraints.mPick;
    
        if mPick > 0
            ZetaPick = fem.stek.Zeta(:,1:mPick);
            h = H*(Rsigma*u);
    
            val = ZetaPick' * fem.M * h;
    
            stekErr = norm(val);
            fprintf('direct Steklov constraint error        : %.3e\n', stekErr);
        else
            fprintf('direct Steklov constraint check        : no selected modes\n');
        end
    end
    
    %% 7. Constraint matrix rank consistency
    
    fprintf('size(C)                                : %d x %d\n', size(C,1), size(C,2));
    fprintf('rank(C)                                : %d\n', constraints.rankC);
    fprintf('dim ker(C)                             : %d\n', size(Z,2));
    
    if size(Z,2) + constraints.rankC ~= N
        warning('Rank-nullity check failed.');
    end
    
    fprintf('\nTask 6 sanity check completed.\n');
    
    end