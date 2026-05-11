function fem = main_task5_local_dirichlet_eigenpairs()
    %% TASK 5: Local Dirichlet eigenpairs on each subdomain
    %
    % For each subdomain Omega_i, solve
    %
    %     K_II^i y = lambda M_II^i y,
    %
    % where I denotes the local interior DOFs of Omega_i.
    %
    % We select modes satisfying
    %
    %     lambda < C0*k^2.
    %
    % These modes correspond to the local Dirichlet eigenmodes that should be
    % removed from the zero-trace part w_i(v).
    %
    % Author: ChatGPT
    % Style: close to Chen Long's iFEM convention
    
    clear; clc;
    
    %% Parameters
    
    C0 = 4;              % threshold constant: lambda < C0*k^2
    useFullEig = true;   % recommended for the current moderate-size tests
    plotFlag = true;
    
    %% Build data from previous tasks
    
    fem = main_task4_weighted_steklov_eigenpairs();
    
    node = fem.node;
    elem = fem.elem;
    dd   = fem.dd;
    k    = fem.stek.k;
    
    %% Solve local Dirichlet eigenproblems
    
    dirichlet = solveLocalDirichletEigenpairs(fem,C0,useFullEig);
    
    fem.dirichlet = dirichlet;
    
    %% Print summary
    
    fprintf('\n===== Local Dirichlet eigenpair summary =====\n');
    fprintf('Wave number k                         : %.6g\n', k);
    fprintf('Threshold C0*k^2                      : %.6e\n', C0*k^2);
    fprintf('Number of subdomains                  : %d\n', dd.nSub);
    
    totalPick = sum(dirichlet.nPick);
    
    fprintf('Total selected Dirichlet modes         : %d\n', totalPick);
    
    fprintf('\nPer-subdomain selected modes:\n');
    for s = 1:dd.nSub
        fprintf('  subdomain %3d: nI = %4d, selected = %4d', ...
            s, dirichlet.nInterior(s), dirichlet.nPick(s));
    
        if dirichlet.nInterior(s) > 0
            fprintf(', lambda_1 = %.4e', dirichlet.lambda{s}(1));
            if dirichlet.nPick(s) < numel(dirichlet.lambda{s})
                fprintf(', lambda_{s+1} after threshold = %.4e', ...
                    dirichlet.lambda{s}(dirichlet.nPick(s)+1));
            else
                fprintf(', lambda_{s+1} after threshold = not computed');
            end
        end
        fprintf('\n');
    end
    
    %% Sanity checks
    
    checkLocalDirichletEigenpairs(fem,dirichlet);
    
    %% Plot
    
    if plotFlag
        plotLocalDirichletEigenvalues(fem,dirichlet,C0);
    end
    
    end
    
    
    %% ------------------------------------------------------------------------
    function dirichlet = solveLocalDirichletEigenpairs(fem,C0,useFullEig)
    %%SOLVELOCALDIRICHLETEIGENPAIRS Compute local Dirichlet eigenpairs.
    %
    % Input:
    %     fem        : structure from previous tasks
    %     C0         : threshold lambda < C0*k^2
    %     useFullEig : if true, solve full local generalized eigenproblems
    %
    % Output:
    %     dirichlet.lambda{s}     : all computed eigenvalues on subdomain s
    %     dirichlet.PhiI{s}       : selected eigenvectors on local interior DOFs
    %     dirichlet.PhiLocal{s}   : selected eigenvectors in full local numbering
    %     dirichlet.PhiGlobal{s}  : selected eigenvectors in global numbering
    %     dirichlet.nPick(s)      : number of selected modes
    %     dirichlet.MII{s}        : local interior mass matrix
    %     dirichlet.KII{s}        : local interior stiffness matrix
    %     dirichlet.threshold     : C0*k^2
    
    node = fem.node;
    elem = fem.elem;
    dd = fem.dd;
    schur = fem.schur;
    
    k = fem.stek.k;
    threshold = C0*k^2;
    
    nSub = dd.nSub;
    N = size(node,1);
    
    lambda = cell(nSub,1);
    PhiI = cell(nSub,1);
    PhiLocal = cell(nSub,1);
    PhiGlobal = cell(nSub,1);
    
    Mloc = cell(nSub,1);
    KIIcell = cell(nSub,1);
    MIIcell = cell(nSub,1);
    
    nInterior = zeros(nSub,1);
    nPick = zeros(nSub,1);
    
    for s = 1:nSub
        %% Local stiffness and local mass
    
        Kloc = schur.localK{s};
        locNode = schur.localNode{s};
        lI = schur.localI{s};
        gI = schur.globalI{s};
    
        Mloc{s} = assembleLocalP1Mass(node,elem,dd.subElem{s},locNode);
    
        nloc = numel(locNode);
        nI = numel(lI);
        nInterior(s) = nI;
    
        if nI == 0
            lambda{s} = [];
            PhiI{s} = sparse(0,0);
            PhiLocal{s} = sparse(nloc,0);
            PhiGlobal{s} = sparse(N,0);
            KIIcell{s} = sparse(0,0);
            MIIcell{s} = sparse(0,0);
            continue;
        end
    
        KII = Kloc(lI,lI);
        MII = Mloc{s}(lI,lI);
    
        KII = 0.5*(KII+KII');
        MII = 0.5*(MII+MII');
    
        KIIcell{s} = KII;
        MIIcell{s} = MII;
    
        %% Solve generalized eigenproblem
    
        if useFullEig || nI <= 300
            [V,D] = eig(full(KII),full(MII));
            lam = real(diag(D));
    
            [lam,idx] = sort(lam,'ascend');
            V = V(:,idx);
    
            % Remove tiny roundoff negatives.
            scaleLam = max(1,max(abs(lam)));
            lam(abs(lam) < 1e-12*scaleLam) = 0;
    
        else
            % Sparse fallback. We compute enough low modes near the threshold.
            %
            % For robust threshold counting, full eig is safer. This branch is
            % intended only for larger experiments.
            nev = min(nI-2, max(20,ceil(0.5*nI)));
    
            opts.isreal = true;
            opts.issym = true;
            opts.tol = 1e-10;
            opts.maxit = 1000;
            opts.disp = 0;
    
            try
                [V,D] = eigs(KII,MII,nev,'smallestreal',opts);
            catch
                [V,D] = eigs(KII,MII,nev,'sm',opts);
            end
    
            lam = real(diag(D));
            [lam,idx] = sort(lam,'ascend');
            V = V(:,idx);
        end
    
        %% Normalize eigenvectors in the MII inner product
    
        for j = 1:size(V,2)
            nj = sqrt(abs(V(:,j)'*MII*V(:,j)));
            V(:,j) = V(:,j)/nj;
        end
    
        % Mild re-orthogonalization.
        G = V'*MII*V;
        errG = norm(G-eye(size(G)),'fro');
    
        if errG > 1e-8
            [R,flag] = chol(G);
            if flag == 0
                V = V/R;
            else
                warning('Subdomain %d: M-orthogonalization failed.',s);
            end
        end
    
        %% Select modes below threshold
    
        idxPick = find(lam < threshold);
        nPick(s) = numel(idxPick);
    
        VP = V(:,idxPick);
    
        % Full local representation, zero on local boundary DOFs.
        PhiLoc = sparse(nloc,nPick(s));
        PhiLoc(lI,:) = VP;
    
        % Global representation, supported in Omega_i and zero on partial Omega_i.
        PhiG = sparse(N,nPick(s));
        PhiG(gI,:) = VP;
    
        lambda{s} = lam;
        PhiI{s} = sparse(VP);
        PhiLocal{s} = sparse(PhiLoc);
        PhiGlobal{s} = sparse(PhiG);
    end
    
    %% Store
    
    dirichlet.lambda = lambda;
    dirichlet.PhiI = PhiI;
    dirichlet.PhiLocal = PhiLocal;
    dirichlet.PhiGlobal = PhiGlobal;
    
    dirichlet.Mloc = Mloc;
    dirichlet.KII = KIIcell;
    dirichlet.MII = MIIcell;
    
    dirichlet.nInterior = nInterior;
    dirichlet.nPick = nPick;
    
    dirichlet.threshold = threshold;
    dirichlet.C0 = C0;
    dirichlet.k = k;
    
    end
    
    
    %% ------------------------------------------------------------------------
    function Mloc = assembleLocalP1Mass(node,elem,elemId,locNode)
    %%ASSEMBLELOCALP1MASS Assemble local P1 mass matrix on one subdomain.
    %
    % Input:
    %     node    : global node coordinates
    %     elem    : global element connectivity
    %     elemId  : elements belonging to the subdomain
    %     locNode : global node indices used by this subdomain
    %
    % Output:
    %     Mloc    : local mass matrix in local numbering
    
    nloc = numel(locNode);
    N = size(node,1);
    
    g2l = zeros(N,1);
    g2l(locNode) = 1:nloc;
    
    elemSub = elem(elemId,:);
    locElem = g2l(elemSub);
    
    Mloc = sparse(nloc,nloc);
    
    localMassRef = [2 1 1; 1 2 1; 1 1 2]/12;
    
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
    
        localM = area*localMassRef;
    
        Mloc(vidLocal,vidLocal) = Mloc(vidLocal,vidLocal) + localM;
    end
    
    end
    
    
    %% ------------------------------------------------------------------------
    function checkLocalDirichletEigenpairs(fem,dirichlet)
    %%CHECKLOCALDIRICHLETEIGENPAIRS Sanity checks for local Dirichlet modes.
    
    dd = fem.dd;
    nSub = dd.nSub;
    
    fprintf('\n===== Task 5 sanity checks =====\n');
    
    maxResidual = 0;
    maxOrthErr = 0;
    maxSymErrK = 0;
    maxSymErrM = 0;
    
    for s = 1:nSub
        lam = dirichlet.lambda{s};
        Phi = dirichlet.PhiI{s};
    
        KII = dirichlet.KII{s};
        MII = dirichlet.MII{s};
    
        if isempty(lam) || isempty(Phi)
            continue;
        end
    
        % Symmetry.
        symK = norm(KII-KII','fro')/max(1,norm(KII,'fro'));
        symM = norm(MII-MII','fro')/max(1,norm(MII,'fro'));
    
        maxSymErrK = max(maxSymErrK,symK);
        maxSymErrM = max(maxSymErrM,symM);
    
        % Check selected eigenpairs only.
        np = size(Phi,2);
    
        if np > 0
            lamPick = lam(1:np);
    
            Res = KII*Phi - MII*Phi*diag(lamPick);
            denom = max(1,norm(KII*Phi,'fro') + norm(MII*Phi*diag(lamPick),'fro'));
            res = norm(Res,'fro')/denom;
    
            G = Phi'*MII*Phi;
            orthErr = norm(G-eye(size(G)),'fro');
    
            maxResidual = max(maxResidual,res);
            maxOrthErr = max(maxOrthErr,orthErr);
        end
    end
    
    fprintf('max symmetry error KII                 : %.3e\n', maxSymErrK);
    fprintf('max symmetry error MII                 : %.3e\n', maxSymErrM);
    fprintf('max selected eigenpair residual        : %.3e\n', maxResidual);
    fprintf('max selected MII-orthonormality error  : %.3e\n', maxOrthErr);
    
    % Check threshold logic.
    badCount = 0;
    for s = 1:nSub
        lam = dirichlet.lambda{s};
        np = dirichlet.nPick(s);
    
        if isempty(lam)
            continue;
        end
    
        if np > 0
            badCount = badCount + nnz(lam(1:np) >= dirichlet.threshold);
        end
    
        if np < numel(lam)
            if lam(np+1) < dirichlet.threshold
                badCount = badCount + 1;
            end
        end
    end
    
    fprintf('threshold-selection inconsistencies    : %d\n', badCount);
    
    if maxResidual > 1e-8
        warning('Some local Dirichlet eigenpair residuals are larger than expected.');
    end
    
    if maxOrthErr > 1e-8
        warning('Some selected Dirichlet eigenvectors are not well M-orthonormalized.');
    end
    
    if badCount > 0
        warning('Threshold selection has inconsistencies.');
    end
    
    fprintf('\nTask 5 sanity check completed.\n');
    
    end
    
    
    %% ------------------------------------------------------------------------
    function plotLocalDirichletEigenvalues(fem,dirichlet,C0)
    %%PLOTLOCALDIRICHLETEIGENVALUES Plot local Dirichlet eigenvalues.
    
    dd = fem.dd;
    k = dirichlet.k;
    threshold = C0*k^2;
    
    figure;
    hold on;
    
    for s = 1:dd.nSub
        lam = dirichlet.lambda{s};
        if isempty(lam)
            continue;
        end
    
        plot(1:numel(lam),lam,'o-','LineWidth',1.0);
    end
    
    yline(threshold,'k--','LineWidth',1.5);
    set(gca,'YScale','log');
    grid on;
    
    xlabel('local mode index j');
    ylabel('\lambda_j^i');
    title(sprintf('Local Dirichlet eigenvalues, threshold %.1f k^2',C0));
    
    hold off;
    
    figure;
    bar(1:dd.nSub,dirichlet.nPick);
    grid on;
    xlabel('subdomain index i');
    ylabel('number of selected modes');
    title('Number of selected local Dirichlet modes per subdomain');
    
    end