addpath(genpath(pwd));
[nd,el,bd]=squaremesh([0,1,0,1],0.125);
A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el); K=A+M;
fprintf('A: sym=%d, min(diag)=%.4e, max(diag)=%.4e\n', issymmetric(A), full(min(diag(A))), full(max(diag(A))));
fprintf('M: sym=%d, min(diag)=%.4e, max(diag)=%.4e\n', issymmetric(M), full(min(diag(M))), full(max(diag(M))));
fprintf('K: sym=%d\n', issymmetric(K));
% Check few eigenvalues
NE=size(A,1);
[~,eidx,~]=edgeMesh2D(el); bc=[2,3,1]; NT=size(el,1);
eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
% Get boundary
bdFlag_to_e=[2,3,1]; bdE=[];
for k=1:3, isBd=bd(:,k)==1; if any(isBd), bdE=[bdE;eidx(isBd,bdFlag_to_e(k))]; end; end
bdE=unique(bdE); freeE=setdiff(1:NE,bdE)';
Kff=K(freeE,freeE);
% Check eigenvalues via eigs
opts.issym=1; opts.tol=1e-4;
[V,D]=eigs(Kff,5,'smallestabs',opts);
fprintf('Smallest eigenvalues of K_free:\n'); disp(diag(D)');
fprintf('Condition number (est): %.2e\n', condest(Kff));

% Check if A+M*u produces zero for constant curl fields
% For NE_1, constant curl = gradient of a linear function.
% The curl-curl kernel contains gradients of H^1_0 functions.
% Check near-nullspace
[V2,D2]=eigs(A(freeE,freeE),3,'smallestabs',opts);
fprintf('Smallest eigenvalues of A_free:\n'); disp(diag(D2)');
exit(0);
