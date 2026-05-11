addpath(genpath(pwd));
node=[0,0;1,0;0,1]; elem=[1,2,3];
[~,eidx,es]=edgeMesh2D(elem);
A=assembleNed2CurlCurl2D(node,elem); M=assembleNed2Mass2D(node,elem);
fprintf('A: %dx%d nnz=%d sym=%d\n',size(A,1),size(A,2),nnz(A),issymmetric(A));
fprintf('M: %dx%d nnz=%d sym=%d\n',size(M,1),size(M,2),nnz(M),issymmetric(M));
fprintf('A diag: '); fprintf('%.4e ',full(diag(A))); fprintf('\n');
fprintf('M diag: '); fprintf('%.4e ',full(diag(M))); fprintf('\n');
% Check nullspace of A
rs=sum(A,2); fprintf('A row sums: '); fprintf('%.2e ',full(rs)); fprintf('\n');
% Check interior DOFs (7,8) vs edge DOFs
fprintf('Interior DOF coupling: A(7,1:6)='); fprintf('%.4e ',full(A(7,1:6))); fprintf('\n');
exit(0);
