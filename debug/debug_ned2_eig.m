addpath(genpath(pwd));
node=[0,0;1,0;0,1]; elem=[1,2,3];
A=assembleNed2CurlCurl2D(node,elem); M=assembleNed2Mass2D(node,elem);
K=A+M;
e=eig(full(K)); fprintf('Eigenvalues of K:\n'); disp(sort(e)');
fprintf('min eig=%e max eig=%e cond=%e\n', min(e), max(e), max(e)/min(e));
% Also check singular values
s=svd(full(K)); fprintf('min sv=%e\n', min(s));
exit(0);
