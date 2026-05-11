addpath(genpath(pwd));
[nd,el]=squaremesh([0,1,0,1],0.25);
M=assembleNedMass2D(nd,el);
fprintf('M diag range: [%.4e, %.4e]\n', full(min(diag(M))), full(max(diag(M))));
% Compare: old (factor 1) vs new (factor 2)?
% For h=0.25, diag entries should be O(0.5 * 2) = O(1)
fprintf('Expected diag ~ 0.3-0.7 (with factor 2 fix)\n');
exit(0);
