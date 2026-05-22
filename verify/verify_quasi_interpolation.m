% VERIFY_QUASI_INTERPOLATION  Verify weighted Clement and Scott-Zhang P1 transfers.

fprintf('========== Quasi-Interpolation Verification ==========\n\n');

fprintf('Test 1: 2D constants and Scott-Zhang linear reproduction ... ');
[cNode2, cElem2, cBd2] = squaremesh([0, 1, 0, 1], 0.5);
[fNode2, fElem2] = squaremesh([0, 1, 0, 1], 0.25);
Qc2 = weightedClementP1(fNode2, fElem2, cNode2, cElem2);
Qs2 = scottZhangP1(fNode2, fElem2, cNode2, cElem2);
assert(norm(Qc2 * ones(size(fNode2,1),1) - 1, inf) < 1e-12, ...
    '2D weighted Clement must reproduce constants.');
uFine = 1 + 2*fNode2(:,1) - 3*fNode2(:,2);
uCoarse = 1 + 2*cNode2(:,1) - 3*cNode2(:,2);
assert(norm(Qs2*uFine - uCoarse, inf) < 1e-11, ...
    '2D Scott-Zhang must reproduce coarse linear functions.');
fprintf('PASSED\n');

fprintf('Test 2: 2D Scott-Zhang preserves homogeneous boundary data ... ');
u0 = fNode2(:,1).*(1 - fNode2(:,1)).*fNode2(:,2).*(1 - fNode2(:,2));
uc0 = Qs2 * u0;
bd = getBoundaryNodes2D(cElem2, cBd2);
assert(norm(uc0(bd), inf) < 1e-12, ...
    '2D Scott-Zhang boundary coefficients must vanish for zero boundary traces.');
fprintf('PASSED\n');

fprintf('Test 3: 3D constants and Scott-Zhang linear reproduction ... ');
[cNode3, cElem3] = cubemesh([0, 1, 0, 1, 0, 1], 0.5);
[fNode3, fElem3] = cubemesh([0, 1, 0, 1, 0, 1], 0.25);
Qc3 = weightedClementP1(fNode3, fElem3, cNode3, cElem3);
Qs3 = scottZhangP1(fNode3, fElem3, cNode3, cElem3);
assert(norm(Qc3 * ones(size(fNode3,1),1) - 1, inf) < 1e-12, ...
    '3D weighted Clement must reproduce constants.');
uFine3 = 1 + fNode3(:,1) - 2*fNode3(:,2) + 0.5*fNode3(:,3);
uCoarse3 = 1 + cNode3(:,1) - 2*cNode3(:,2) + 0.5*cNode3(:,3);
assert(norm(Qs3*uFine3 - uCoarse3, inf) < 1e-10, ...
    '3D Scott-Zhang must reproduce coarse linear functions.');
fprintf('PASSED\n');

fprintf('Test 4: nested P1 prolongation consistency ... ');
P2 = prolongateNestedP1(cNode2, cElem2, fNode2);
assert(norm(P2*uCoarse - uFine, inf) < 1e-12, ...
    '2D nested P1 prolongation must reproduce linear functions.');
P3 = prolongateNestedP1(cNode3, cElem3, fNode3);
assert(norm(P3*uCoarse3 - uFine3, inf) < 1e-12, ...
    '3D nested P1 prolongation must reproduce linear functions.');
fprintf('PASSED\n');

fprintf('\n========== Quasi-Interpolation tests PASSED ==========\n');
