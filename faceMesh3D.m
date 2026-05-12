function [face, faceIdx, faceTrans] = faceMesh3D(node, elem)
% FACEMESH3D  Build global face list with orientation transformations.
%
%   [face, faceIdx, faceTrans] = FACEMESH3D(node, elem)
%
%   faceTrans{t,f} is a 2×2 matrix mapping global → local face DOFs.
%   For assembly: A_global += T' * A_local(face_block) * T.

NT = size(elem, 1);
faceDefs = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
nLocal = 4;

allFaces = zeros(NT*nLocal, 3);
for k = 1:nLocal
    allFaces((k-1)*NT+(1:NT), :) = elem(:, faceDefs{k});
end
[face, ~, ifa] = unique(sort(allFaces,2), 'rows');
faceIdx = reshape(ifa, NT, nLocal);

faceTrans = cell(NT, nLocal);

for t = 1:NT
    v = node(elem(t,:), :);
    e12=v(2,:)-v(1,:); e13=v(3,:)-v(1,:); e14=v(4,:)-v(1,:);
    J = [e12; e13; e14]';
    invJ = inv(J);  % J^{-T} rows give ∇λ_2, ∇λ_3, ∇λ_4
    g2=invJ(1,:); g3=invJ(2,:); g4=invJ(3,:); g1=-(g2+g3+g4);
    Gphys = {g1, g2, g3, g4};

    for f = 1:nLocal
        lfIdx = faceDefs{f};                   % local vertex indices [i,j,k]
        localV = elem(t, lfIdx);               % global vertex numbers
        globalV = face(faceIdx(t,f), :);       % sorted global [a,b,c]

        % Map global→local: which local index (1-4) matches each global vertex
        % perm(1)=local index of sorted global vertex a, etc.
        [~, perm] = ismember(globalV, localV); % perm maps global→local_index

        % Sample: use 5 points for robust least-squares fit
        lam_vals = [1/3 1/3 1/3; 2/3 1/3 0; 0 2/3 1/3; 1/2 1/2 0; 1/4 1/2 1/4];
        nPt = 5;
        A_mat = zeros(6*nPt, 4); B_vec = zeros(6*nPt, 1);
        row = 0;

        for pt = 1:nPt
            lam = zeros(1,4);
            lam(lfIdx) = lam_vals(pt,:);
            li=lam(lfIdx(1)); lj=lam(lfIdx(2)); lk=lam(lfIdx(3));
            L1 = li*lj * Gphys{lfIdx(3)};  % ∇λ_k
            L2 = lj*lk * Gphys{lfIdx(1)};  % ∇λ_i
            la = lam(lfIdx(perm(1))); lb = lam(lfIdx(perm(2))); lc = lam(lfIdx(perm(3)));
            G1 = la*lb * Gphys{lfIdx(perm(3))};
            G2 = lb*lc * Gphys{lfIdx(perm(1))};
            for comp = 1:3
                row = row + 1;
                A_mat(row,:) = [G1(comp), G2(comp), 0, 0];
                B_vec(row) = L1(comp);
                row = row + 1;
                A_mat(row,:) = [0, 0, G1(comp), G2(comp)];
                B_vec(row) = L2(comp);
            end
        end
        T_vec = pinv(A_mat) * B_vec;       % least-squares, robust to rank deficiency
        faceTrans{t,f} = reshape(T_vec, [2,2])';
    end
end
end
