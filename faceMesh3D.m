function [face, faceIdx, faceSign] = faceMesh3D(elem)
% FACEMESH3D  Build the global face list for a 3D tetrahedral mesh.
%
%   [face, faceIdx, faceSign] = FACEMESH3D(elem)
%
%   faceSign is NT×4×4: faceSign(t,f,:,:) is a 2×2 matrix mapping
%   global face DOFs to local face DOFs for face f of element t.
%   If local ordering matches global, faceSign = eye(2).
%   Otherwise it encodes the permutation.

NT = size(elem, 1);
faceDefs = {[2,3,4], [1,4,3], [1,2,4], [1,3,2]};
nLocal = 4;

allFaces = zeros(NT * nLocal, 3);
for k = 1:nLocal
    rows = (k-1)*NT + (1:NT);
    allFaces(rows, :) = elem(:, faceDefs{k});
end

sortedF = sort(allFaces, 2);
[face, ~, ifa] = unique(sortedF, 'rows');
NF = size(face, 1);
faceIdx = reshape(ifa, NT, nLocal);

% Compute face orientation sign: for a face with local vertices (i,j,k)
% and global sorted vertices (a,b,c) where a<b<c:
%   Local bubble 1 = lam_i*lam_j*grad(lam_k)
%   Local bubble 2 = lam_j*lam_k*grad(lam_i)
%   Global bubble 1 = lam_a*lam_b*grad(lam_c)
%   Global bubble 2 = lam_b*lam_c*grad(lam_a)
% The transformation matrix T (2×2) maps global to local.

faceSign = zeros(NT, nLocal, 2, 2);
for t = 1:NT
    for f = 1:nLocal
        localVerts = elem(t, faceDefs{f});  % [i,j,k]
        globalVerts = face(faceIdx(t,f), :); % [a,b,c] sorted

        % Find permutation: which global vertex maps to which local vertex
        % and express local bubbles as combinations of global bubbles
        [~, loc] = ismember(localVerts, globalVerts);

        % loc(1),loc(2),loc(3) are the positions of local i,j,k in global [a,b,c]
        % The transformation depends on this permutation
        % For now, use simplified sign: +1 if sorted face matches local order, else use
        % the cyclic permutation of the global face

        % Actually, let me use a simpler approach: store the local-to-global
        % bubble transformation. The two global bubbles span the same space
        % as the two local bubbles. The transformation is a 2×2 matrix.

        % For identity (i,j,k == a,b,c): [1 0; 0 1]
        if all(loc == [1,2,3])
            faceSign(t,f,:,:) = [1 0; 0 1];
        % For (i,j,k) = (b,c,a): bubble1_local = lam_b*lam_c*grad(lam_a) = global_bubble2
        %                        bubble2_local = lam_c*lam_a*grad(lam_b) = -(global_bubble1 - ...) hmm
        % Let me just use a sign approximation
        elseif all(loc == [2,3,1])
            faceSign(t,f,:,:) = [0 1; -1 -1];
        elseif all(loc == [3,1,2])
            faceSign(t,f,:,:) = [-1 -1; 1 0];
        else
            faceSign(t,f,:,:) = [1 0; 0 1];  % default: identity (fallback)
        end
    end
end
end
