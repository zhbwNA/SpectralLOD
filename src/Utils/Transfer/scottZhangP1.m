function Q = scottZhangP1(fineNode, fineElem, coarseNode, coarseElem, opts)
% SCOTTZHANGP1  Matrix Scott-Zhang quasi-interpolation, P1 fine to P1 coarse.
%
%   uc = Q*uf maps fine P1 nodal values to coarse P1 nodal values.
%   Boundary coarse nodes use boundary edges/faces so homogeneous boundary
%   data are preserved.

if nargin < 5
    opts = struct();
end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = 4;
end

dim = size(coarseNode, 2);
Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);

entities = selectScottZhangEntities(coarseNode, coarseElem);

ii = zeros(64 * Nc, 1);
jj = zeros(64 * Nc, 1);
ss = zeros(64 * Nc, 1);
idx = 0;

for z = 1:Nc
    ent = entities(z);
    [lambdaE, wRef] = entityQuadrature(dim, ent.entityDim, opts.quadOrder);
    m = size(lambdaE, 2);
    verts = ent.vertices(:).';
    xVerts = coarseNode(verts, :);
    jacScale = entityJacobianScale(xVerts, ent.entityDim);

    G = zeros(m, m);
    for q = 1:length(wRef)
        phi = lambdaE(q, :);
        G = G + jacScale * wRef(q) * (phi.' * phi);
    end

    rhs = zeros(m, 1);
    rhs(ent.localPos) = 1;
    dualCoeff = G \ rhs;

    for q = 1:length(wRef)
        xq = lambdaE(q, :) * xVerts;
        psi = lambdaE(q, :) * dualCoeff;
        [tf, lambdaF] = locateSimplexP1(fineNode, fineElem, xq, 1e-10);
        if tf == 0
            error('scottZhangP1:notNested', ...
                'A Scott-Zhang quadrature point was not found in the fine mesh.');
        end
        fv = fineElem(tf, :);
        val = jacScale * wRef(q) * psi * lambdaF;
        need = idx + (1:length(fv));
        if need(end) > numel(ii)
            grow = max(numel(ii), 1024);
            ii(end+grow) = 0; %#ok<AGROW>
            jj(end+grow) = 0; %#ok<AGROW>
            ss(end+grow) = 0; %#ok<AGROW>
        end
        ii(need) = z;
        jj(need) = fv(:);
        ss(need) = val(:);
        idx = need(end);
    end
end

Q = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Nc, Nf);
end


function entities = selectScottZhangEntities(node, elem)
dim = size(node, 2);
Nc = size(node, 1);
entities = repmat(struct('vertices', [], 'entityDim', dim, 'localPos', 1), Nc, 1);

if dim == 2
    [edge, ~, ~, edge2elem] = edgeMesh2D(elem);
    bdEntity = edge(edge2elem(:,2) == 0, :);
elseif dim == 3
    bdEntity = boundaryFaces3D(elem);
else
    error('scottZhangP1:dim', 'Only 2D and 3D meshes are supported.');
end

isBdNode = false(Nc, 1);
isBdNode(unique(bdEntity(:))) = true;

incident = cell(Nc, 1);
for t = 1:size(elem, 1)
    for a = 1:size(elem, 2)
        incident{elem(t, a)}(end+1) = t; %#ok<AGROW>
    end
end

for z = 1:Nc
    if isBdNode(z)
        rows = find(any(bdEntity == z, 2), 1);
        verts = bdEntity(rows, :);
        entities(z).vertices = verts;
        entities(z).entityDim = dim - 1;
        entities(z).localPos = find(verts == z, 1);
    else
        t = incident{z}(1);
        verts = elem(t, :);
        entities(z).vertices = verts;
        entities(z).entityDim = dim;
        entities(z).localPos = find(verts == z, 1);
    end
end
end


function face = boundaryFaces3D(elem)
faceDefs = {[2 3 4], [1 3 4], [1 2 4], [1 2 3]};
NT = size(elem, 1);
allFaces = zeros(4 * NT, 3);
for f = 1:4
    allFaces((f - 1) * NT + (1:NT), :) = sort(elem(:, faceDefs{f}), 2);
end
[faceAll, ~, ic] = unique(allFaces, 'rows');
counts = accumarray(ic, 1);
face = faceAll(counts == 1, :);
end


function [lambda, weight] = entityQuadrature(dim, entityDim, order)
if entityDim == 1
    [t, weight] = gauss1D(max(2, min(6, order)));
    lambda = [1 - t, t];
elseif entityDim == 2
    [lambda, weight] = quadtriangle(order);
elseif entityDim == 3 && dim == 3
    [lambda, weight] = quadtet(order);
else
    error('scottZhangP1:entityDim', 'Unsupported entity dimension.');
end
end


function s = entityJacobianScale(v, entityDim)
if entityDim == 1
    s = norm(v(2,:) - v(1,:));
elseif entityDim == 2
    if size(v, 2) == 2
        area = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:)])) / 2;
    else
        area = norm(cross(v(2,:) - v(1,:), v(3,:) - v(1,:))) / 2;
    end
    s = 2 * area;
elseif entityDim == 3
    vol = abs(det([v(2,:) - v(1,:); v(3,:) - v(1,:); v(4,:) - v(1,:)])) / 6;
    s = 6 * vol;
else
    error('scottZhangP1:entityDim', 'Unsupported entity dimension.');
end
end
