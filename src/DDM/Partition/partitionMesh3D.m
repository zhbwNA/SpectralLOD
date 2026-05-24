function parts = partitionMesh3D(node, elem, bdFlag, nSub, varargin)
% PARTITIONMESH3D  Partition a 3D tetrahedral mesh into subdomains.
%
%   parts = PARTITIONMESH3D(node, elem, bdFlag, nSub)
%   parts = PARTITIONMESH3D(node, elem, bdFlag, nSub, 'overlap', delta)
%
%   nSub: scalar → slab partition (along x)
%         [nx, ny, nz] → 3D grid partition

p = inputParser;
p.addParameter('overlap', 0, @(x) isnumeric(x) && x >= 0);
p.parse(varargin{:});
delta = p.Results.overlap;

if isscalar(nSub)
    nx = nSub;  ny = 1;  nz = 1;  mode = 'slab';
else
    nx = nSub(1);  ny = nSub(2);  nz = nSub(3);  mode = 'grid';
end
nTotal = nx * ny * nz;

xMin = min(node(:,1));  xMax = max(node(:,1));
yMin = min(node(:,2));  yMax = max(node(:,2));
zMin = min(node(:,3));  zMax = max(node(:,3));
NT = size(elem, 1);
Hx = (xMax - xMin) / nx;
Hy = (yMax - yMin) / ny;
Hz = (zMax - zMin) / nz;

% Element centroids
xC = (node(elem(:,1), 1) + node(elem(:,2), 1) + ...
      node(elem(:,3), 1) + node(elem(:,4), 1)) / 4;
yC = (node(elem(:,1), 2) + node(elem(:,2), 2) + ...
      node(elem(:,3), 2) + node(elem(:,4), 2)) / 4;
zC = (node(elem(:,1), 3) + node(elem(:,2), 3) + ...
      node(elem(:,3), 3) + node(elem(:,4), 3)) / 4;

bdNodes = getBoundaryNodes3D(elem, bdFlag);

% Linear index: sub2lin(i,j,k) = (k-1)*nx*ny + (j-1)*nx + i
sub2lin = @(i,j,k) (k-1)*nx*ny + (j-1)*nx + i;

% ---- Non-overlapping base --------------------------------------------------
baseElem = cell(nTotal, 1);
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            s = sub2lin(i, j, k);
            xL = xMin + (i-1)*Hx;  xR = xMin + i*Hx;
            yL = yMin + (j-1)*Hy;  yR = yMin + j*Hy;
            zL = zMin + (k-1)*Hz;  zR = zMin + k*Hz;
            if i == 1,  xL = xMin - 1e-12; end
            if i == nx, xR = xMax + 1e-12; end
            if j == 1,  yL = yMin - 1e-12; end
            if j == ny, yR = yMax + 1e-12; end
            if k == 1,  zL = zMin - 1e-12; end
            if k == nz, zR = zMax + 1e-12; end
            baseElem{s} = find(xC >= xL & xC < xR & yC >= yL & yC < yR & zC >= zL & zC < zR);
        end
    end
end

% ---- Overlapping -----------------------------------------------------------
if delta > 0
    extElem = cell(nTotal, 1);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                s = sub2lin(i, j, k);
                xL = max(xMin, xMin + (i-1)*Hx - delta);
                xR = min(xMax, xMin + i*Hx + delta);
                yL = max(yMin, yMin + (j-1)*Hy - delta);
                yR = min(yMax, yMin + j*Hy + delta);
                zL = max(zMin, zMin + (k-1)*Hz - delta);
                zR = min(zMax, zMin + k*Hz + delta);
                extElem{s} = find(xC >= xL & xC <= xR & yC >= yL & yC <= yR & zC >= zL & zC <= zR);
            end
        end
    end
    useElem = extElem;
else
    useElem = baseElem;
end

% ---- Build partition structs -----------------------------------------------
faceVertMapping = {[2 3 4], [1 3 4], [1 2 4], [1 2 3]};
elemIdsByNode = repmat((1:NT)', 1, 4);
nodeElemInc = sparse(elem(:), elemIdsByNode(:), true, size(node, 1), NT);
nodeIncidentCount = full(sum(nodeElemInc, 2));

parts = struct();
for s = 1:nTotal
    eIdx = useElem{s};
    parts(s).elemIdx = eIdx;
    parts(s).nodeIdx = unique(elem(eIdx, :));

    nLocal = length(parts(s).nodeIdx);
    g2l = zeros(size(node, 1), 1);
    g2l(parts(s).nodeIdx) = (1:nLocal)';
    parts(s).localNode = node(parts(s).nodeIdx, :);
    parts(s).localElem = g2l(elem(eIdx, :));
    parts(s).global2local = g2l;

    isBd = ismember(parts(s).nodeIdx, bdNodes);
    parts(s).bdIdx   = find(isBd);
    parts(s).freeIdx = find(~isBd);

    localIncidentCount = full(sum(nodeElemInc(parts(s).nodeIdx, eIdx), 2));
    isInterior = localIncidentCount == nodeIncidentCount(parts(s).nodeIdx);
    parts(s).interiorNodeIdx = parts(s).nodeIdx(isInterior);
    parts(s).boundaryNodeIdx = parts(s).nodeIdx(~isInterior);
end

% ---- Interface detection ---------------------------------------------------
for s = 1:nTotal
    parts(s).nIfaces = 0;
    parts(s).ifaceNodes = {};
    parts(s).ifaceEdges = {};
    parts(s).ifaceNeighbor = [];
    parts(s).ifaceNodeMap = {};
end

if delta == 0
    % Generate neighbor pairs: all (i,j,k) ↔ (i±1,j,k), (i,j±1,k), (i,j,k±1)
    neighborPairs = zeros(0, 2);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                s = sub2lin(i, j, k);
                if i < nx, neighborPairs(end+1, :) = [s, sub2lin(i+1, j, k)]; end %#ok<AGROW>
                if j < ny, neighborPairs(end+1, :) = [s, sub2lin(i, j+1, k)]; end %#ok<AGROW>
                if k < nz, neighborPairs(end+1, :) = [s, sub2lin(i, j, k+1)]; end %#ok<AGROW>
            end
        end
    end

    for p = 1:size(neighborPairs, 1)
        sL = neighborPairs(p, 1);  sR = neighborPairs(p, 2);
        sharedNodes = intersect(parts(sL).nodeIdx, parts(sR).nodeIdx);
        if isempty(sharedNodes), continue; end

        for side = 1:2
            if side == 1, sCur = sL; sNbr = sR;
            else,         sCur = sR; sNbr = sL; end

            g2l = parts(sCur).global2local;
            locIfcNodes = g2l(sharedNodes);
            locIfcNodes = locIfcNodes(locIfcNodes > 0);
            locIfcSet = false(size(parts(sCur).localNode, 1), 1);
            locIfcSet(locIfcNodes) = true;

            locElem = parts(sCur).localElem;
            nLocElem = size(locElem, 1);
            locFaces = [locElem(:, faceVertMapping{1});
                        locElem(:, faceVertMapping{2});
                        locElem(:, faceVertMapping{3});
                        locElem(:, faceVertMapping{4})];
            elemFaceRows = repmat((1:nLocElem)', 4, 1);
            localFaceIds = [ones(nLocElem, 1);
                            2 * ones(nLocElem, 1);
                            3 * ones(nLocElem, 1);
                            4 * ones(nLocElem, 1)];
            isIfaceFace = all(locIfcSet(locFaces), 2);
            ifaceFaces = [elemFaceRows(isIfaceFace), localFaceIds(isIfaceFace)];

            g2lNbr = parts(sNbr).global2local;
            locToNbr = zeros(length(locIfcNodes), 1);
            for i = 1:length(locIfcNodes)
                glb = parts(sCur).nodeIdx(locIfcNodes(i));
                locToNbr(i) = g2lNbr(glb);
            end

            nIf = parts(sCur).nIfaces + 1;
            parts(sCur).nIfaces = nIf;
            parts(sCur).ifaceNodes{nIf} = locIfcNodes;
            parts(sCur).ifaceEdges{nIf} = ifaceFaces;
            parts(sCur).ifaceNeighbor(nIf) = sNbr;
            parts(sCur).ifaceNodeMap{nIf} = locToNbr;
        end
    end
end
end
