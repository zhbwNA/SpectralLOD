function applyPrecon = orasPMLHelmholtz2D(node, elem, k, parts, pml, solverMode, useParfor)
% ORASPMLHELMHOLTZ2D  Additive Schwarz/ORAS-style preconditioner with PML local solves.
%
%   applyPrecon = ORASPMLHELMHOLTZ2D(node,elem,k,parts,pml)
%
%   This first version is P1-only. Each subdomain matrix is assembled with
%   assembleHelmholtzPML2D on the local mesh. Local PML boxes are taken from
%   parts(j).coreBox and parts(j).extendedBox when present; otherwise the
%   local mesh bounding box is used.

if nargin < 6 || isempty(solverMode), solverMode = 'lu'; end
if nargin < 7 || isempty(useParfor), useParfor = false; end
if nargin < 5 || isempty(pml), pml = struct(); end
if size(elem, 2) ~= 3
    error('orasPMLHelmholtz2D:p1Only', 'PML ORAS currently supports P1 triangles only.');
end

N = size(node, 1);
nSub = length(parts);

useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(1).weightFun);
nodeWeight = zeros(N, 1);
for j = 1:nSub
    idx = unique(elem(parts(j).elemIdx, :));
    if useWeightFun
        raw = parts(j).weightFun(node(idx,1), node(idx,2));
        nodeWeight(idx) = nodeWeight(idx) + max(raw(:), 0);
    else
        nodeWeight(idx) = nodeWeight(idx) + 1;
    end
end
nodeWeight(nodeWeight == 0) = 1;

locSolvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
wgt = cell(nSub, 1);
freeLocal = cell(nSub, 1);

if useParfor
    parfor j = 1:nSub
        [locSolvers{j}, gIdx{j}, wgt{j}, freeLocal{j}] = setupPMLSubdomain( ...
            j, node, elem, parts, k, pml, solverMode, nodeWeight, useWeightFun);
    end
else
    for j = 1:nSub
        [locSolvers{j}, gIdx{j}, wgt{j}, freeLocal{j}] = setupPMLSubdomain( ...
            j, node, elem, parts, k, pml, solverMode, nodeWeight, useWeightFun);
    end
end

    function x = applyImpl(r)
        x = zeros(N, 1);
        for s = 1:nSub
            rloc = r(gIdx{s});
            zloc = zeros(size(rloc));
            fd = freeLocal{s};
            if isempty(fd)
                continue;
            end
            if strcmpi(solverMode, 'direct')
                zloc(fd) = locSolvers{s} \ rloc(fd);
            else
                S = locSolvers{s};
                zfd = zeros(length(fd), 1);
                zfd(S{4}) = S{2} \ (S{1} \ rloc(fd(S{3})));
                zloc(fd) = zfd;
            end
            x(gIdx{s}) = x(gIdx{s}) + zloc .* wgt{s};
        end
    end

applyPrecon = @applyImpl;
end


function [solver, gIdx, wgt, freeDof] = setupPMLSubdomain( ...
    j, node, elem, parts, k, pml, solverMode, nodeWeight, useWeightFun)

eIdx = parts(j).elemIdx;
gIdx = unique(elem(eIdx, :));
g2l = zeros(size(node, 1), 1);
g2l(gIdx) = (1:length(gIdx)).';

localNode = node(gIdx, :);
localElem = g2l(elem(eIdx, :));
localPML = localPMLStruct(parts(j), localNode, k, pml);

[A_loc, ~] = assembleHelmholtzPML2D(localNode, localElem, k, localPML, []);
bdDof = localBoundaryNodes2D(localElem);
freeDof = setdiff((1:size(localNode,1)).', bdDof(:));
A_free = A_loc(freeDof, freeDof);

if strcmpi(solverMode, 'direct')
    solver = A_free;
else
    [L, U, p, q] = lu(A_free, 'vector');
    solver = {L, U, p(:), q(:)};
end


function bd = localBoundaryNodes2D(elem)
edges = [elem(:, [1 2]); elem(:, [2 3]); elem(:, [3 1])];
edgesS = sort(edges, 2);
[uEdges, ~, ic] = unique(edgesS, 'rows'); %#ok<ASGLU>
counts = accumarray(ic, 1);
bdEdges = edges(counts(ic) == 1, :);
bd = unique(bdEdges(:));
end

if useWeightFun
    raw = max(parts(j).weightFun(node(gIdx,1), node(gIdx,2)), 0);
    wgt = raw(:) ./ nodeWeight(gIdx);
else
    wgt = 1 ./ nodeWeight(gIdx);
end
end


function localPML = localPMLStruct(part, localNode, k, pml)
localPML = pml;
if ~isfield(localPML, 'sigmaMax') || isempty(localPML.sigmaMax)
    localPML.sigmaMax = k;
end
if ~isfield(localPML, 'sigmaOrder') || isempty(localPML.sigmaOrder)
    localPML.sigmaOrder = 2;
end
if isfield(part, 'coreBox') && isfield(part, 'extendedBox')
    localPML.physicalBox = part.coreBox;
    localPML.pmlBox = part.extendedBox;
else
    box = [min(localNode(:,1)), max(localNode(:,1)), min(localNode(:,2)), max(localNode(:,2))];
    localPML.physicalBox = box;
    localPML.pmlBox = box;
end
end
