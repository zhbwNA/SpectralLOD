function [node,elem,mesh,dd] = main_task1_mesh_partition()
%% TASK 1: Mesh generation and nonoverlapping subdomain partition
%
% Domain:
%     Omega = (0,1)^2
%
% Subdomain decomposition:
%     square subdomains of size approximately H0.
%
% Fine mesh:
%     Each subdomain is divided into nFinePerSub intervals in each direction.
%
% Output:
%     node, elem : iFEM-style mesh arrays
%     mesh       : mesh auxiliary data
%     dd         : domain-decomposition data
%
% Author: ChatGPT
% Style: close to Chen Long's iFEM convention

clear; clc;

%% Adjustable parameters

domain = [0, 1, 0, 1];   % [xmin xmax ymin ymax]

H0 = 1/4;                % target subdomain size

nFinePerSub = 8;         % number of fine intervals per subdomain edge

plotFlag = true;

%% Build mesh and domain decomposition structure

[node,elem,mesh,dd] = squareDDMesh(domain,H0,nFinePerSub);

%% Print basic information

fprintf('\n===== Mesh information =====\n');
fprintf('Number of nodes             : %d\n', size(node,1));
fprintf('Number of elements          : %d\n', size(elem,1));
fprintf('Number of edges             : %d\n', size(mesh.edge,1));

fprintf('\n===== Subdomain information =====\n');
fprintf('Number of subdomains in x   : %d\n', dd.nSubX);
fprintf('Number of subdomains in y   : %d\n', dd.nSubY);
fprintf('Total number of subdomains  : %d\n', dd.nSub);
fprintf('Actual subdomain size Hx    : %.6g\n', dd.Hx);
fprintf('Actual subdomain size Hy    : %.6g\n', dd.Hy);

fprintf('\n===== DOF classification =====\n');
fprintf('Physical boundary nodes     : %d\n', numel(dd.bdNode));
fprintf('Skeleton nodes              : %d\n', numel(dd.skeletonNode));
fprintf('Full skeleton Sigma nodes   : %d\n', numel(dd.sigmaNode));
fprintf('Physical boundary edges     : %d\n', numel(dd.bdEdgeIdx));
fprintf('Internal skeleton edges     : %d\n', numel(dd.skeletonEdgeIdx));
fprintf('Full skeleton Sigma edges   : %d\n', numel(dd.sigmaEdgeIdx));

%% Sanity check

checkDDMesh(node,elem,mesh,dd);

%% Plot

if plotFlag
    plotDDMesh(node,elem,mesh,dd);
end

end


%% ------------------------------------------------------------------------
function [node,elem,mesh,dd] = squareDDMesh(domain,H0,nFinePerSub)
%%SQUAREDDMESH Generate a uniform triangular mesh and subdomain partition.
%
% Input:
%     domain       = [xmin xmax ymin ymax]
%     H0           = target subdomain size
%     nFinePerSub  = number of fine intervals per subdomain edge
%
% Output:
%     node, elem   : iFEM-style mesh arrays
%     mesh         : auxiliary mesh structure
%     dd           : domain decomposition structure

xmin = domain(1); xmax = domain(2);
ymin = domain(3); ymax = domain(4);

Lx = xmax - xmin;
Ly = ymax - ymin;

% Number of subdomains.
% We round to the nearest integer so that subdomains exactly fit the domain.
nSubX = round(Lx/H0);
nSubY = round(Ly/H0);

if nSubX < 1 || nSubY < 1
    error('H0 is too large for the given domain.');
end

Hx = Lx/nSubX;
Hy = Ly/nSubY;

% Fine mesh resolution.
Nx = nSubX*nFinePerSub;
Ny = nSubY*nFinePerSub;

hx = Lx/Nx;
hy = Ly/Ny;

% Node coordinates.
[xgrid,ygrid] = meshgrid(linspace(xmin,xmax,Nx+1), ...
                         linspace(ymin,ymax,Ny+1));
node = [xgrid(:), ygrid(:)];

% Node indexing: id(i,j), i=0,...,Nx, j=0,...,Ny.
id = @(i,j) j*(Nx+1) + i + 1;

% Triangulation: split each rectangle into two triangles.
elem = zeros(2*Nx*Ny,3);
t = 0;
for j = 0:Ny-1
    for i = 0:Nx-1
        n1 = id(i,j);
        n2 = id(i+1,j);
        n3 = id(i,j+1);
        n4 = id(i+1,j+1);

        % Counterclockwise orientation.
        t = t + 1;
        elem(t,:) = [n1 n2 n4];

        t = t + 1;
        elem(t,:) = [n1 n4 n3];
    end
end

%% Build auxiliary mesh structure

mesh = auxstructure(node,elem);

mesh.domain = domain;
mesh.Nx = Nx;
mesh.Ny = Ny;
mesh.hx = hx;
mesh.hy = hy;

%% Build domain decomposition structure

dd = partitionSubdomains(node,elem,mesh,domain,H0,nSubX,nSubY,Hx,Hy);

end


%% ------------------------------------------------------------------------
function mesh = auxstructure(node,elem)
%%AUXSTRUCTURE Build edge, elem2edge, bdEdge data.
%
% This follows the usual iFEM-style auxiliary structure.

NT = size(elem,1);

% Local edges opposite to local vertices:
% edge 1: [2 3], edge 2: [3 1], edge 3: [1 2].
allEdge = [elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])];
allEdge = sort(allEdge,2);

[edge,~,idx] = unique(allEdge,'rows');

elem2edge = reshape(idx,NT,3);

edgeCount = accumarray(idx,1);
bdEdgeIdx = find(edgeCount == 1);
bdEdge = edge(bdEdgeIdx,:);

mesh.edge = edge;
mesh.elem2edge = elem2edge;
mesh.bdEdgeIdx = bdEdgeIdx;
mesh.bdEdge = bdEdge;

mesh.bdNode = unique(bdEdge(:));

end


%% ------------------------------------------------------------------------
function dd = partitionSubdomains(node,elem,mesh,domain,H0,nSubX,nSubY,Hx,Hy)
%%PARTITIONSUBDOMAINS Classify elements, nodes, edges by subdomains.

xmin = domain(1); xmax = domain(2);
ymin = domain(3); ymax = domain(4);

tol = 1e-12;

nSub = nSubX*nSubY;
NT = size(elem,1);

%% Assign each element to one subdomain by centroid

centroid = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;

ix = floor((centroid(:,1)-xmin)/Hx) + 1;
iy = floor((centroid(:,2)-ymin)/Hy) + 1;

ix = min(max(ix,1),nSubX);
iy = min(max(iy,1),nSubY);

elemSubId = ix + (iy-1)*nSubX;

subElem = cell(nSub,1);
for s = 1:nSub
    subElem{s} = find(elemSubId == s);
end

%% Physical boundary nodes

x = node(:,1);
y = node(:,2);

isBdNode = abs(x-xmin) < tol | abs(x-xmax) < tol | ...
           abs(y-ymin) < tol | abs(y-ymax) < tol;
bdNode = find(isBdNode);

%% Internal skeleton nodes

isInternalVertical = false(size(x));
for p = 1:nSubX-1
    xp = xmin + p*Hx;
    isInternalVertical = isInternalVertical | abs(x-xp) < tol;
end

isInternalHorizontal = false(size(y));
for q = 1:nSubY-1
    yq = ymin + q*Hy;
    isInternalHorizontal = isInternalHorizontal | abs(y-yq) < tol;
end

isSkeletonNode = (isInternalVertical | isInternalHorizontal) & ~isBdNode;
skeletonNode = find(isSkeletonNode);

% Full skeleton Sigma = internal skeleton plus physical boundary.
sigmaNode = union(skeletonNode,bdNode);

%% Edge classification

edge = mesh.edge;
e1 = edge(:,1);
e2 = edge(:,2);

x1 = node(e1,1); y1 = node(e1,2);
x2 = node(e2,1); y2 = node(e2,2);

% Physical boundary edges are already available from auxstructure.
bdEdgeIdx = mesh.bdEdgeIdx;

% Internal skeleton edges: edges lying exactly on internal subdomain lines.
isSameX = abs(x1-x2) < tol;
isSameY = abs(y1-y2) < tol;

isOnInternalVerticalEdge = false(size(e1));
for p = 1:nSubX-1
    xp = xmin + p*Hx;
    isOnInternalVerticalEdge = isOnInternalVerticalEdge | ...
        (isSameX & abs(x1-xp) < tol & abs(x2-xp) < tol);
end

isOnInternalHorizontalEdge = false(size(e1));
for q = 1:nSubY-1
    yq = ymin + q*Hy;
    isOnInternalHorizontalEdge = isOnInternalHorizontalEdge | ...
        (isSameY & abs(y1-yq) < tol & abs(y2-yq) < tol);
end

isSkeletonEdge = isOnInternalVerticalEdge | isOnInternalHorizontalEdge;

% Exclude physical boundary edges.
isBdEdge = false(size(e1));
isBdEdge(bdEdgeIdx) = true;

skeletonEdgeIdx = find(isSkeletonEdge & ~isBdEdge);

sigmaEdgeIdx = union(bdEdgeIdx,skeletonEdgeIdx);

%% Per-subdomain node and DOF classification

subNode = cell(nSub,1);
subBdNode = cell(nSub,1);
subIntNode = cell(nSub,1);
subBdEdgeIdx = cell(nSub,1);
subSigmaNode = cell(nSub,1);

for sy = 1:nSubY
    for sx = 1:nSubX
        s = sx + (sy-1)*nSubX;

        x0 = xmin + (sx-1)*Hx;
        x1s = xmin + sx*Hx;
        y0 = ymin + (sy-1)*Hy;
        y1s = ymin + sy*Hy;

        elems = subElem{s};
        nodes = unique(elem(elems,:));

        xs = node(nodes,1);
        ys = node(nodes,2);

        isOnSubBd = abs(xs-x0) < tol | abs(xs-x1s) < tol | ...
                    abs(ys-y0) < tol | abs(ys-y1s) < tol;

        subNode{s} = nodes;
        subBdNode{s} = nodes(isOnSubBd);
        subIntNode{s} = setdiff(nodes,subBdNode{s});

        % Subdomain boundary edges.
        midx = 0.5*(x1+x2);
        midy = 0.5*(y1+y2);

        edgeInBox = midx >= x0-tol & midx <= x1s+tol & ...
                    midy >= y0-tol & midy <= y1s+tol;

        edgeOnSubBd = edgeInBox & ...
            (abs(midx-x0) < tol | abs(midx-x1s) < tol | ...
             abs(midy-y0) < tol | abs(midy-y1s) < tol);

        subBdEdgeIdx{s} = find(edgeOnSubBd);

        % Local full skeleton nodes for this subdomain:
        % all nodes on partial Omega_i.
        subSigmaNode{s} = subBdNode{s};
    end
end

%% Store data

dd.domain = domain;

dd.H0_input = H0;
dd.Hx = Hx;
dd.Hy = Hy;

dd.nSubX = nSubX;
dd.nSubY = nSubY;
dd.nSub = nSub;

dd.elemSubId = elemSubId;
dd.subElem = subElem;

dd.bdNode = bdNode;
dd.skeletonNode = skeletonNode;
dd.sigmaNode = sigmaNode;

dd.bdEdgeIdx = bdEdgeIdx;
dd.skeletonEdgeIdx = skeletonEdgeIdx;
dd.sigmaEdgeIdx = sigmaEdgeIdx;

dd.subNode = subNode;
dd.subBdNode = subBdNode;
dd.subIntNode = subIntNode;
dd.subBdEdgeIdx = subBdEdgeIdx;
dd.subSigmaNode = subSigmaNode;

end


%% ------------------------------------------------------------------------
function checkDDMesh(node,elem,mesh,dd)
%%CHECKDDMESH Basic consistency checks.

N = size(node,1);
NT = size(elem,1);

% Every element should belong to exactly one subdomain.
assigned = zeros(NT,1);
for s = 1:dd.nSub
    assigned(dd.subElem{s}) = assigned(dd.subElem{s}) + 1;
end

if any(assigned ~= 1)
    error('Some elements are not assigned to exactly one subdomain.');
end

% Subdomain interior and boundary nodes should be disjoint.
for s = 1:dd.nSub
    if ~isempty(intersect(dd.subIntNode{s},dd.subBdNode{s}))
        error('Subdomain %d has overlapping interior and boundary DOFs.',s);
    end
end

% Skeleton nodes and physical boundary nodes are disjoint by definition.
if ~isempty(intersect(dd.skeletonNode,dd.bdNode))
    error('Internal skeleton nodes overlap physical boundary nodes.');
end

% Sigma nodes should be the union.
sigmaCheck = union(dd.skeletonNode,dd.bdNode);
if numel(setdiff(sigmaCheck,dd.sigmaNode)) > 0 || ...
   numel(setdiff(dd.sigmaNode,sigmaCheck)) > 0
    error('sigmaNode is not union(skeletonNode,bdNode).');
end

% Basic index range check.
if max(elem(:)) > N || min(elem(:)) < 1
    error('elem contains invalid node indices.');
end

if max(mesh.edge(:)) > N || min(mesh.edge(:)) < 1
    error('edge contains invalid node indices.');
end

fprintf('\nSanity check passed.\n');

end


%% ------------------------------------------------------------------------
function plotDDMesh(node,elem,mesh,dd)
%%PLOTDDMESH Plot mesh, physical boundary, and skeleton.

figure;
triplot(elem,node(:,1),node(:,2),'Color',[0.75 0.75 0.75]);
axis equal tight;
hold on;

% Plot internal skeleton edges.
edge = mesh.edge;

skelEdge = edge(dd.skeletonEdgeIdx,:);
for k = 1:size(skelEdge,1)
    p = skelEdge(k,:);
    plot(node(p,1),node(p,2),'r-','LineWidth',1.8);
end

% Plot physical boundary edges.
bdEdge = edge(dd.bdEdgeIdx,:);
for k = 1:size(bdEdge,1)
    p = bdEdge(k,:);
    plot(node(p,1),node(p,2),'k-','LineWidth',1.8);
end

% Plot subdomain numbers.
for s = 1:dd.nSub
    elems = dd.subElem{s};
    c = mean((node(elem(elems,1),:) + node(elem(elems,2),:) + node(elem(elems,3),:))/3,1);
    text(c(1),c(2),num2str(s), ...
        'HorizontalAlignment','center', ...
        'FontSize',12, ...
        'FontWeight','bold', ...
        'Color','b');
end

title('Task 1: fine mesh, physical boundary, and internal skeleton');
legend({'fine mesh','skeleton \Gamma','physical boundary \partial\Omega'}, ...
       'Location','bestoutside');

hold off;

end