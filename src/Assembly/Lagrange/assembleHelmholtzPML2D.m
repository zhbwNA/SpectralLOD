function [A, b, freeDof, bdDof, info] = assembleHelmholtzPML2D(node, elem, k, pml, f)
% ASSEMBLEHELMHOLTZPML2D  Assemble P1 Helmholtz operator with Cartesian PML.
%
%   [A,b,freeDof,bdDof] = ASSEMBLEHELMHOLTZPML2D(node,elem,k,pml,f)
%
%   The bilinear form is
%       int A_pml grad u . grad v - k^2 int b_pml u v,
%   with A_pml = diag(s2/s1, s1/s2), b_pml = s1*s2.
%
%   Homogeneous Dirichlet data are imposed by the caller on bdDof, usually
%   by solving A(freeDof,freeDof) u = b(freeDof).

if nargin < 5
    f = [];
end
if size(elem, 2) ~= 3
    error('assembleHelmholtzPML2D:p1Only', 'PML assembly currently supports P1 triangles only.');
end
if nargin < 4 || isempty(pml)
    pml = struct();
end
if ~isfield(pml, 'physicalBox') || isempty(pml.physicalBox)
    pml.physicalBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
if ~isfield(pml, 'pmlBox') || isempty(pml.pmlBox)
    pml.pmlBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
if ~isfield(pml, 'quadOrder') || isempty(pml.quadOrder)
    pml.quadOrder = 4;
end

N = size(node, 1);
NT = size(elem, 1);
[lambda, weight] = quadtriangle(pml.quadOrder);
nQuad = length(weight);

x1 = node(elem(:,1), 1);  y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);  y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);  y3 = node(elem(:,3), 2);

area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area = abs(area2) / 2;

g1x = (y2 - y3) ./ area2;  g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2;  g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2;  g3y = (x2 - x1) ./ area2;
gx = [g1x, g2x, g3x];
gy = [g1y, g2y, g3y];

[aa, bb] = ndgrid(1:3, 1:3);
aa = aa(:).';  bb = bb(:).';
S = zeros(NT, 9);
rhs = zeros(NT, 3);

for q = 1:nQuad
    lq = lambda(q, :);
    xq = lq(1)*x1 + lq(2)*x2 + lq(3)*x3;
    yq = lq(1)*y1 + lq(2)*y2 + lq(3)*y3;
    [a11, a22, bcoef] = pmlCoefficients2D(xq, yq, k, pml);
    wphys = 2 * weight(q) * area;

    stiff = a11 .* gx(:,aa) .* gx(:,bb) + a22 .* gy(:,aa) .* gy(:,bb);
    mass = bcoef .* lq(aa) .* lq(bb);
    S = S + wphys .* (stiff - k^2 * mass);

    if nargout > 1 && ~isempty(f)
        fq = evalLoad(f, xq, yq);
        rhs = rhs + wphys .* fq .* lq;
    end
end

ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(S, [], 1);
A = sparse(ii, jj, ss, N, N);

if nargout > 1
    b = accumarray(reshape(elem, [], 1), reshape(rhs, [], 1), [N, 1], @sum, 0);
end

bdDof = outerBoundaryNodes2D(node, pml.pmlBox);
freeDof = setdiff((1:N).', bdDof(:));

info = struct();
info.physicalBox = pml.physicalBox;
info.pmlBox = pml.pmlBox;
info.dirichletDof = bdDof;
info.freeDof = freeDof;
end


function val = evalLoad(f, x, y)
if isnumeric(f)
    val = f * ones(size(x));
else
    val = f(x, y);
end
end


function bd = outerBoundaryNodes2D(node, box)
tol = 100 * eps(max(1, max(abs(box))));
onX = abs(node(:,1) - box(1)) <= tol | abs(node(:,1) - box(2)) <= tol;
onY = abs(node(:,2) - box(3)) <= tol | abs(node(:,2) - box(4)) <= tol;
bd = find(onX | onY);
end
