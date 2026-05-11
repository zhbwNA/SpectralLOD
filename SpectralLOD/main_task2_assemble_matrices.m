function fem = main_task2_assemble_matrices()
%% TASK 2: Assemble P1 FEM matrices
%
% Assemble:
%     K      = stiffness matrix
%     M      = volume mass matrix
%     Mbd    = boundary mass matrix on physical boundary partial Omega
%     Mskel  = edge mass matrix on internal skeleton Gamma
%     Msigma = edge mass matrix on Sigma = Gamma union partial Omega
%
% Also assemble the impedance Helmholtz matrix
%
%     Ahelm = K - k^2 M - 1i*k*Mbd.
%
% This script calls Task 1:
%
%     [node,elem,mesh,dd] = main_task1_mesh_partition();

clear; clc;

%% Parameters

k = 20;       % wave number, adjustable later

%% Generate mesh and decomposition from Task 1

[node,elem,mesh,dd] = main_task1_mesh_partition();

%% Assemble volume matrices

[K,M,area] = assembleP1VolumeMatrices(node,elem);

%% Assemble edge mass matrices

Mbd    = assembleP1EdgeMass(node,mesh.edge,dd.bdEdgeIdx);
Mskel  = assembleP1EdgeMass(node,mesh.edge,dd.skeletonEdgeIdx);
Msigma = assembleP1EdgeMass(node,mesh.edge,dd.sigmaEdgeIdx);

%% Helmholtz impedance matrix

Ahelm = K - k^2*M - 1i*k*Mbd;

%% Store outputs

fem.node = node;
fem.elem = elem;
fem.mesh = mesh;
fem.dd = dd;

fem.K = K;
fem.M = M;
fem.Mbd = Mbd;
fem.Mskel = Mskel;
fem.Msigma = Msigma;
fem.Ahelm = Ahelm;
fem.area = area;
fem.k = k;

%% Print and test

fprintf('\n===== FEM matrix information =====\n');
fprintf('Number of DOFs                : %d\n', size(node,1));
fprintf('Total area from elements      : %.16e\n', sum(area));
fprintf('nnz(K)                        : %d\n', nnz(K));
fprintf('nnz(M)                        : %d\n', nnz(M));
fprintf('nnz(Mbd)                      : %d\n', nnz(Mbd));
fprintf('nnz(Mskel)                    : %d\n', nnz(Mskel));
fprintf('nnz(Msigma)                   : %d\n', nnz(Msigma));

checkAssembledMatrices(node,K,M,Mbd,Mskel,Msigma,area);

end


%% ------------------------------------------------------------------------
function [K,M,area] = assembleP1VolumeMatrices(node,elem)
%%ASSEMBLEP1VOLUMEMATRICES Assemble stiffness and mass matrices for P1 FEM.
%
% Input:
%     node : N x 2 node coordinates
%     elem : NT x 3 triangle connectivity
%
% Output:
%     K    : stiffness matrix
%     M    : volume mass matrix
%     area : NT x 1 element areas

N = size(node,1);
NT = size(elem,1);

K = sparse(N,N);
M = sparse(N,N);

area = zeros(NT,1);

localMassRef = [2 1 1; 1 2 1; 1 1 2]/12;

for t = 1:NT
    vid = elem(t,:);
    p = node(vid,:);

    x1 = p(1,1); y1 = p(1,2);
    x2 = p(2,1); y2 = p(2,2);
    x3 = p(3,1); y3 = p(3,2);

    % Signed double area.
    detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
    area(t) = abs(detJ)/2;

    if area(t) <= 0
        error('Element %d has nonpositive area.',t);
    end

    % Gradients of barycentric basis:
    %
    % grad phi_i = [b_i, c_i]/(2*area).
    b = [y2-y3; y3-y1; y1-y2];
    c = [x3-x2; x1-x3; x2-x1];

    localK = (b*b' + c*c')/(4*area(t));
    localM = area(t)*localMassRef;

    K(vid,vid) = K(vid,vid) + localK;
    M(vid,vid) = M(vid,vid) + localM;
end

end


%% ------------------------------------------------------------------------
function Medge = assembleP1EdgeMass(node,edge,edgeIdx)
%%ASSEMBLEP1EDGEMASS Assemble P1 edge mass matrix on selected edges.
%
% For an edge e = [i,j], the local P1 edge mass is
%
%     |e|/6 * [2 1; 1 2].
%
% Input:
%     node    : N x 2 node coordinates
%     edge    : NE x 2 edge connectivity
%     edgeIdx : selected edge indices
%
% Output:
%     Medge   : sparse N x N edge mass matrix

N = size(node,1);
Medge = sparse(N,N);

if isempty(edgeIdx)
    return;
end

selectedEdge = edge(edgeIdx,:);

for e = 1:size(selectedEdge,1)
    vid = selectedEdge(e,:);
    p1 = node(vid(1),:);
    p2 = node(vid(2),:);

    len = norm(p2-p1);

    localM = len/6*[2 1; 1 2];

    Medge(vid,vid) = Medge(vid,vid) + localM;
end

end


%% ------------------------------------------------------------------------
function checkAssembledMatrices(node,K,M,Mbd,Mskel,Msigma,area)
%%CHECKASSEMBLEDMATRICES Sanity checks for assembled matrices.

N = size(node,1);
one = ones(N,1);

fprintf('\n===== Matrix sanity checks =====\n');

% Symmetry checks.
symK = norm(K-K','fro')/max(1,norm(K,'fro'));
symM = norm(M-M','fro')/max(1,norm(M,'fro'));
symMbd = norm(Mbd-Mbd','fro')/max(1,norm(Mbd,'fro'));
symMskel = norm(Mskel-Mskel','fro')/max(1,norm(Mskel,'fro'));
symMsigma = norm(Msigma-Msigma','fro')/max(1,norm(Msigma,'fro'));

fprintf('symmetry error K              : %.3e\n', symK);
fprintf('symmetry error M              : %.3e\n', symM);
fprintf('symmetry error Mbd            : %.3e\n', symMbd);
fprintf('symmetry error Mskel          : %.3e\n', symMskel);
fprintf('symmetry error Msigma         : %.3e\n', symMsigma);

% Stiffness annihilates constants for pure Neumann volume matrix.
Kone = norm(K*one);
fprintf('||K*1||                       : %.3e\n', Kone);

% Mass of constant one equals area.
volMass = one'*M*one;
trueArea = sum(area);
fprintf('1^T M 1                       : %.16e\n', volMass);
fprintf('sum element areas             : %.16e\n', trueArea);
fprintf('volume mass error             : %.3e\n', abs(volMass-trueArea));

% Boundary mass of constant one should equal perimeter of unit square.
bdMass = one'*Mbd*one;
fprintf('1^T Mbd 1                     : %.16e\n', bdMass);
fprintf('expected perimeter            : %.16e\n', 4.0);
fprintf('boundary mass error           : %.3e\n', abs(bdMass-4.0));

% Sigma mass should contain boundary and skeleton.
sigmaMass = one'*Msigma*one;
skelMass = one'*Mskel*one;
fprintf('1^T Mskel 1                   : %.16e\n', skelMass);
fprintf('1^T Msigma 1                  : %.16e\n', sigmaMass);
fprintf('check Msigma = Mbd + Mskel    : %.3e\n', ...
    norm(Msigma-(Mbd+Mskel),'fro'));

tol = 1e-10;

if symK > tol || symM > tol || symMbd > tol || symMskel > tol || symMsigma > tol
    warning('Some assembled matrices have symmetry errors larger than tolerance.');
end

if abs(volMass-trueArea) > 1e-10
    warning('Volume mass check failed.');
end

if abs(bdMass-4.0) > 1e-10
    warning('Boundary mass check failed. This check assumes Omega=(0,1)^2.');
end

fprintf('\nSanity check completed.\n');

end