% VERIFY_NED2_3D  Convergence test for NE_2 in 3D.
%   Element-local face DOFs (non-conforming, simpler).
%   20 DOFs/tet: 2 per edge (shared) + 8 interior face-bubbles (element-local).

fprintf('========== 3D NE_2 Convergence Study ==========\n\n');
fprintf('NOTE: element-local face DOFs (non-conforming approximation)\n\n');

u_ex = @(x,y,z) deal(zeros(size(x)),zeros(size(x)),x.*(1-x).*y.*(1-y));
curl_ex = @(x,y,z) deal(x.*(1-x).*(1-2*y),-y.*(1-y).*(1-2*x),zeros(size(x)));
f_rhs = @(x,y,z) deal(zeros(size(x)),zeros(size(x)),2*y.*(1-y)+2*x.*(1-x)+x.*(1-x).*y.*(1-y));

nRefine = 2;
fmt = '%-6s  %-8s  %-12s  %-8s  %-12s  %-8s\n';
fprintf(fmt, 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-',1,66));

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = cubemesh([0,1,0,1,0,1], hk);
    [~, eidx, es] = edgeMesh3D(el);
    [~, fidx] = faceMesh3D(nd, el);
    NE = max(eidx(:)); NT = size(el,1);
    Ntot = 2*NE + 8*NT;                   % 2 per edge (shared) + 8 per element

    [gIdx, gSign] = ned2dofElem3D(el, eidx, es, fidx, NE, NT);

    A = ned2curl3D(nd, el, gIdx, gSign, Ntot);
    M = ned2mass3D(nd, el, gIdx, gSign, Ntot);
    K = A + M;
    b = ned2rhs3D(nd, el, f_rhs, gIdx, gSign, Ntot);

    % Boundary edges only (face DOFs are element-local, not constrained)
    faceEdges = {[4,6,5],[2,6,3],[1,5,3],[1,4,2]};
    bdE = [];
    for f = 1:4
        if any(bd(:,f)==1)
            bdE = [bdE; eidx(bd(:,f)==1, faceEdges{f})];
        end
    end
    bdE=unique(bdE);
    bdDOFs = [2*(bdE-1)+1; 2*(bdE-1)+2];
    freeDOFs = setdiff(1:Ntot, bdDOFs)';

    uh = zeros(Ntot,1);
    uh(freeDOFs) = K(freeDOFs, freeDOFs) \ b(freeDOFs);

    [eL2, eHc] = ned2err3D(nd, el, uh, u_ex, curl_ex, gIdx, gSign);

    if k > 1
        rL2=log(eL2/eL2p)/log(hk/hp); rHc=log(eHc/eHcp)/log(hk/hp);
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',Ntot), sprintf('%.4e',eL2), sprintf('%.2f',rL2), sprintf('%.4e',eHc), sprintf('%.2f',rHc));
    else
        fprintf(fmt, sprintf('%.4f',hk), sprintf('%d',Ntot), sprintf('%.4e',eL2), '-', sprintf('%.4e',eHc), '-');
    end
    eL2p=eL2; eHcp=eHc; hp=hk;
end
fprintf('\nExpected: NE_2 L2~O(h^2), H(curl)~O(h^2)\n');
fprintf('========== Done ==========\n');


% ===== DOF indexing (element-local face DOFs) =============================
function [gIdx, gSign] = ned2dofElem3D(elem, eidx, es, fidx, NE, NT)
nLocal = 20;
gIdx = zeros(NT, nLocal); gSign = zeros(NT, nLocal);
edges = [1 2;1 3;1 4;2 3;2 4;3 4];
for k = 1:6
    gIdx(:,2*(k-1)+1)=2*(eidx(:,k)-1)+1; gSign(:,2*(k-1)+1)=es(:,k);  % odd
    gIdx(:,2*(k-1)+2)=2*(eidx(:,k)-1)+2; gSign(:,2*(k-1)+2)=1;         % even
end
% Face DOFs: element-local (NOT shared)
for t = 1:NT
    for f = 1:4
        d0 = 12 + 2*(f-1) + 1; d1 = d0 + 1;
        gIdx(t, d0) = 2*NE + 8*(t-1) + 2*(f-1) + 1;
        gIdx(t, d1) = 2*NE + 8*(t-1) + 2*(f-1) + 2;
    end
end
gSign(:,13:20) = 1;
end

% ===== Geometry ============================================================
function [G, vol] = tetGeom(node, elem)
NT=size(elem,1); v1=node(elem(:,1),:);v2=node(elem(:,2),:);v3=node(elem(:,3),:);v4=node(elem(:,4),:);
e12=v2-v1;e13=v3-v1;e14=v4-v1;
detJ=e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2))+e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3))+e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
vol=abs(detJ)/6; invJ=1./detJ;
g2=cross(e13,e14).*invJ; g3=cross(e14,e12).*invJ; g4=cross(e12,e13).*invJ; g1=-(g2+g3+g4);
G=cell(4,1); G{1}=g1; G{2}=g2; G{3}=g3; G{4}=g4;
end

% ===== Basis evaluation ===================================================
function [phix,phiy,phiz,curlV] = ned2basis3D(l, G)
NT=size(G{1},1); nLocal=20;
phix=zeros(NT,nLocal); phiy=zeros(NT,nLocal); phiz=zeros(NT,nLocal);
curlV=zeros(NT,nLocal,3);
edges=[1 2;1 3;1 4;2 3;2 4;3 4];
for k=1:6
    i=edges(k,1);j=edges(k,2);gi=G{i};gj=G{j};li=l(i);lj=l(j);d0=2*(k-1)+1;d1=d0+1;
    phix(:,d0)=li*gj(:,1)-lj*gi(:,1);phiy(:,d0)=li*gj(:,2)-lj*gi(:,2);phiz(:,d0)=li*gj(:,3)-lj*gi(:,3);
    curlV(:,d0,1)=2*(gi(:,2).*gj(:,3)-gi(:,3).*gj(:,2));
    curlV(:,d0,2)=2*(gi(:,3).*gj(:,1)-gi(:,1).*gj(:,3));
    curlV(:,d0,3)=2*(gi(:,1).*gj(:,2)-gi(:,2).*gj(:,1));
    cij=li-lj;dgi=gi-gj;
    phix(:,d1)=cij.*phix(:,d0);phiy(:,d1)=cij.*phiy(:,d0);phiz(:,d1)=cij.*phiz(:,d0);
    curlV(:,d1,1)=(dgi(:,2).*phiz(:,d0)-dgi(:,3).*phiy(:,d0))+cij.*curlV(:,d0,1);
    curlV(:,d1,2)=(dgi(:,3).*phix(:,d0)-dgi(:,1).*phiz(:,d0))+cij.*curlV(:,d0,2);
    curlV(:,d1,3)=(dgi(:,1).*phiy(:,d0)-dgi(:,2).*phix(:,d0))+cij.*curlV(:,d0,3);
end
fb={[2,3,4;3,4,2],[1,4,3;4,3,1],[1,2,4;2,4,1],[1,3,2;3,2,1]};
for f=1:4
    d0=12+2*(f-1)+1; d1=d0+1; fb0=fb{f};
    a=fb0(1,1);b=fb0(1,2);c=fb0(1,3);ga=G{a};gb=G{b};gc=G{c};la=l(a);lb=l(b);
    phix(:,d0)=la.*lb.*gc(:,1);phiy(:,d0)=la.*lb.*gc(:,2);phiz(:,d0)=la.*lb.*gc(:,3);
    gx=lb.*ga(:,1)+la.*gb(:,1);gy=lb.*ga(:,2)+la.*gb(:,2);gz=lb.*ga(:,3)+la.*gb(:,3);
    curlV(:,d0,1)=gy.*gc(:,3)-gz.*gc(:,2);curlV(:,d0,2)=gz.*gc(:,1)-gx.*gc(:,3);curlV(:,d0,3)=gx.*gc(:,2)-gy.*gc(:,1);
    a2=fb0(2,1);b2=fb0(2,2);c2=fb0(2,3);ga2=G{a2};gb2=G{b2};gc2=G{c2};la2=l(a2);lb2=l(b2);
    phix(:,d1)=la2.*lb2.*gc2(:,1);phiy(:,d1)=la2.*lb2.*gc2(:,2);phiz(:,d1)=la2.*lb2.*gc2(:,3);
    gx2=lb2.*ga2(:,1)+la2.*gb2(:,1);gy2=lb2.*ga2(:,2)+la2.*gb2(:,2);gz2=lb2.*ga2(:,3)+la2.*gb2(:,3);
    curlV(:,d1,1)=gy2.*gc2(:,3)-gz2.*gc2(:,2);curlV(:,d1,2)=gz2.*gc2(:,1)-gx2.*gc2(:,3);curlV(:,d1,3)=gx2.*gc2(:,2)-gy2.*gc2(:,1);
end
end

% ===== Stiffness ===========================================================
function A = ned2curl3D(node, elem, gIdx, gSign, Ntot)
NT=size(elem,1); nLocal=20; [G,vol]=tetGeom(node,elem);
[lq,wq]=quadtet(4); nQ=length(wq);
nEntries=NT*nLocal*nLocal*nQ;
ii=zeros(nEntries,1);jj=zeros(nEntries,1);ss=zeros(nEntries,1);idx=0;
for q=1:nQ
    l=lq(q,:); [~,~,~,curlV]=ned2basis3D(l,G); w=6*wq(q)*vol;
    for p=1:nLocal
        gp=gIdx(:,p); sp=gSign(:,p);
        for qq=1:nLocal
            gq=gIdx(:,qq); sq=gSign(:,qq);
            dp=curlV(:,p,1).*curlV(:,qq,1)+curlV(:,p,2).*curlV(:,qq,2)+curlV(:,p,3).*curlV(:,qq,3);
            s=sp.*sq.*w.*dp;
            nxt=idx+1;idx=idx+NT;
            ii(nxt:idx)=gp;jj(nxt:idx)=gq;ss(nxt:idx)=s;
        end
    end
end
A = sparse(ii(1:idx),jj(1:idx),ss(1:idx),Ntot,Ntot);
end

% ===== Mass ================================================================
function M = ned2mass3D(node, elem, gIdx, gSign, Ntot)
NT=size(elem,1); nLocal=20; [G,vol]=tetGeom(node,elem);
[lq,wq]=quadtet(4); nQ=length(wq);
nEntries=NT*nLocal*nLocal*nQ;
ii=zeros(nEntries,1);jj=zeros(nEntries,1);ss=zeros(nEntries,1);idx=0;
for q=1:nQ
    l=lq(q,:); [phix,phiy,phiz]=ned2basis3D(l,G); w=6*wq(q)*vol;
    for p=1:nLocal
        gp=gIdx(:,p);sp=gSign(:,p);
        for qq=1:nLocal
            gq=gIdx(:,qq);sq=gSign(:,qq);
            dp=phix(:,p).*phix(:,qq)+phiy(:,p).*phiy(:,qq)+phiz(:,p).*phiz(:,qq);
            s=sp.*sq.*w.*dp;
            nxt=idx+1;idx=idx+NT;
            ii(nxt:idx)=gp;jj(nxt:idx)=gq;ss(nxt:idx)=s;
        end
    end
end
M = sparse(ii(1:idx),jj(1:idx),ss(1:idx),Ntot,Ntot);
end

% ===== RHS =================================================================
function b = ned2rhs3D(node, elem, f_rhs, gIdx, gSign, Ntot)
NT=size(elem,1); nLocal=20; [G,vol]=tetGeom(node,elem);
[lq,wq]=quadtet(4); b=zeros(Ntot,1);
v1=node(elem(:,1),:);v2=node(elem(:,2),:);v3=node(elem(:,3),:);v4=node(elem(:,4),:);
for q=1:length(wq)
    l=lq(q,:); px=l(1)*v1(:,1)+l(2)*v2(:,1)+l(3)*v3(:,1)+l(4)*v4(:,1);
    py=l(1)*v1(:,2)+l(2)*v2(:,2)+l(3)*v3(:,2)+l(4)*v4(:,2);
    pz=l(1)*v1(:,3)+l(2)*v2(:,3)+l(3)*v3(:,3)+l(4)*v4(:,3);
    [fx,fy,fz]=f_rhs(px,py,pz);
    [phix,phiy,phiz]=ned2basis3D(l,G);
    w=6*wq(q)*vol;
    for kk=1:nLocal
        c=w.*(fx.*phix(:,kk)+fy.*phiy(:,kk)+fz.*phiz(:,kk));
        b=b+accumarray(gIdx(:,kk),gSign(:,kk).*c,[Ntot,1]);
    end
end
end
% Note: RHS face transform not applied (face DOFs use identity mapping for RHS)

% ===== Error ===============================================================
function [eL2,eHc] = ned2err3D(node,elem,uh,u_ex,curl_ex,gIdx,gSign)
NT=size(elem,1);nLocal=20;[G,vol]=tetGeom(node,elem);
[lq,wq]=quadtet(6);eL2=0;eHc=0;
v1=node(elem(:,1),:);v2=node(elem(:,2),:);v3=node(elem(:,3),:);v4=node(elem(:,4),:);
for q=1:length(wq)
    l=lq(q,:);px=l(1)*v1(:,1)+l(2)*v2(:,1)+l(3)*v3(:,1)+l(4)*v4(:,1);
    py=l(1)*v1(:,2)+l(2)*v2(:,2)+l(3)*v3(:,2)+l(4)*v4(:,2);
    pz=l(1)*v1(:,3)+l(2)*v2(:,3)+l(3)*v3(:,3)+l(4)*v4(:,3);
    [uex,uey,uez]=u_ex(px,py,pz);[cex,cey,cez]=curl_ex(px,py,pz);
    [phix,phiy,phiz,curlV]=ned2basis3D(l,G);
    uhx=zeros(NT,1);uhy=zeros(NT,1);uhz=zeros(NT,1);cuhx=zeros(NT,1);cuhy=zeros(NT,1);cuhz=zeros(NT,1);
    for kk=1:nLocal
        uv=uh(gIdx(:,kk));ss=gSign(:,kk);
        uhx=uhx+ss.*uv.*phix(:,kk);uhy=uhy+ss.*uv.*phiy(:,kk);uhz=uhz+ss.*uv.*phiz(:,kk);
        cuhx=cuhx+ss.*uv.*curlV(:,kk,1);cuhy=cuhy+ss.*uv.*curlV(:,kk,2);cuhz=cuhz+ss.*uv.*curlV(:,kk,3);
    end
    ex=uhx-uex;ey=uhy-uey;ez=uhz-uez;ecx=cuhx-cex;ecy=cuhy-cey;ecz=cuhz-cez;
    wa=6*wq(q)*vol;
    eL2=eL2+sum(wa.*(ex.^2+ey.^2+ez.^2));
    eHc=eHc+sum(wa.*(ex.^2+ey.^2+ez.^2+ecx.^2+ecy.^2+ecz.^2));
end
eL2=sqrt(eL2);eHc=sqrt(eHc);
end

