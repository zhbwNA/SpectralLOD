addpath(genpath(pwd));
hh=0.125;
[nd,el,bd]=squaremesh([0,1,0,1],hh);
[~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
bc=[2,3,1];
eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el);
% Solve using verify code path
b_v = assembleNedRHS2D(nd,el,@(x,y)(2+y.*(1-y)));
bdE_v = findBoundaryEdges2D(el,bd);
freeE_v = setdiff(1:NE,bdE_v)';
K=A+M;
uh_v = zeros(NE,1); uh_v(freeE_v)=K(freeE_v,freeE_v)\b_v(freeE_v);
% Solve using debug code path
[lq,wq]=quadtriangle(2); b_d=zeros(NE,1);
x1=nd(el(:,1),:);x2=nd(el(:,2),:);x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
for q=1:length(wq)
    l=lq(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=2+py.*(1-py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p3x=l(1)*g2(:,1)-l(2)*g1(:,1);
    c1=2*wq(q)*area.*fx.*p1x; c2=2*wq(q)*area.*fx.*p2x; c3=2*wq(q)*area.*fx.*p3x;
    b_d=b_d+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b_d=b_d+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b_d=b_d+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
end
bdFlag_to_e=[2,3,1]; bdE_d=[];
for k=1:3, isBd=bd(:,k)==1; if any(isBd), bdE_d=[bdE_d;eidx(isBd,bdFlag_to_e(k))]; end; end
bdE_d=unique(bdE_d); freeE_d=setdiff(1:NE,bdE_d)';
uh_d=zeros(NE,1); uh_d(freeE_d)=K(freeE_d,freeE_d)\b_d(freeE_d);
% Compare
fprintf('||b_v||=%.4e ||b_d||=%.4e ||b_v-b_d||=%.4e\n',norm(b_v),norm(b_d),norm(b_v-b_d));
fprintf('||uh_v||=%.4e ||uh_d||=%.4e ||uh_v-uh_d||=%.4e\n',norm(uh_v),norm(uh_d),norm(uh_v-uh_d));
% Check boundary detection
fprintf('bdE_v: %d edges, bdE_d: %d edges\n',length(bdE_v),length(bdE_d));
fprintf('freeE_v: %d, freeE_d: %d\n',length(freeE_v),length(freeE_d));
fprintf('setdiff(bdE_v,bdE_d): %d edges differ\n',length(setdiff(bdE_v,bdE_d)));
exit(0);

function b=assembleNedRHS2D(node,elem,f_rhs)
[~,edgeIdx,edgeSign]=edgeMesh2D(elem); NE=max(edgeIdx(:)); NT=size(elem,1);
[lambda_q,weight]=quadtriangle(2);
x1=node(elem(:,1),:);x2=node(elem(:,2),:);x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
bc=[2,3,1];
eid=[edgeIdx(:,bc(1)),edgeIdx(:,bc(2)),edgeIdx(:,bc(3))];
sig=[edgeSign(:,bc(1)),edgeSign(:,bc(2)),edgeSign(:,bc(3))];
b=zeros(NE,1);
for q=1:length(weight)
    l=lambda_q(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py); fy=zeros(size(fx));
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    c1=2*weight(q)*area.*(fx.*p1x+fy.*p1y);
    c2=2*weight(q)*area.*(fx.*p2x+fy.*p2y);
    c3=2*weight(q)*area.*(fx.*p3x+fy.*p3y);
    b=b+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b=b+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b=b+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
end
end

function bdEdges=findBoundaryEdges2D(elem,bdFlag)
[~,edgeIdx]=edgeMesh2D(elem);
bdFlag_to_edgeIdx=[2,3,1]; bdEdges=[];
for k=1:3, isBd=bdFlag(:,k)==1; if any(isBd), bdEdges=[bdEdges;edgeIdx(isBd,bdFlag_to_edgeIdx(k))]; end; end
bdEdges=unique(bdEdges);
end
