addpath(genpath(pwd));
hh=0.125;
[nd,el,bd]=squaremesh([0,1,0,1],hh);
[~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
bc=[2,3,1];
eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el); K=A+M;
b_v = assembleNedRHS2D(nd,el,@(x,y)(2+y.*(1-y)));
bdE_v = findBoundaryEdges2D(el,bd);
freeE_v = setdiff(1:NE,bdE_v)';
uh=zeros(NE,1); uh(freeE_v)=K(freeE_v,freeE_v)\b_v(freeE_v);

% Error via verify function
[eL2_v, eHc_v] = computeNedError2D(nd,el,uh,@(x,y)deal(y.*(1-y),zeros(size(x))),@(x,y)(2*y-1));

% Error via debug function
x1=nd(el(:,1),:);x2=nd(el(:,2),:);x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
c1=2*(g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1));
c2=2*(g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1));
c3=2*(g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1));
uv1=uh(eid(:,1)); uv2=uh(eid(:,2)); uv3=uh(eid(:,3));
[lq4,wq4]=quadtriangle(4);
eL2_d=0; eHc_d=0;
for q=1:length(wq4)
    l=lq4(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    uex=py.*(1-py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    uhx=sig(:,1).*uv1.*p1x+sig(:,2).*uv2.*p2x+sig(:,3).*uv3.*p3x;
    uhy=sig(:,1).*uv1.*p1y+sig(:,2).*uv2.*p2y+sig(:,3).*uv3.*p3y;
    cur=2*py-1; cuh=sig(:,1).*uv1.*c1+sig(:,2).*uv2.*c2+sig(:,3).*uv3.*c3;
    wa=2*wq4(q)*area;
    eL2_d=eL2_d+sum(wa.*((uhx-uex).^2+uhy.^2));
    eHc_d=eHc_d+sum(wa.*((uhx-uex).^2+uhy.^2+(cuh-cur).^2));
end
fprintf('Verify: L2=%.4e Hc=%.4e\n',eL2_v,eHc_v);
fprintf('Debug:  L2=%.4e Hc=%.4e\n',sqrt(eL2_d),sqrt(eHc_d));
exit(0);

function bdEdges=findBoundaryEdges2D(elem,bdFlag)
[~,eidx]=edgeMesh2D(elem); m=[2,3,1]; bdEdges=[];
for k=1:3, isBd=bdFlag(:,k)==1; if any(isBd), bdEdges=[bdEdges;eidx(isBd,m(k))]; end; end
bdEdges=unique(bdEdges);
end

function b=assembleNedRHS2D(node,elem,f_rhs)
[~,ei,es]=edgeMesh2D(elem); NE=max(ei(:)); NT=size(elem,1);
[lq,wq]=quadtriangle(2);
x1=node(elem(:,1),:);x2=node(elem(:,2),:);x3=node(elem(:,3),:);
a2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
a=abs(a2)/2; iA=1./a2;
g1=[(x2(:,2)-x3(:,2)).*iA,(x3(:,1)-x2(:,1)).*iA];
g2=[(x3(:,2)-x1(:,2)).*iA,(x1(:,1)-x3(:,1)).*iA];
g3=[(x1(:,2)-x2(:,2)).*iA,(x2(:,1)-x1(:,1)).*iA];
bc=[2,3,1]; eid=[ei(:,bc(1)),ei(:,bc(2)),ei(:,bc(3))]; sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
b=zeros(NE,1);
for q=1:length(wq)
    l=lq(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py); fy=zeros(size(fx));
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    c1=2*wq(q)*a.*(fx.*p1x+fy.*p1y); c2=2*wq(q)*a.*(fx.*p2x+fy.*p2y); c3=2*wq(q)*a.*(fx.*p3x+fy.*p3y);
    b=b+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b=b+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b=b+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
end
end

function [eL2,eHc]=computeNedError2D(node,elem,uh,u_exact,curl_exact)
[~,ei,es]=edgeMesh2D(elem); NT=size(elem,1);
[lq,wq]=quadtriangle(4);
x1=node(elem(:,1),:);x2=node(elem(:,2),:);x3=node(elem(:,3),:);
a2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
a=abs(a2)/2; iA=1./a2;
g1=[(x2(:,2)-x3(:,2)).*iA,(x3(:,1)-x2(:,1)).*iA];
g2=[(x3(:,2)-x1(:,2)).*iA,(x1(:,1)-x3(:,1)).*iA];
g3=[(x1(:,2)-x2(:,2)).*iA,(x2(:,1)-x1(:,1)).*iA];
c1=2*(g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1));
c2=2*(g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1));
c3=2*(g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1));
bc=[2,3,1]; eid=[ei(:,bc(1)),ei(:,bc(2)),ei(:,bc(3))]; sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
uv1=uh(eid(:,1)); uv2=uh(eid(:,2)); uv3=uh(eid(:,3));
eL2=0; eHc=0;
for q=1:length(wq)
    l=lq(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    [uex,uey]=u_exact(px,py); curlex=curl_exact(px,py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    uhx=sig(:,1).*uv1.*p1x+sig(:,2).*uv2.*p2x+sig(:,3).*uv3.*p3x;
    uhy=sig(:,1).*uv1.*p1y+sig(:,2).*uv2.*p2y+sig(:,3).*uv3.*p3y;
    cuh=sig(:,1).*uv1.*c1+sig(:,2).*uv2.*c2+sig(:,3).*uv3.*c3;
    ex=uhx-uex; ey=uhy-uey; ec=cuh-curlex;
    w_area=2*wq(q)*a;
    eL2=eL2+sum(w_area.*(ex.^2+ey.^2));
    eHc=eHc+sum(w_area.*(ex.^2+ey.^2+ec.^2));
end
eL2=sqrt(eL2); eHc=sqrt(eHc);
end
