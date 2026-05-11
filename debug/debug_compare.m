addpath(genpath(pwd));
hh=0.125;
[nd,el,bd]=squaremesh([0,1,0,1],hh);
[~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
bc=[2,3,1];
eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el); K=A+M;

% RHS from verify_ned1_2D function
f_rhs=@(x,y)(2+y.*(1-y));
b_verify=assembleNedRHS2D(nd,el,f_rhs);

% RHS from debug_conv_fixed
[lq,wq]=quadtriangle(2); b_debug=zeros(NE,1);
x1=nd(el(:,1),:);x2=nd(el(:,2),:);x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
for q=1:length(wq)
    l=lq(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p3x=l(1)*g2(:,1)-l(2)*g1(:,1);
    c1=2*wq(q)*area.*fx.*p1x; c2=2*wq(q)*area.*fx.*p2x; c3=2*wq(q)*area.*fx.*p3x;
    b_debug=b_debug+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b_debug=b_debug+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b_debug=b_debug+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
end
fprintf('||b_verify||=%.4e ||b_debug||=%.4e ||diff||=%.4e\n',norm(b_verify),norm(b_debug),norm(b_verify-b_debug));
exit(0);

function b=assembleNedRHS2D(node,elem,f_rhs)
[~,edgeIdx,edgeSign]=edgeMesh2D(elem); NE=max(edgeIdx(:)); NT=size(elem,1);
[lambda_q,weight]=quadtriangle(2); nQuad=length(weight);
x1=node(elem(:,1),:);x2=node(elem(:,2),:);x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
bc=[2,3,1];
eid=zeros(NT,3); sig=zeros(NT,3);
eid(:,1)=edgeIdx(:,bc(1)); eid(:,2)=edgeIdx(:,bc(2)); eid(:,3)=edgeIdx(:,bc(3));
sig(:,1)=edgeSign(:,bc(1)); sig(:,2)=edgeSign(:,bc(2)); sig(:,3)=edgeSign(:,bc(3));
b=zeros(NE,1);
for q=1:nQuad
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
