addpath(genpath(pwd));
for hh = [0.25, 0.125, 0.0625]
  [nd,el,bd]=squaremesh([0,1,0,1],hh);
  [~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
  bc=[2,3,1];
  eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
  sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
  A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el); K=A+M;
  % RHS
  [lq,wq]=quadtriangle(2); b=zeros(NE,1);
  x1=nd(el(:,1),:);x2=nd(el(:,2),:);x3=nd(el(:,3),:);
  area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
  area=abs(area2)/2; invA2=1./area2;
  g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
  g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
  g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
  f_rhs=@(x,y)(2+y.*(1-y));
  for q=1:length(wq)
    l=lq(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p3x=l(1)*g2(:,1)-l(2)*g1(:,1);
    c1=wq(q)*area.*fx.*p1x; c2=wq(q)*area.*fx.*p2x; c3=wq(q)*area.*fx.*p3x;
    b=b+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b=b+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b=b+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
  end
  % Boundary
  bdFlag_to_e=[2,3,1]; bdE=[];
  for k=1:3, isBd=bd(:,k)==1; if any(isBd), bdE=[bdE;eidx(isBd,bdFlag_to_e(k))]; end; end
  bdE=unique(bdE); freeE=setdiff(1:NE,bdE)';
  % Solve
  uh=zeros(NE,1); uh(freeE)=K(freeE,freeE)\b(freeE);
  % Residual
  res=K*uh-b;
  fprintf('h=%.4f (NE=%d): ||K*uh-b||/||b|| = %.2e\n',hh,NE,norm(res)/norm(b));
end
exit(0);
