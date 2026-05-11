addpath(genpath(pwd));
for hh = [0.25, 0.125, 0.0625]
  [nd,el,bd]=squaremesh([0,1,0,1],hh);
  [~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
  bc=[2,3,1];
  eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
  sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
  A=assembleCurlCurl2D(nd,el); M=assembleNedMass2D(nd,el); K=A+M;
  % Compute exact DOFs
  u_ex=zeros(NE,1); cnt=zeros(NE,1);
  for t=1:NT
    v=nd(el(t,:),:);
    for k=1:3
      switch k, case 1, va=2;vb=3; case 2, va=3;vb=1; case 3, va=1;vb=2; end
      a2=v(va,:); b2=v(vb,:); L=norm(b2-a2); tvec=(b2-a2)/L;
      xi=[-sqrt(3/5);0;sqrt(3/5)]; ww=[5/9;8/9;5/9]; val=0;
      for qi=1:3, s=(xi(qi)+1)/2; pt=a2+s*(b2-a2); val=val+ww(qi)/2*L*pt(2)*(1-pt(2))*tvec(1); end
      u_ex(eid(t,k))=u_ex(eid(t,k))+sig(t,k)*val; cnt(eid(t,k))=cnt(eid(t,k))+1;
    end
  end
  u_ex=u_ex./max(cnt,1);
  % Compute RHS via mass: b_M = (A+M)*u_ex
  b_M = K*u_ex;
  % Compute RHS via quadrature
  b_Q = zeros(NE,1);
  [lq,wq]=quadtriangle(2);
  x1=nd(el(:,1),:); x2=nd(el(:,2),:); x3=nd(el(:,3),:);
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
    b_Q=b_Q+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b_Q=b_Q+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b_Q=b_Q+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
  end
  r=norm(b_M-b_Q)/norm(b_Q);
  fprintf('h=%.4f: ||K*u_ex-b_Q||/||b_Q|| = %.2e\n',hh,r);
end
exit(0);
