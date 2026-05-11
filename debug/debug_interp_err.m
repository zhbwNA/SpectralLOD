addpath(genpath(pwd));
for hh = [0.5, 0.25, 0.125, 0.0625]
  [nd,el,~]=squaremesh([0,1,0,1],hh);
  [~,eidx,es]=edgeMesh2D(el); NE=max(eidx(:)); NT=size(el,1);
  bc=[2,3,1];
  eid=[eidx(:,bc(1)),eidx(:,bc(2)),eidx(:,bc(3))];
  sig=[es(:,bc(1)),es(:,bc(2)),es(:,bc(3))];
  % Interpolant DOFs
  uI=zeros(NE,1); cnt=zeros(NE,1);
  for t=1:NT
    v=nd(el(t,:),:);
    for k=1:3
      switch k, case 1, va=2;vb=3; case 2, va=3;vb=1; case 3, va=1;vb=2; end
      a2=v(va,:); b2=v(vb,:); L=norm(b2-a2); tvec=(b2-a2)/L;
      xi=[-sqrt(3/5);0;sqrt(3/5)]; ww=[5/9;8/9;5/9]; val=0;
      for qi=1:3, s=(xi(qi)+1)/2; pt=a2+s*(b2-a2); val=val+ww(qi)/2*L*pt(2)*(1-pt(2))*tvec(1); end
      uI(eid(t,k))=uI(eid(t,k))+sig(t,k)*val; cnt(eid(t,k))=cnt(eid(t,k))+1;
    end
  end
  uI=uI./max(cnt,1);
  % L2 error of interpolant
  x1=nd(el(:,1),:);x2=nd(el(:,2),:);x3=nd(el(:,3),:);
  area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
  area=abs(area2)/2; invA2=1./area2;
  g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
  g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
  g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
  [lq4,wq4]=quadtriangle(4);
  uv1I=uI(eid(:,1)); uv2I=uI(eid(:,2)); uv3I=uI(eid(:,3));
  eI=0;
  for q=1:length(wq4)
    l=lq4(q,:); px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1); py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    uex=py.*(1-py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    uhx=sig(:,1).*uv1I.*p1x+sig(:,2).*uv2I.*p2x+sig(:,3).*uv3I.*p3x;
    uhy=sig(:,1).*uv1I.*p1y+sig(:,2).*uv2I.*p2y+sig(:,3).*uv3I.*p3y;
    wa=wq4(q)*area;
    eI=eI+sum(wa.*((uhx-uex).^2+uhy.^2));
  end
  fprintf('h=%.4f: ||u-uI||_L2=%.4e\n',hh,sqrt(eI));
end
exit(0);
