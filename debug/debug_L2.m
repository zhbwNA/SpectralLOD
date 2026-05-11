addpath(genpath(pwd));
[nd, el, bd] = squaremesh([0,1,0,1], 0.125);
[~, edgeIdx, edgeSign] = edgeMesh2D(el);
NE = max(edgeIdx(:));
bc=[2,3,1]; NT=size(el,1);
eid=zeros(NT,3); sig=zeros(NT,3);
eid(:,1)=edgeIdx(:,bc(1)); eid(:,2)=edgeIdx(:,bc(2)); eid(:,3)=edgeIdx(:,bc(3));
sig(:,1)=edgeSign(:,bc(1)); sig(:,2)=edgeSign(:,bc(2)); sig(:,3)=edgeSign(:,bc(3));

% Compute NE_1 interpolant DOFs
u_interp = zeros(NE,1);
cnt = zeros(NE,1);
for t=1:NT
    v=nd(el(t,:),:);
    for k=1:3
        switch k, case 1, va=2;vb=3; case 2, va=3;vb=1; case 3, va=1;vb=2; end
        a=v(va,:); b_pt=v(vb,:); L=norm(b_pt-a); tvec=(b_pt-a)/L;
        % Line integral of u·t_local
        [xi,w1d]=deal([-sqrt(3/5);0;sqrt(3/5)],[5/9;8/9;5/9]);
        val=0;
        for qi=1:3
            s=(xi(qi)+1)/2; pt=a+s*(b_pt-a);
            ux=pt(2)*(1-pt(2));  % y(1-y)
            val=val+w1d(qi)/2*L*ux*tvec(1);
        end
        % Global DOF contribution: sig * val (val uses LOCAL tangent)
        u_interp(eid(t,k)) = u_interp(eid(t,k)) + sig(t,k)*val;
        cnt(eid(t,k)) = cnt(eid(t,k)) + 1;
    end
end
u_interp = u_interp ./ max(cnt,1);

% Compute FE solution
A = assembleCurlCurl2D(nd, el);
M = assembleNedMass2D(nd, el);
K = A + M;

f_rhs=@(x,y)(2+y.*(1-y));
b=zeros(NE,1);
[lambda_q,wq]=quadtriangle(2);
x1=nd(el(:,1),:); x2=nd(el(:,2),:); x3=nd(el(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
for q=1:length(wq)
    l=lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx=f_rhs(px,py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p3x=l(1)*g2(:,1)-l(2)*g1(:,1);
    c1=wq(q)*area.*fx.*p1x; c2=wq(q)*area.*fx.*p2x; c3=wq(q)*area.*fx.*p3x;
    b=b+accumarray(eid(:,1),sig(:,1).*c1,[NE,1]);
    b=b+accumarray(eid(:,2),sig(:,2).*c2,[NE,1]);
    b=b+accumarray(eid(:,3),sig(:,3).*c3,[NE,1]);
end

bdFlag_to_edgeIdx=[2,3,1]; bdEdges=[];
for k=1:3, isBd=bd(:,k)==1; if any(isBd), bdEdges=[bdEdges;edgeIdx(isBd,bdFlag_to_edgeIdx(k))]; end; end
bdEdges=unique(bdEdges);
freeEdges=setdiff(1:NE,bdEdges)';
u_f=K(freeEdges,freeEdges)\b(freeEdges);
uh=zeros(NE,1); uh(freeEdges)=u_f;

% Compare interpolant vs FE solution
fprintf('||u_interp||=%.4e ||uh||=%.4e ||u_interp-uh||=%.4e\n',...
    norm(u_interp), norm(uh), norm(u_interp-uh));

% Check residual: does (A+M)*uh = b?
res = K*uh - b;
fprintf('||(A+M)uh - b|| = %.4e, ||b|| = %.4e\n', norm(res), norm(b));

% Check: does (A+M)*u_interp = b?
res_interp = K*u_interp - b;
fprintf('||(A+M)u_interp - b|| = %.4e\n', norm(res_interp));

% Compare L2 norms
[lambda_q4,wq4]=quadtriangle(4);
c1=2*(g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1));
c2=2*(g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1));
c3=2*(g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1));
uv1=uh(eid(:,1)); uv2=uh(eid(:,2)); uv3=uh(eid(:,3));
iv1=u_interp(eid(:,1)); iv2=u_interp(eid(:,2)); iv3=u_interp(eid(:,3));

err_L2_sq=0; u_L2_sq=0; uh_L2_sq=0; ui_L2_sq=0;
for q=1:length(wq4)
    l=lambda_q4(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    uex=py.*(1-py);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    
    uhx=sig(:,1).*uv1.*p1x+sig(:,2).*uv2.*p2x+sig(:,3).*uv3.*p3x;
    uhy=sig(:,1).*uv1.*p1y+sig(:,2).*uv2.*p2y+sig(:,3).*uv3.*p3y;
    uix=sig(:,1).*iv1.*p1x+sig(:,2).*iv2.*p2x+sig(:,3).*iv3.*p3x;
    uiy=sig(:,1).*iv1.*p1y+sig(:,2).*iv2.*p2y+sig(:,3).*iv3.*p3y;
    
    w_area=wq4(q)*area;
    err_L2_sq=err_L2_sq+sum(w_area.*((uhx-uex).^2+uhy.^2));
    u_L2_sq=u_L2_sq+sum(w_area.*(uex.^2));
    uh_L2_sq=uh_L2_sq+sum(w_area.*(uhx.^2+uhy.^2));
    ui_L2_sq=ui_L2_sq+sum(w_area.*((uix-uex).^2+uiy.^2));
end
fprintf('\n||u||_L2=%.4f ||uh||_L2=%.4f ||u-uh||_L2=%.4f ||u-ui||_L2=%.4f\n',...
    sqrt(u_L2_sq), sqrt(uh_L2_sq), sqrt(err_L2_sq), sqrt(ui_L2_sq));
fprintf('Relative L2 error: FE=%.2f%%, Interpolant=%.2f%%\n',...
    100*sqrt(err_L2_sq/u_L2_sq), 100*sqrt(ui_L2_sq/u_L2_sq));

% Key question: is the FE solution close to the interpolant?
diff_L2_sq=0;
for q=1:length(wq4)
    l=lambda_q4(q,:);
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);
    dux=sig(:,1).*(uv1-iv1).*p1x+sig(:,2).*(uv2-iv2).*p2x+sig(:,3).*(uv3-iv3).*p3x;
    duy=sig(:,1).*(uv1-iv1).*p1y+sig(:,2).*(uv2-iv2).*p2y+sig(:,3).*(uv3-iv3).*p3y;
    diff_L2_sq=diff_L2_sq+sum(wq4(q)*area.*(dux.^2+duy.^2));
end
fprintf('||uh - u_interp||_L2 = %.4e\n', sqrt(diff_L2_sq));
exit(0);
