addpath(genpath(pwd));
node=[0,0;1,0;0,1]; elem=[1,2,3];
[~,eidx,es]=edgeMesh2D(elem); NE=max(eidx(:)); NT=1; nLocal=8;
[lq,wq]=quadtriangle(4);
x1=node(elem(:,1),:);x2=node(elem(:,2),:);x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];
fprintf('g3=(%f,%f)\n',g3(1),g3(2));

% Manually trace M(7,7)
M77=0;
for q=1:length(wq)
    l=lq(q,:);
    phix7=l(1)*l(2)*g3(1); phiy7=l(1)*l(2)*g3(2);
    dotprod=phix7^2+phiy7^2;
    w=2*wq(q)*area;
    contrib=w*dotprod;
    M77=M77+contrib;
    fprintf('q%d: l=(%.4f,%.4f,%.4f) phix7=%.6f phiy7=%.6f dot=%.6f w=%.6f contrib=%.6f cum=%.6f\n',...
        q,l(1),l(2),l(3),phix7,phiy7,dotprod,w,contrib,M77);
end
fprintf('Manual M(7,7)=%.10f\n',M77);

% Now trace what the actual assembly code does
% Simulate the assembly loop
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
[~,~,es]=edgeMesh2D(elem);
localEdgeCols=[2,3,1]; localEdgeBases=[3,5,1];
gIdx_actual=zeros(NT,nLocal); gSign_actual=zeros(NT,nLocal);
for kk=1:3
    col=localEdgeCols(kk); eidk=eidx(:,col); sigk=es(:,col); b0=localEdgeBases(kk);
    gIdx_actual(:,b0)=2*(eidk-1)+1; gIdx_actual(:,b0+1)=2*(eidk-1)+2;
    gSign_actual(:,b0)=sigk; gSign_actual(:,b0+1)=sigk;
end
gIdx_actual(1,7)=2*NE+2*(1-1)+1; gIdx_actual(1,8)=2*NE+2*(1-1)+2;
fprintf('\nActual gIdx(1,:)='); fprintf('%d ',gIdx_actual(1,:)); fprintf('\n');
fprintf('Actual gSign(1,:)='); fprintf('%d ',gSign_actual(1,:)); fprintf('\n');

% Simulate the actual assembly
Msim=zeros(Ntot,Ntot); Ntot=2*NE+2*NT;
for q=1:length(wq)
    l=lq(q,:);
    phix=zeros(NT,nLocal); phiy=zeros(NT,nLocal);
    phix(:,1)=l(1)*g2(:,1)-l(2)*g1(:,1); phiy(:,1)=l(1)*g2(:,2)-l(2)*g1(:,2);
    c12=l(1)-l(2); phix(:,2)=c12.*phix(:,1); phiy(:,2)=c12.*phiy(:,1);
    phix(:,3)=l(2)*g3(:,1)-l(3)*g2(:,1); phiy(:,3)=l(2)*g3(:,2)-l(3)*g2(:,2);
    c23=l(2)-l(3); phix(:,4)=c23.*phix(:,3); phiy(:,4)=c23.*phiy(:,3);
    phix(:,5)=l(3)*g1(:,1)-l(1)*g3(:,1); phiy(:,5)=l(3)*g1(:,2)-l(1)*g3(:,2);
    c31=l(3)-l(1); phix(:,6)=c31.*phix(:,5); phiy(:,6)=c31.*phiy(:,5);
    phix(:,7)=l(1).*l(2).*g3(:,1); phiy(:,7)=l(1).*l(2).*g3(:,2);
    phix(:,8)=l(2).*l(3).*g1(:,1); phiy(:,8)=l(2).*l(3).*g1(:,2);
    w=2*wq(q)*area;
    for p=1:nLocal
        gp=gIdx_actual(1,p); sp=gSign_actual(1,p);
        for qq=1:nLocal
            gq=gIdx_actual(1,qq); sq=gSign_actual(1,qq);
            s=sp*sq*w*(phix(1,p)*phix(1,qq)+phiy(1,p)*phiy(1,qq));
            Msim(gp,gq)=Msim(gp,gq)+s;
        end
    end
end
fprintf('\nMsim(7,7)=%.10f\n',Msim(7,7));
% Check where DOF 7 gets contributions from
fprintf('\nNonzero M(7,:) from simulation:\n');
for j=1:8, if abs(Msim(7,j))>1e-15, fprintf('M(7,%d)=%.10f\n',j,Msim(7,j)); end; end
exit(0);
