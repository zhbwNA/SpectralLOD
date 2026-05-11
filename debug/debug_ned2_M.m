addpath(genpath(pwd));
node=[0,0;1,0;0,1]; elem=[1,2,3];
M=assembleNed2Mass2D(node,elem);
fprintf('M(7,7)=%.10f  M(8,8)=%.10f\n', full(M(7,7)), full(M(8,8)));
fprintf('M(1,1)=%.10f\n', full(M(1,1)));
% Manually compute M(7,7)
x1=0;y1=0;x2=1;y2=0;x3=0;y3=1;
area2=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
area=abs(area2)/2;
g3=[(y1-y2)/area2,(x2-x1)/area2];
[lq,wq]=quadtriangle(4);
M77=0;
for q=1:length(wq)
    l=lq(q,:);
    phix=l(1)*l(2)*g3(1); phiy=l(1)*l(2)*g3(2);
    M77=M77+2*wq(q)*area*(phix^2+phiy^2);
end
fprintf('Manual M(7,7)=%.10f\n', M77);
% Check DOF indexing
[~,eidx,es]=edgeMesh2D(elem);
fprintf('NE=%d\n', max(eidx(:)));
% What are gIdx(1,7) and gIdx(1,8)?
NE=max(eidx(:)); NT=1;
localEdgeCols=[2,3,1]; localEdgeBases=[3,5,1];
gIdx=zeros(NT,8);
for kk=1:3
    col=localEdgeCols(kk); eidk=eidx(:,col); sigk=es(:,col); b0=localEdgeBases(kk);
    gIdx(:,b0)=2*(eidk-1)+1; gIdx(:,b0+1)=2*(eidk-1)+2;
end
gIdx(1,7)=2*NE+2*(1-1)+1; gIdx(1,8)=2*NE+2*(1-1)+2;
fprintf('gIdx(1,7:8) = %d %d\n', gIdx(1,7), gIdx(1,8));
exit(0);
