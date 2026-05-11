addpath(genpath(pwd));
node=[0,0;1,0;0,1]; elem=[1,2,3];
M=assembleNedMass2D(node,elem);
A=assembleCurlCurl2D(node,elem);
fprintf('Fixed NE_1 single-element matrices:\n');
fprintf('Stiffness:\n'); disp(full(A));
fprintf('Mass:\n'); disp(full(M));
% Check manually
x1=node(1,:);x2=node(2,:);x3=node(3,:);
area2=(x2(1)-x1(1))*(x3(2)-x1(2))-(x3(1)-x1(1))*(x2(2)-x1(2));
area=abs(area2)/2;
g1=[(x2(2)-x3(2))/area2,(x3(1)-x2(1))/area2];
g2=[(x3(2)-x1(2))/area2,(x1(1)-x3(1))/area2];
g3=[(x1(2)-x2(2))/area2,(x2(1)-x1(1))/area2];
[lambda_q,wq]=quadtriangle(2);
M_man=zeros(3);
for q=1:length(wq)
    l=lambda_q(q,:);
    p1=l(2)*g3-l(3)*g2; p2=l(3)*g1-l(1)*g3; p3=l(1)*g2-l(2)*g1;
    M_man(1,1)=M_man(1,1)+wq(q)*area*dot(p1,p1);
    M_man(1,2)=M_man(1,2)+wq(q)*area*dot(p1,p2);
    M_man(1,3)=M_man(1,3)+wq(q)*area*dot(p1,p3);
    M_man(2,2)=M_man(2,2)+wq(q)*area*dot(p2,p2);
    M_man(2,3)=M_man(2,3)+wq(q)*area*dot(p2,p3);
    M_man(3,3)=M_man(3,3)+wq(q)*area*dot(p3,p3);
end
M_man(2,1)=M_man(1,2);M_man(3,1)=M_man(1,3);M_man(3,2)=M_man(2,3);
fprintf('Manual mass:\n'); disp(M_man);
fprintf('Diff:\n'); disp(full(M)-M_man);
% Check which global edge is which
[~,edgeIdx,edgeSign]=edgeMesh2D(elem);
fprintf('edgeIdx: '); fprintf('%d ', edgeIdx); fprintf('\n');
fprintf('edgeSign: '); fprintf('%d ', edgeSign); fprintf('\n');
exit(0);
