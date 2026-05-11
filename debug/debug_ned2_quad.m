addpath(genpath(pwd));
[l,w]=quadtriangle(4);
sum_l12l22 = sum(w .* l(:,1).^2 .* l(:,2).^2);
fprintf('Sum w*l1^2*l2^2 = %.10f (exact = 1/180 = %.10f)\n', sum_l12l22, 1/180);
% Check at each point
for q=1:size(l,1)
    fprintf('q%d: (%f,%f,%f) w=%f  l1^2*l2^2=%e\n',q,l(q,1),l(q,2),l(q,3),w(q),l(q,1)^2*l(q,2)^2);
end
exit(0);
