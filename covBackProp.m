function p_cov = covBackProp(p_v, X3d_v, x2d_cov)
% % This function computes the theoretical covariance matrix of estimated
% % parameters 'p' at the MLE 'p_v'.
% % input p_v: 9x1 or 1x9, where to evaluate
% % input X3d_v: nx3, 3D landmarks, homogeneous coord.
% % input x2d_cov: 2nx2n, the uncertainty covariance of 2d landmark measures.
p = sym('p', [9,1]);
K = [p(1) 0 p(2);
     0 p(1) p(3);
     0 0 1];
Rx = [1 0 0;
    0 cos(p(4)*pi/180) -sin(p(4)*pi/180);
    0 sin(p(4)*pi/180) cos(p(4)*pi/180)];
Ry = [cos(p(5)*pi/180) 0 sin(p(5)*pi/180);
    0 1 0;
    -sin(p(5)*pi/180) 0 cos(p(5)*pi/180)];
Rz = [cos(p(6)*pi/180) -sin(p(6)*pi/180) 0;
    sin(p(6)*pi/180) cos(p(6)*pi/180) 0;
    0 0 1];
R = Rx*Ry*Rz;
t = p(7:9);
P = K*[R,t];

X3d_v = [X3d_v, ones(size(X3d_v,1), 1)];
x2dP = (P*X3d_v')';
x2dP = x2dP./repmat(x2dP(:,3), 1, 3);
x2dP = x2dP(:, 1:2);
x2dP = x2dP';
x2dP = x2dP(:);
J = jacobian(x2dP, p);
for i = 1:9
    eval(['p', num2str(i), '=p_v(i);']);
end
J_v = eval(J);
p_cov = inv(J_v'/x2d_cov*J_v);