function P = para2Proj(p)
% % This function computes the projection matrix given parameters 'p'
% % input p: 9x1 or 1x9
p = p(:);
p = double(p);
f = p(1);
cx = p(2); cy = p(3);
ea = p(4:6);
t = p(7:9);
K = [f 0 cx;
     0 f cy;
     0 0 1];
R = rotm('xyz', ea);
P = K*[R,t];