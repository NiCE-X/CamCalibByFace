function d = D_fov(fov,x,X)
% % This fuction computes the camera z-distance given the wanted fov and the
% % pixel distance between two eyes: x. The 3D distance between two eyes
% % is known as: X
w = 1024;
f = w/2/tand(fov/2);
d = f*X/x;