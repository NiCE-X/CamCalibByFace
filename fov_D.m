function fov = fov_D(d,x,X)
% % This fuction computes the fov given the wanted z-distance:d and the
% % pixel distance between two eyes: x. The 3D distance between two eyes
% % is known as: X
w = 512;
f = d*x/X;
fov = 2*atand(w/2/f);