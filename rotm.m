function R = rotm(str, theta)
% % input str: rotation axis order, e.g. 'xyz'
% % input theta: 1x3 or 3x1 roation angle, in degree, not radian!
R = eul2rotm(str(1), theta(1))*eul2rotm(str(2), theta(2))*eul2rotm(str(3), theta(3));