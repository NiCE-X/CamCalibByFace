function R = eul2rotm(axis, theta)
% % eulerian angle to rotation matrix
% % input axis: 'x', 'y' or 'z' axis
% % input theta: rotation angle, in degree, not radian!!
if strcmp(axis, 'x')
    R = [1 0 0;
        0 cosd(theta) -sind(theta);
        0 sind(theta) cosd(theta)];
elseif strcmp(axis, 'y')
    R = [cosd(theta) 0 sind(theta);
        0 1 0;
        -sind(theta) 0 cosd(theta)];
elseif strcmp(axis, 'z')
    R = [cosd(theta) -sind(theta) 0;
        sind(theta) cosd(theta) 0;
        0 0 1];
else
    disp('axis must be one of x, y, z!');
    assert(1 == 0);
end