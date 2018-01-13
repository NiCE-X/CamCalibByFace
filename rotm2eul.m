function angles = rotm2eul(m)
% % This function converts a rotation matrix to euler angles in 'xyz' order
% % algorithm: explicitly write m(alpha, beta, gama) and observe. According
% % to online resources, the solution in not unique and there may be
% % singular value errors. Watch out. This just provides one of the
% % possible correct solutions. Has been tested for correctness.
temp = m*m'-eye(3);
assert(sum(abs(temp(:))) < 1e-6);
assert(abs(det(m)-1) < 1e-6);    % assert on rotation matrix ( reflective not OK).

gama = atan2d(-m(1,2), m(1,1));
alpha = atan2d(-m(2,3), m(3,3));
temp = sqrt(m(1,2)^2 + m(1,1)^2);
assert(abs(temp) > 1e-6);   % sigularity
beta = atan2d(m(1,3), temp);
angles = [alpha, beta, gama];