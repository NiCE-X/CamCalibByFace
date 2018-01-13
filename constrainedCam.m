function [P,K,R,t,rmse,p]=constrainedCam(X3d, x2d, K,R,t,rmse, OptDisplay)
% % This function puts constraints on camera skewness to 0, and fx=fy.
% % First gentle constraint then hard constraint.
if nargin < 7
    OptDisplay = 'Final';
end
ea = rotm2eul(R);
fx = K(1,1);
fy = K(2,2);
skew = K(1,2);
cx = K(1,3);
cy = K(2,3);
e_skew = abs(skew);
e_focal = abs(fx-fy);
lambda = 3;     % small lambda more gentle.
options = optimset('Algorithm','levenberg-marquardt', 'Display', OptDisplay);
p = [fx, fy, cx, cy, skew, ea, t'];
p = double(p);
count = 0;
flag = 0;
for i = 1:50     
    count = count+1;
    if ~strcmp(OptDisplay, 'off')
        display(count);
    end
    if e_skew/rmse > lambda || e_focal/rmse > lambda    % gentlely add constraints on skew and focal
        w_skew = double(rmse/e_skew*lambda);
        w_focal = double(rmse/e_focal*lambda);
        if abs(e_skew)<1e-4     % in case of one of these error approaches 0 and w turns INF
            w_skew = 0;
        end
        if abs(e_focal)<1e-4
            w_focal = 0;
        end
        p = lsqnonlin(@(p) errProjectionNonlin(p,x2d,X3d,w_skew, w_focal), p,[],[],options);
    else        % hard constraint at last
        f = (p(1)+p(2))/2;
        skew = 0;
        p([2,5]) = [];
        p(1) = f;
        p = lsqnonlin(@(p) errProjectionNonlin2(p,x2d,X3d), p,[],[],options);
        p = [p(1), p(1), p(2:3), skew, p(4:9)];
        flag = 1;   % finished.
    end
    
    fx = p(1); fy = p(2);
    cx = p(3); cy = p(4);
    skew = p(5);
    ea = p(6:8);
    t = p(9:11)';
    K = [fx skew cx;
         0 fy cy;
         0 0 1];
    R = rotm('xyz', ea);
    P = K*[R,t];
    
    x2dP = (P*X3d')';
    x2dP = x2dP./repmat(x2dP(:,3), 1, 3);
    x2d2 = x2d(:, 1:2);
    x2dP = x2dP(:, 1:2);
    F = x2d2(:)-x2dP(:);
    rmse=sqrt(sum((F).^2)/size(x2d,1));
    e_skew = abs(skew);
    e_focal = abs(fx-fy);
    if(flag == 1)
        break;
    end
end
p([2,5]) = [];
p = p(:);