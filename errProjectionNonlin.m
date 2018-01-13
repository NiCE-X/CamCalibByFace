function F = errProjectionNonlin(p,x2d,X3d,w_skew,w_focal)
% % non-linear projection system. Also constrain on skew and focal length.
    X3d = double(X3d);
    x2d = double(x2d);
    p = double(p);
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
    F = [F; w_skew*skew; w_focal*(fx-fy)];
end