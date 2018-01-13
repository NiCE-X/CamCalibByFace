function F = errProjectionNonlin3(p,x2d,X3d,x2d_w,X3d_w,w)
% % non-linear projection system. skew and focal length are clamped.
% % This is the weighted version
    X3d = double(X3d);
    x2d = double(x2d);
    P = para2Proj(p);  
    x2dP = (P*X3d')';
    x2dP = x2dP./repmat(x2dP(:,3), 1, 3);
    x2d2 = x2d(:, 1:2);
    x2dP = x2dP(:, 1:2);
    
    X3d_w = double(X3d_w);
    x2d_w = double(x2d_w);
    x2dP_w = (P*X3d_w')';
    x2dP_w = x2dP_w./repmat(x2dP_w(:,3), 1, 3);
    x2d2_w = x2d_w(:, 1:2);
    x2dP_w = x2dP_w(:, 1:2);
    F = [x2d2(:)-x2dP(:); w*(x2d2_w(:)-x2dP_w(:))];
end