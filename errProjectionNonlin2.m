function F = errProjectionNonlin2(p,x2d,X3d)
% % non-linear projection system. skew and focal length are clamped
    X3d = double(X3d);
    x2d = double(x2d);
    P = para2Proj(p);  
    x2dP = (P*X3d')';
    x2dP = x2dP./repmat(x2dP(:,3), 1, 3);
    x2d2 = x2d(:, 1:2);
    x2dP = x2dP(:, 1:2);
    F = x2d2(:)-x2dP(:);
end