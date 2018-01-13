function F = errProjection(p,x2d,X3d)
% % linear projection system parameterized by P
    X3d = double(X3d);
    x2d = double(x2d);
    P = (reshape(p, 4, 3))';
    x2dP = (P*X3d')';
    x2dP = x2dP./repmat(x2dP(:,3), 1, 3);
    x2d = x2d(:, 1:2);
    x2dP = x2dP(:, 1:2);
    F = x2d(:)-x2dP(:);
end