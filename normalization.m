function [T2,T3,xn,Xn]=normalization(x2d,X3d)
x=x2d;
X=X3d;
mean2d = mean(x(:, 1:2), 1);    %归一化
mean3d = mean(X(:, 1:3), 1);
norm2d = zeros(size(x,1), 1);
norm3d = zeros(size(X,1), 1);
xmean = x(:, 1:2)-repmat(mean2d, size(x,1), 1);
Xmean = X(:, 1:3)-repmat(mean3d, size(X,1), 1);
for i = 1:size(x,1) %计算每个点到中心的距离
    norm2d(i,1) = norm(xmean(i,:));
    norm3d(i,1) = norm(Xmean(i,:));
end
% %使点到中心距离的均方根为1 from 高伟课件
% norm2 = sqrt(mean(norm2d.^2));%点到中心的距离看作方差
% norm3 = sqrt(mean(norm3d.^2));

% 使2d点到中心的平均距离为sqrt2, 3d点为sqrt3. from Hartley
norm2 = mean(norm2d)/sqrt(2);
norm3 = mean(norm3d)/sqrt(3);

T2 = inv([norm2, 0, mean2d(1);...
         0, norm2, mean2d(2);...
         0, 0, 1]);
T3 = inv([norm3, 0, 0, mean3d(1);...
         0, norm3, 0, mean3d(2);...
         0, 0, norm3, mean3d(3);...
         0, 0, 0, 1]);
% T2 = eye(3); T3 = eye(4);
xn = transpose(T2*x');%减均值除方差
Xn = transpose(T3*X');
