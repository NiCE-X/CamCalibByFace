function [P,K,R,t,rmse]=resectioning(X3d0, x2d0, OptDisplay)
% % This function uses Direct Linear Transform (DLT) to compute projectiong matrix. 
% % Input points are in homogeneous coord form.
% % input X3d: nx4
% % input x2d: nx3
% % input display: whether to display optimization information.
    if nargin < 3
        OptDisplay = 'Final';
    end
    [T2,T3,x2d,X3d]=normalization(x2d0,X3d0);   %归一化  
% Hartley's method   
% Direct Linear Transform (DLT)
    A=zeros(2*size(x2d,1),12);
    for i=1:size(x2d, 1)
        A(2*i-1,:)=[0 0 0 0, -X3d(i,:), x2d(i,2)*X3d(i,:)];
        A(2*i,:)=[X3d(i,:), 0 0 0 0, -x2d(i,1)*X3d(i,:)];
    end
    [U, D, V] = svd(A, 0);
    p0 = V(:, end);
    options = optimset('Algorithm','levenberg-marquardt', 'Display', OptDisplay);
    p = lsqnonlin(@(p) errProjection(p,x2d,X3d), p0,[],[],options);
    P = (reshape(p, 4, 3))';
    P = T2\P*T3;        % 去归一化
    P = P/P(3,4);       % This is very important. Although P and -P will produce the same 2D projection points.
%     The following QR decomposition will result in multually negative
%     Rotation matrix, which are different in handness. I think this
%     normalizing operation forces the camera Coord's z axis outwards.
    
    M=P(:,1:3);
    [R,K]=qr(inv(M));
    f = diag(K);
    ind=find(f<0);
    F=[1 0 0;0 1 0;0 0 1];
    for i=1:length(ind)
        F(ind(i),:)=-F(ind(i),:);
    end
    R=R*F';
    K=F*K;
    R=R';
    K=inv(K);
    K=K/K(3,3);
    t=R/M*P(:,4);
    
    if abs(det(R)+1)<1e-6     % sometimes it still happens
        P = -P;
        M=P(:,1:3);
        [R,K]=qr(inv(M));
        f = diag(K);
        ind=find(f<0);
        F=[1 0 0;0 1 0;0 0 1];
        for i=1:length(ind)
            F(ind(i),:)=-F(ind(i),:);
        end
        R=R*F';
        K=F*K;
        R=R';
        K=inv(K);
        K=K/K(3,3);
        t=R/M*P(:,4);
    end
        
    xx=(P*X3d0')';
    xx=xx./repmat(xx(:,3),1,3);
    rmse=sqrt(sum(sum((x2d0(:,1:2)-xx(:,1:2)).^2))/size(x2d0,1));
end
