function p = contourRefine(shp, tl, edgePts, X3d, x2d, p, K_ICP, contourA, hImOverlay, show)
% % Use contour observations to refine the parameter estimation
% % input shp: nx3, 3D model points
% % input tl: nx3, 3D triangle list
% % input edgePts: nx2
% % input X3d: nx3, 3D landmarks
% % input x2d: nx2, 2D landmarks
% % input p: 9x1
% % input K_ICP: number of iterations
% % input contourA: angle for determining contour points
dThresh = 0.03*range(x2d(:,2));    % threshold for outlier matching
w0 = 0.1;  % the weight ratio of landmarks to contour points. normalized by number of points.
OptDisplay = 'off';
for k = 1:K_ICP
    ea = p(4:6);
    R = rotm('xyz', ea);
    t = p(7:9);
    shpC = (R*shp'+repmat(t,1,size(shp,1)))';
    normalC = compute_normal(shpC', tl');
    normalC = normalC';
    cosine = zeros(size(shpC,1),1);
    for i = 1:size(shpC,1)
        cosine(i) = shpC(i,:)/norm(shpC(i,:))*normalC(i,:)';
    end
    idxProjC = cosine<cosd(contourA)&cosine>0;
    
    P = para2Proj(p);
    shpH = [shp, ones(size(shp,1),1)];
    shpP = (P*shpH')';   % project model
    shpP = shpP./repmat(shpP(:,3), 1, 3);
    shpP = shpP(:, 1:2);
%     temp = boundary(shpP(:,1), shpP(:,2));
%     idxBoundary = zeros(size(shpP,1),1);
%     idxBoundary(temp) = 1;
%     idxBoundary = logical(idxBoundary);
%     idx_CB = idxProjC|idxBoundary;
    idx_CB = idxProjC;
    
    projC2d = shpP(idx_CB,:);
    projC3d = shp(idx_CB,:);
%     faces = delaunay(edgePts(:,1), edgePts(:,2));
%     [idx, d] = dsearchn(edgePts, faces, projC2d); 
%     idxOutlier = d>dThresh;
%     x2d_edge = edgePts(idx(~idxOutlier),:);
%     X3d_edge = projC3d(~idxOutlier,:);

    warning('off', 'MATLAB:delaunay:DupPtsDelaunayWarnId');
    faces = delaunay(projC2d(:,1), projC2d(:,2));
    [idx, d] = dsearchn(projC2d, faces, edgePts);
    idxOutlier = d>dThresh;
    x2d_edge = edgePts(~idxOutlier,:);
    X3d_edge = projC3d(idx(~idxOutlier),:);
    
    if show
        X3dH = [X3d, ones(size(X3d,1),1)];
        x2d_hat = (P*X3dH')';   % project 3D landmars
        x2d_hat = x2d_hat./repmat(x2d_hat(:,3), 1, 3);
        figure(hImOverlay);
        if exist('hProj', 'var')
            delete([hProj, hProjC, hProjL, hProjC_outlier]);
        end
        hProj = plot(shpP(:,1), shpP(:,2), 'y.');
        hProjL = plot(x2d_hat(:,1), x2d_hat(:,2), 'g.', 'MarkerSize', 10); 
        hProjC = plot(projC2d(:,1), projC2d(:,2), 'k.', 'MarkerSize', 8);
%         hProjC_outlier = plot(projC2d(idxOutlier,1), projC2d(idxOutlier,2), 'r.');
        hProjC_outlier = plot(edgePts(idxOutlier,1), edgePts(idxOutlier,2), 'r.');
        drawnow;
%         display(p);
    end
   
    X3dH = [X3d, ones(size(X3d,1),1)];
    x2dH = [x2d, ones(size(x2d,1),1)];
    X3d_edgeH = [X3d_edge, ones(size(X3d_edge,1),1)];
    x2d_edgeH = [x2d_edge, ones(size(x2d_edge,1),1)];
    w = w0*size(X3d_edge,1)/size(X3d,1);
    options = optimset('Algorithm','levenberg-marquardt', 'Display', OptDisplay);
    p = lsqnonlin(@(p) errProjectionNonlin3(p,x2d_edgeH,X3d_edgeH,x2dH,X3dH,w), p,[],[],options);
end
if exist('hProj', 'var')
    delete([hProj, hProjC, hProjL, hProjC_outlier]);
end