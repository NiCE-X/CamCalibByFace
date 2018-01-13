% % show the result of one image from a batchEvaluation.m result
% % First, you need to load e.g. 'data\FaceGen1\images5\workspaceVars.mat'
% % to the workspace.
% close all;

for i = 1
    im = imread(fullfile(dataP, imList{i}));
    im = im(:,:,1:3);
    annot = load(fullfile(dataP, [imListNE{i}, '_annot.mat']));
    x2d0 = annot.landmark;
    x2d = x2d0(sum(x2d0, 2)~=0, :);
    edgePts = annot.contour;

    % M distance
    fThresh = max(size(im))/2/tand(12/2);
    pAlli = squeeze(pAll(1:3,:,i));
    idx = pAlli(1,:)<fThresh;
    pAlli = pAlli(:, idx);
    Sigma = cov(transpose(pAlli));
    mu = mean(pAlli, 2);
    gt = [fList(2*i), size(im,2)/2, size(im,1)/2]';
    d = sqrt((gt-mu)'/Sigma*(gt-mu));
    display(d);
    
    % % show
    hImOverlay = figure();
    imshow(im);
    hold on;
    plot(edgePts(:,1), edgePts(:,2), 'w.');
    plot(x2d(:,1), x2d(:,2), 'r.', 'MarkerSize', 10);

    shpH = [shp, ones(size(shp,1),1)];
    pAlli = pAll(:,:,i);
    pAlli = pAlli(:, idx);
    p = pAlli(:,1);  % the un-disturbed estimation
    P = para2Proj(p);
    shpP = (P*shpH')';   % project model
    shpP = shpP./repmat(shpP(:,3), 1, 3);
    shpP = shpP(:, 1:2);
    figure(hImOverlay);
    plot(shpP(:,1), shpP(:,2), 'y.');
    X3dH = [X3d, ones(size(X3d,1),1)];
    x2d_hat = (P*X3dH')';   % project 3D landmars
    x2d_hat = x2d_hat./repmat(x2d_hat(:,3), 1, 3);
    plot(x2d_hat(:,1), x2d_hat(:,2), 'g.');
    plot3(size(im,2)/2, size(im,1)/2, fList(2*i), 'ro');
    plot3([size(im,2)/2,1], [size(im,1)/2,1], [f,0], 'r');
    plot3([size(im,2)/2,size(im,2)], [size(im,1)/2,1], [f,0], 'r');
    plot3([size(im,2)/2,1], [size(im,1)/2,size(im,1)], [f,0], 'r');
    plot3([size(im,2)/2,size(im,2)], [size(im,1)/2,size(im,1)], [f,0], 'r');
    plot3(pAlli(2,:), pAlli(3,:), pAlli(1,:), 'bx');
%     plot3(pAlli(2,1), pAlli(3,1), pAlli(1,1), 'rx');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis on;

    Nim = 1;
    c = [size(im,2)/2, size(im,1)/2];
    fEList = zeros(Ns, Nim);
    cxEList = zeros(Ns, Nim);
    cyEList = zeros(Ns, Nim);
    for j = 1:Nim
        for n = 1:size(pAlli,2)
            p = pAlli(:,n);
            fEList(n,j) = p(1);
            cxEList(n,j) = p(2);
            cyEList(n,j) = p(3);
        end
    end
    % % plot
    figure('Name', 'Estimations');
    subplot(1,3,1);
    plot(1:Nim, fEList(1,:), 'rx');    % the one without perturbation    
    hold on;
    plot(1:Nim, fList(2*i), 'ro');
    boxplot(fEList);
    grid on;
    legend('MLE0', 'GT')
    title('f Estimation');
    xlabel('Image ID');
    ylabel('f /pixel');
    axis auto;
    title('f Estimation');

    subplot(1,3,2);
    plot(1:Nim, cxEList(1,:), 'rx');    % the one without perturbation
    hold on;
    plot(1:Nim, c(1)*ones(Nim,1), 'ro');
    boxplot(cxEList);
    grid on;
    legend('MLE0', 'GT')
    title('c_x Estimation');
    xlabel('Image ID');
    ylabel('c_x /pixel');
    axis auto;
    title('c_x Estimation');

    subplot(1,3,3);
    plot(1:Nim, cyEList(1,:), 'rx');    % the one without perturbation
    hold on;
    plot(1:Nim, c(2)*ones(Nim,1), 'ro');
    boxplot(cyEList);
    grid on;
    legend('MLE0', 'GT')
    title('c_y Estimation');
    xlabel('Image ID');
    ylabel('c_y /pixel');
    axis auto;
    title('c_y Estimation');
end