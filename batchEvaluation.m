% % estimation resutls for all images in a folder. First use sampling
% % strategy, and then calculate the theoretical uncertainty.
% close all;
addpath('D:\allProjects\toolBox\toolbox_graph');
addpath('..\toolBox\xml_io_tools\');

% % paths
plyP = 'data\FaceGen1\meshLab\model_trim.ply';
pickpointP = 'data\FaceGen1\meshLab\model_picked_points.pp';
dataP = 'data\FaceGen1\images5\';
fovList = load(fullfile(dataP, 'fovList'));
fovList = fovList.fovList;
dirs = dir(dataP);
dirs = dirs(3:end); % skip './' and '../'
imList = {};
imListNE = {};
for i = 1:numel(dirs)
    idx = regexpi(dirs(i).name, '.tiff');
    if ~isempty(idx)
        imList = [imList; dirs(i).name];
        imListNE = [imListNE; dirs(i).name(1:idx-1)];
    end
end

Ns = 100;
Nim = numel(imList);
contourA = 75;      % angle for determining contour points
% contourA = 80;
K_ICP = 5;
show = false;

rmseAll = zeros(Nim, 1);
[shp, tl] = read_ply(plyP);
xml = xml_read(pickpointP);
X3d0 = zeros(size(xml.point, 1), 3); 
for i = 1:size(X3d0, 1)
    X3d0(i, 1) = xml.point(i).ATTRIBUTE.x;
    X3d0(i, 2) = xml.point(i).ATTRIBUTE.y;
    X3d0(i, 3) = xml.point(i).ATTRIBUTE.z;
end

pAll = zeros(9, Ns, Nim);
pAll_land = zeros(9, Ns, Nim);
fList = zeros(1, Nim*2);
fEList = zeros(Ns, Nim*2);
cList = zeros(2, Nim*2);
cxEList = zeros(Ns, Nim*2);
cyEList = zeros(Ns, Nim*2);
dList = zeros(Nim*2, 1);
dListR = [];
dListF = [];
for i = 1:Nim
    display([imList{i}, ' ...']);
    im = imread(fullfile(dataP, imList{i}));
    im = im(:,:,1:3);
    annot = load(fullfile(dataP, [imListNE{i}, '_annot.mat']));
%     annot = load(fullfile(dataP, [imListNE{i}, '_annot_partial.mat']));
    x2d0 = annot.landmark;
    x2d = x2d0(sum(x2d0, 2)~=0, :);
    edgePts = annot.contour;
    X3d = X3d0(sum(x2d0, 2)~=0, :);
    X3dH = [X3d, ones(size(X3d,1),1)];
    if show
        hImOverlay = figure();
        imshow(im);
        hold on;
        plot(x2d(:,1), x2d(:,2), 'r.');
        plot(edgePts(:,1), edgePts(:,2), 'w.');
    end

    pAll0 = zeros(9, Ns);
    pAll0_land = zeros(9, Ns);
    for n = 1:Ns
        if n == 1
            x2d_n = x2d;  % the first time keep un-disturbed
        else
%             x2d_n = x2d + sqrt(2)/2*sigma*randn(size(x2d,1), size(x2d,2));  % random sampling of 2D landmarks
            x2d_n = x2d + sigma*randn(size(x2d,1), size(x2d,2));  % random sampling of 2D landmarks
        end
        x2dH_n = [x2d_n, ones(size(x2d_n, 1), 1)];
        [P,K,R,t,rmse_n]=resectioning(X3dH, x2dH_n, 'off');
        [P,K,R,t,rmse_n,p]=constrainedCam(X3dH, x2dH_n, K,R,t,rmse_n, 'off');
        pAll0_land(:,n) = p(:);
        if n == 1
            rmse = rmse_n;
            sigma = rmse_n;
        end
        p = contourRefine(shp, tl, edgePts, X3d, x2d_n, p, K_ICP, contourA, 0, show);
        pAll0(:,n) = p(:);
    end
    pAll(:,:,i) = pAll0;
    pAll_land(:,:,i) = pAll0_land;
    rmseAll(i) = rmse;
    
    fov = fovList(i);
    fList(2*i-1) = max(size(im))/2/tand(fov/2);
    fList(2*i) = max(size(im))/2/tand(fov/2);
    cList(:, 2*i-1) = [size(im,2)/2, size(im,1)/2]';
    cList(:, 2*i) = [size(im,2)/2, size(im,1)/2]';
    for n = 1:Ns
        p = squeeze(pAll(:,n,i));
        p_land = squeeze(pAll_land(:,n,i));
        fEList(n,2*i-1) = p_land(1);     % before contour refinement
        cxEList(n,2*i-1) = p_land(2);
        cyEList(n,2*i-1) = p_land(3);
        fEList(n,2*i) = p(1);     % after contour refinement
        cxEList(n,2*i) = p(2);
        cyEList(n,2*i) = p(3);
    end
    
    fThresh = max(size(im))/2/tand(12/2);
    pAlli = squeeze(pAll_land(1:3,:,i));
    idx = pAlli(1,:)<fThresh;
    pAlli = pAlli(:, idx);
    Sigma = cov(transpose(pAlli));
    mu = mean(pAlli, 2);
    gt = [fList(2*i-1), size(im,2)/2, size(im,1)/2]';
    dList(2*i-1) = sqrt((gt-mu)'/Sigma*(gt-mu));
    pAlli = squeeze(pAll(1:3,:,i));
    idx = pAlli(1,:)<fThresh;
    pAlli = pAlli(:, idx);
    Sigma = cov(transpose(pAlli));
    mu = mean(pAlli, 2);
    gt = [fList(2*i), size(im,2)/2, size(im,1)/2]';
    dList(2*i) = sqrt((gt-mu)'/Sigma*(gt-mu));
    display(dList(2*i));
    
    pAlli = squeeze(pAll(1:3,:,i));
    idx = pAlli(1,:)<fThresh;
    pAlli = pAlli(:, idx);
    Sigma = cov(transpose(pAlli));
    mu = mean(pAlli, 2);
    for j = 1:Nim
        fov = fovList(j);
        f = max(size(im))/2/tand(fov/2);
        gt = [f, size(im,2)/2, size(im,1)/2]';
        d = sqrt((gt-mu)'/Sigma*(gt-mu));
        if j == i
            dListR = [dListR; d];
        else
            dListF = [dListF; d];
        end
    end
    
% %     save(fullfile(dataP, 'workspaceVars_partial.mat'));
end

vMax = max([dListR; dListF]);
edge = 0:0.5:vMax;
[countR, centerR] = hist(dListR, edge);
[countF, centerF] = hist(dListF, edge);

fs = 12;
ms = 7;
lw = 1;
figure;
plot(centerR, countR/numel(dListR), 'b', 'marker', 'x','LineWidth', lw, 'MarkerSize', ms); hold on;
plot(centerF, countF/numel(dListF), 'r', 'marker', 'o','LineWidth', lw, 'MarkerSize', ms);
ylabel('Percentage');
xlabel('Distance');
legend('Real', 'Fake')
title('Comparison of distance distributions');
grid on;
set(gca, 'fontsize', fs);

FAR = sum(dListF<max(dListR))/numel(dListF);
display(FAR);
% % ROC curve
label = zeros(Nim^2,1);
label(Nim+1:Nim^2) = 1;
[FA,DR,T,AUC] = perfcurve(label, [dListR; dListF], 1);       % false alarm and detection rate
display(AUC);
figure;
% res = numel(T)/50;
res = 1;
plot(FA(1:res:end), DR(1:res:end), 'b', 'LineWidth', lw);
xlabel('False Alarm Rate'); 
ylabel('Detection Rate');
title('ROC curve');
grid on;
set(gca, 'fontsize', fs);

figure('Name', 'f estimation');
x = 1:3*Nim;
x(3:3:3*Nim) = [];
hMLE0 = plot(x, fEList(1,:), 'rx', 'MarkerSize', ms);    % the one without perturbation    
hold on;
hGT = plot(x, fList, 'ro', 'MarkerSize', ms);
hBox1 = boxplot(fEList(:, 1:2:Nim*2), 'color', 'b', 'positions', x(1:2:Nim*2),...
    'labels', x(1:2:Nim*2), 'width', 0.7);
hBox2 = boxplot(fEList(:, 2:2:Nim*2), 'color', 'm', 'positions', x(2:2:Nim*2),...
    'labels', x(2:2:Nim*2), 'width', 0.7);
set(hBox1(:), 'linewidth', lw);
set(hBox2(:), 'linewidth', lw);
delete(findobj(gca,'tag','Outliers'));
set(gca, 'xtick', x);
set(gca, 'xticklabel', ceil(x/3));
grid on;
legend([hMLE0, hGT, hBox1(1), hBox2(1)], 'MLE0', 'GT', 'No refinement', 'Refinement');
title('f Estimation');
xlabel('Image ID');
ylabel('f /pixels');
set(gca, 'fontsize', fs);
axis auto;

figure('Name', 'c_x estimation');
x = 1:3*Nim;
x(3:3:3*Nim) = [];
hMLE0 = plot(x, cxEList(1,:), 'rx', 'MarkerSize', ms);    % the one without perturbation    
hold on;
hGT = plot(x, cList(1,:), 'ro', 'MarkerSize', ms);
hBox1 = boxplot(cxEList(:, 1:2:Nim*2), 'color', 'b', 'positions', x(1:2:Nim*2),...
    'labels', x(1:2:Nim*2), 'width', 0.7);
hBox2 = boxplot(cxEList(:, 2:2:Nim*2), 'color', 'm', 'positions', x(2:2:Nim*2),...
    'labels', x(2:2:Nim*2), 'width', 0.7);
set(hBox1(:), 'linewidth', lw);
set(hBox2(:), 'linewidth', lw);
delete(findobj(gca,'tag','Outliers'));
set(gca, 'xtick', x);
set(gca, 'xticklabel', ceil(x/3));
grid on;
legend([hMLE0, hGT, hBox1(1), hBox2(1)], 'MLE0', 'GT', 'No refinement', 'Refinement');
title('c_x Estimation');
xlabel('Image ID');
ylabel('c_x /pixels');
set(gca, 'fontsize', fs);
axis auto;

figure('Name', 'c_y estimation');
x = 1:3*Nim;
x(3:3:3*Nim) = [];
hMLE0 = plot(x, cyEList(1,:), 'rx', 'MarkerSize', ms);    % the one without perturbation    
hold on;
hGT = plot(x, cList(2,:), 'ro', 'MarkerSize', ms);
hBox1 = boxplot(cyEList(:, 1:2:Nim*2), 'color', 'b', 'positions', x(1:2:Nim*2),...
    'labels', x(1:2:Nim*2), 'width', 0.7);
hBox2 = boxplot(cyEList(:, 2:2:Nim*2), 'color', 'm', 'positions', x(2:2:Nim*2),...
    'labels', x(2:2:Nim*2), 'width', 0.7);
set(hBox1(:), 'linewidth', lw);
set(hBox2(:), 'linewidth', lw);
delete(findobj(gca,'tag','Outliers'));
set(gca, 'xtick', x);
set(gca, 'xticklabel', ceil(x/3));
grid on;
legend([hMLE0, hGT, hBox1(1), hBox2(1)], 'MLE0', 'GT', 'No refinement', 'Refinement');
title('c_y Estimation');
xlabel('Image ID');
ylabel('c_y /pixels');
set(gca, 'fontsize', fs);
axis auto;

% % theoritical covariance
% p_sigmaAll = zeros(9,9,Nim);
% for i = 1:Nim
%     idx = sprintf('%03d', i);
%     display([idx, '.tiff ...']);
%     sigma = rmseAll(i);
% %     p_cov = covBackProp(squeeze(pAll_land(:,1,i)), X3d, 1/2*sigma^2*eye(2*size(X3d,1)));
%     p_cov = covBackProp(squeeze(pAll_land(:,1,i)), X3d, sigma^2*eye(2*size(X3d,1)));
%     p_sigma = sqrt(abs(p_cov));
%     p_sigmaAll(:,:,i) = p_sigma;
% end
% 
% figure;
% plot(1:Nim, squeeze(p_sigmaAll(1,1,:)), 'r-o', 'LineWidth', lw, 'MarkerSize', ms);
% hold on;
% plot(1:Nim, std(fEList(:, 1:2:Nim*2), 0, 1), 'b-x', 'LineWidth', lw, 'MarkerSize', ms);
% set(gca, 'xtick', 1:Nim);
% set(gca, 'xticklabel', 1:Nim);
% grid on;
% title('std of MLE for f');
% xlabel('Image ID');
% ylabel('std /pixels');
% legend('Theoretical', 'Random perturbation');
% set(gca, 'fontsize', fs);
% 
% figure;
% plot(1:Nim, squeeze(p_sigmaAll(2,2,:)), 'r-o', 'LineWidth', lw, 'MarkerSize', ms);
% hold on;
% plot(1:Nim, std(cxEList(:, 1:2:Nim*2), 0, 1), 'b-x', 'LineWidth', lw, 'MarkerSize', ms);
% set(gca, 'xtick', 1:Nim);
% set(gca, 'xticklabel', 1:Nim);
% grid on;
% title('std of MLE for c_x');
% xlabel('Image ID');
% ylabel('std /pixels');
% legend('Theoretical', 'Random perturbation');
% set(gca, 'fontsize', fs);
% 
% figure;
% plot(1:Nim, squeeze(p_sigmaAll(3,3,:)), 'r-o', 'LineWidth', lw, 'MarkerSize', ms);
% hold on;
% plot(1:Nim, std(cyEList(:, 1:2:Nim*2), 0, 1), 'b-x', 'LineWidth', lw, 'MarkerSize', ms);
% set(gca, 'xtick', 1:Nim);
% set(gca, 'xticklabel', 1:Nim);
% grid on;
% title('std of MLE for c_y');
% xlabel('Image ID');
% ylabel('std /pixels');
% legend('Theoretical', 'Random perturbation');
% set(gca, 'fontsize', fs);