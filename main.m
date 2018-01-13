% % Run annotation before this.
close all;
addpath('D:\allProjects\toolBox\toolbox_graph');
addpath('..\toolBox\xml_io_tools\');

% % % % % % paths
imP = 'data\FaceGen1\images\004.tiff';
annot = load('data\FaceGen1\images\004_annot.mat');     % annots are obtained using "annotation.m"
plyP = 'data\FaceGen1\meshLab\model_trim.ply';
pickpointP = '.\data\FaceGen1\meshLab\model_picked_points.pp';
fovList = load('data\FaceGen1\images\fovList.mat');
fov = fovList.fovList(4);

% imP = 'data\Jiedong\images\splice\2\DSC_1415.JPG';
% annot = load('data\Jiedong\images\splice\2\DSC_1415_annot.mat');
% plyP = 'data\Pengbo\face_coarse.ply';
% pickpointP = 'data\Pengbo\face_picked_points.pp';
% fov = 44.5;  % GT

% imP = 'data\Jiedong\images\splice\DSC_1415_2.JPG';
% annot = load('data\Jiedong\images\splice\DSC_1415_2_annot.mat');
% plyP = 'data\Jiedong\face_coarse.ply';
% pickpointP = 'data\Jiedong\face_picked_points.pp';
% fov = 44.5;  % GT

% imP = 'data\Jiedong\images\recap\DSC_1406_smooth.JPG';
% annot = load('data\Jiedong\images\recap\DSC_1406_smooth_annot.mat');
% plyP = 'data\Jiedong\face_coarse.ply';
% pickpointP = 'data\Jiedong\face_picked_points.pp';
% fov = 28.8;  % GT

% % parameters
Ns = 100;           % number of repeation for the whole estimation process. For uncertainty evaluation.
contourA = 75;      % angle for determining contour points. 75 degrees better for FaceGen model.
% contourA = 80;        % 80 degrees better for scanned model.
K_ICP =5;           % number of ICP iterations for contour refinement processs

% % show
im = imread(imP);
im = im(:,:,1:3);
x2d0 = double(annot.landmark);
x2d = x2d0(sum(x2d0, 2)~=0, :);
edgePts = annot.contour;
hImOverlay = figure();
imshow(im);
hold on;
plot(edgePts(:,1), edgePts(:,2), 'w.');
plot(x2d(:,1), x2d(:,2), 'k.', 'MarkerSize', 10);
plot(x2d(:,1), x2d(:,2), 'k*', 'MarkerSize', 5);


% % 3D landmarks
[shp, tl] = read_ply(plyP);
xml = xml_read(pickpointP);
X3d = zeros(size(xml.point, 1), 3); 
for i = 1:size(X3d, 1)
    X3d(i, 1) = xml.point(i).ATTRIBUTE.x;
    X3d(i, 2) = xml.point(i).ATTRIBUTE.y;
    X3d(i, 3) = xml.point(i).ATTRIBUTE.z;
end
X3d = X3d(sum(x2d0,2)~=0, :);
X3dH = [X3d, ones(size(X3d,1),1)];
figure();
DrawSolidHead(shp', tl', X3d');

% % Camera calibration
pAll = zeros(9,Ns);
for n = 1:Ns
    if n == 1
        x2d_n = x2d;  % the first time keep un-disturbed
    else
        x2d_n = x2d + sigma*randn(size(x2d,1), size(x2d,2));  % random sampling of 2D landmarks
    end
    x2dH_n = [x2d_n, ones(size(x2d_n, 1), 1)];      % Homogeneous Coordinates
    [P,K,R,t,rmse_n]=resectioning(X3dH, x2dH_n, 'off');
    [P,K,R,t,rmse_n,p]=constrainedCam(X3dH, x2dH_n, K,R,t,rmse_n, 'off');
    if n == 1
        sigma = rmse_n;
    end
    p = contourRefine(shp, tl, edgePts, X3d, x2d_n, p, K_ICP, contourA, hImOverlay, true);
    pAll(:,n) = p(:);
end
fThresh = max(size(im))/2/tand(12/2);  % seldom does a camera has FOV smaller than 12 degrees
fs = pAll(1,:);
pAll(:, fs>fThresh) = [];      % regard these as outliers
Ns = size(pAll, 2);

% % show
shpH = [shp, ones(size(shp,1),1)];
p = pAll(:,1);  % the un-disturbed estimation
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
f = max(size(im))/2/tand(fov/2);    % GT
plot3(size(im,2)/2, size(im,1)/2, f, 'ro');
plot3([size(im,2)/2,1], [size(im,1)/2,1], [f,0], 'r');
plot3([size(im,2)/2,size(im,2)], [size(im,1)/2,1], [f,0], 'r');
plot3([size(im,2)/2,1], [size(im,1)/2,size(im,1)], [f,0], 'r');
plot3([size(im,2)/2,size(im,2)], [size(im,1)/2,size(im,1)], [f,0], 'r');
plot3(pAll(2,:), pAll(3,:), pAll(1,:), 'bx');
plot3(pAll(2,1), pAll(3,1), pAll(1,1), 'rx');
xlabel('x'); ylabel('y'); zlabel('z');
axis on;

Nim = 1;
c = [size(im,2)/2, size(im,1)/2];
fEList = zeros(Ns, Nim);
cxEList = zeros(Ns, Nim);
cyEList = zeros(Ns, Nim);
for i = 1:Nim
    fList = max(size(im))/2/tand(fov/2);
    for n = 1:Ns
        p = pAll(:,n);
        fEList(n,i) = p(1);
        cxEList(n,i) = p(2);
        cyEList(n,i) = p(3);
    end
end
% % plot
figure('Name', 'Estimations');
subplot(1,3,1);
plot(1:Nim, fEList(1,:), 'rx');    % the one without perturbation    
hold on;
plot(1:Nim, fList, 'ro');
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

Sigma = cov(pAll(1:3,:)');
mu = mean(pAll(1:3,:), 2);
gt = [f, size(im,2)/2, size(im,1)/2]';
d = sqrt((gt-mu)'/Sigma*(gt-mu));
display(d);