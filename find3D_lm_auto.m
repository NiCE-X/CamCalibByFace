% % This file finds 3D landmarks automatically using the synthetic image of
% % the 3D model. Depend on back projection of 2D landmarks.
close all;
addpath('Mex\');
addpath('D:\allProjects\toolBox\toolbox_graph');
baseP = pwd();
im = imread('data\scan1\001.tiff');
% im = imread('data\scan1\images\PS.JPG');
im = im(:,:,1:3);
cd('.\intraFace\');
[DM,TM,option] = xx_initialize;
[pred,pose] = xx_track_detect(DM,TM,im,[],option);
cd(baseP);
assert(~isempty(pred));
% % index2 = [20 23 26 29 15 17 19 32 38]';
% index2=[11:14 15 17 19 20:31 32 35 38 41];
index2 = 1:49;
hf = figure();
imshow(im);
hold on;
plot(pred(index2,1), pred(index2,2), 'r.');
for i = 1:5
    [x,y] = ginput(1);
    pred = [pred; [x,y]];
    plot(x, y, 'r.');
end

[shp, tl] = read_ply('.\data\scan1\MeshedReconstruction.ply');
R_World2Cam = [1 0 0;
                0 -1 0;
                0 0 -1];
R_Object2World = rotm('yxz', [10 20 0]);
fov = 20;
t = [0 0 0.8]';
cx = size(im, 2)/2;
cy = size(im, 1)/2;
f = max([cx,cy])/tand(fov/2);
K = [f 0 cx;
     0 f cy
     0 0 1];

shpR = R_Object2World*shp';
shpR = shpR';
opentime = 2;
% [visibility] = VisibilityEstimation(shpR', tl', opentime);    % pb: visibility of each vertex
visibility = logical(ones(size(shpR,1), 1));
shpP = K*(R_World2Cam*R_Object2World*shp'+repmat(t, 1, size(shp, 1)));
shpP = shpP';
shpP = shpP./repmat(shpP(:,3), 1, 3);
shpP = shpP(:, 1:2);
figure(hf);
shpPV = shpP(visibility, :);
plot(shpPV(:,1), shpPV(:,2), '.');
faces = delaunay(shpPV(:,1), shpPV(:,2));
[idx, ~] = dsearchn(shpPV, faces, double(pred));   % closest point search for each 2D landmark
plot(shpPV(idx,1), shpPV(idx,2), 'g.');

temp = cumsum(visibility);
indLandMarks = zeros(size(idx));
for i = 1:numel(idx)
    indLandMarks(i) = find(temp == idx(i), 1);
end
figure;
DrawSolidHead(shp', tl', transpose(shp(indLandMarks, :)));