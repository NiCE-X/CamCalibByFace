% % This script for anotating landmarks and contours
% % The inside landmarks are detected by face alignment algorithms (here, I
% % used IntraFace toolkit). You can use any alignment methods you have.
% The landmarks on left ear top, left ear bottom, right ear top, right ear
% % bottom and chin are annotated by pressing 1,2,3,4,5 respectively on the main keyboard (The
% % small keyboard area of numbers may not work). Then contours are annotated
% % by mouse drag.

close all;
imP = 'data\Jiedong\images\splice\DSC_1415_2.JPG';

global landmark;
global contour;
landmark = zeros(24,2);
contour = [];
baseP = pwd();
im = imread(imP);
im = im(:,:,1:3);
cd('.\intraFace\');
[DM,TM,option] = xx_initialize;
option.min_neighbors = 5;
[pred,pose] = xx_track_detect(DM,TM,im,[],option);
pred = double(pred);
cd(baseP);
assert(~isempty(pred));
% % index2 = 1:49;
index2 = [1 3 5 6 8 10 11:14 15 17 19 20 23 26 29 32 38];
pred = pred(index2,:);
hAnnot = figure();
imshow(im);
hold on;
plot(pred(:,1), pred(:,2), 'r.');
landmark(1:19,:) = pred;
set(hAnnot, 'WindowKeyPressFcn', @mypresskeycallback_hAnnot);
set(hAnnot, 'WindowButtonDownFcn',@mystartDragFcn);
set(hAnnot, 'WindowButtonUpFcn',@mystopDragFcn);