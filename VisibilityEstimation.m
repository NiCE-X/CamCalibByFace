function [visibility] = VisibilityEstimation(ProjectShape, tri, opentime)
% % input ProjectShape, tri: 3xn, 3xm
% % input opentime: int, visibility close times
% % input threshold: don't need, already self-adapted
ProjectShape(1,:) = ProjectShape(1,:) - min(ProjectShape(1,:));
ProjectShape(2,:) = ProjectShape(2,:) - min(ProjectShape(2,:));

width = max(ProjectShape(1,:)) - min (ProjectShape(1,:));
height = max(ProjectShape(2,:)) - min (ProjectShape(2,:));

temp_width = width;
temp_height = height;
temp = randperm(size(tri, 2), 1000);
avrD = 0;
for i = 1:numel(temp)
   avrD = avrD + norm(ProjectShape(:, tri(1, temp(i))) - ProjectShape(:, tri(2, temp(i))))/numel(temp);
end
targSize = 5;
scale = targSize/avrD;     % scale average triangle size to targSize

visibility = Mex_OcclusionDetection(ProjectShape, tri, temp_width, temp_height, 1, scale);      % the threshold is no longer needed, here is 1, any value doesn't matter
visibility = logical(visibility);

for i = 1:opentime         % pb: first dilate, then erode = close operation.
    visibility = MeshDilate(visibility, tri);
end

for i = 1:opentime        % to erode one more time, because the projected points on face profile often includes background
    visibility = MeshErode(visibility, tri);
end

end


function [visibility] = MeshDilate(visibility, tri)     %pb: visiblity propagate. the other two points of a trangle is made visible if one point of it is visible. kind of visibility dilate.
    vis_tri_bin = visibility(tri(1,:)) | visibility(tri(2,:)) | visibility(tri(3,:)); %pb: this coding style is great. I would have used for loops...
    ind = tri(:,vis_tri_bin);
    visibility(ind) = 1;
end

function [visibility] = MeshErode(visibility, tri)      %pb: invisiblity propagate.
    vis_tri_bin = visibility(tri(1,:)) & visibility(tri(2,:)) & visibility(tri(3,:));
    ind = tri(:,~vis_tri_bin);
    visibility(ind) = 0;
end