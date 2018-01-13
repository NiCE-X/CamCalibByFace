% % pick landmark points on a mesh plot, and find their indices
%%
[shp, tl] = read_ply('.\data\scan1\MeshedReconstruction5_trim.ply');
figure;
% trimesh(tl, shp(:,1), shp(:,2), shp(:,3), 'EdgeColor', 'r');
DrawSolidHead(shp', tl');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal vis3d;
%%
% %% output vertex indices. All landmarks have been exported in batch. Set
% % the pickdata call back to 'pickDataCallBack.m' to not display annoying
% % flags.
indLandMarks = [];
for i = numel(cursor_info):-1:1     % data cursor is in reverse order
    temp = findVertex(shp, cursor_info(i).Position);
    indLandMarks = [indLandMarks; temp];
end
figure;
DrawSolidHead(shp', tl');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal vis3d;
hold on;
plot3(shp(indLandMarks, 1), shp(indLandMarks, 2), shp(indLandMarks, 3), 'r*', 'MarkerSize', 10);

%% ouput vertex coordinates. All landmarks have been exported in batch.
% landmarks = zeros(numel(cursor_info), 3);
% for i = 1:numel(cursor_info)
%     landmarks(i,:) = cursor_info(i).Position;
% end
% landmarks = landmarks(end:-1:1, :);     % data cursor is in reverse order
% figure;
% DrawSolidHead(shp', tl');
% hold on;
% plot3(landmarks(:,1), landmarks(:,2), landmarks(:,3), 'r*', 'MarkerSize', 10);