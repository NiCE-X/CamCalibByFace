function pbrt4dataset(tempP, pbrtFile)
% % This function generate synthetic dataset using pbrt
baseP = pwd();
addpath(baseP);
cd(tempP);
[~, ~] = system(['obj2pbrt model.obj model.pbrt']);    % convert obj to pbrt file.
[~, ~] = system(['copy /Y ', pbrtFile, ' dataset.pbrt']);
saveP = 'images';    % folder name to save the rendered images.
[~, ~] = system(['md ', saveP]);
fovList = 20:10:70;
save(fullfile(saveP, 'fovList.mat'), 'fovList');
% xEyes = [150, 100];  % the wanted pixel distance between eyes
xEyes = [150];
XEyes = 69; % the 3D distance between eyes of the model
content = readFile('dataset.pbrt');
N = numel(fovList)*numel(xEyes);
count = 0;
for i = 1:numel(xEyes)
    for j = 1:numel(fovList);
        count = count+1;
        idx = sprintf('%03d', count);
        s = sprintf('rendering %d/%d ...', count, N);
        disp(s);
        d = D_fov(fovList(j), xEyes(i), XEyes);
        content{9} = ['LookAt 0 0 ', num2str(d), '  0 0 0  0 1 0'];
        content{11} = ['Camera "perspective" "float fov" [', num2str(fovList(j)), ']'];
        writeFile('dataset.pbrt', content);
        % execute the altered pbrt file and create tiff images.
        [~, ~] = system(['pbrt dataset.pbrt']);
        s = ['exrtotiff ./', saveP, '/001.exr ./', saveP, '/', idx, '.tiff'];
        [~, ~] = system(s);
    end
end
cd(baseP);
end