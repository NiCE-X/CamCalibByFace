function pbrt4dataset2(tempP, pbrtFile)
% % This function generate synthetic dataset.
baseP = pwd();
addpath(baseP);
cd(tempP);
[~, ~] = system(['obj2pbrt model.obj model.pbrt']);    % convert obj to pbrt file.
[~, ~] = system(['copy /Y ', pbrtFile, ' dataset.pbrt']);
saveP = 'temp';    % folder name to save the rendered images.
[~, ~] = system(['md ', saveP]);
N = 1; % number of images to create
fovList = 20+50*rand(N,1);
rtList = zeros(N, 5); % [ry, rx, tx, ty, tz]
save(fullfile(saveP, 'fovList.mat'), 'fovList');
% xEyes = [150, 100];  % the wanted pixel distance between eyes
XEyes = 69; % the 3D distance between eyes of the model
content = readFile('dataset.pbrt');
count = 0;
for j = 1:numel(fovList);
    count = count+1;
    idx = sprintf('%03d', count);
    s = sprintf('rendering %d/%d ...', count, N);
    disp(s);
    xEyes = 100+50*rand(1);
    d = D_fov(fovList(j), xEyes, XEyes);
    x = -512+1.5*xEyes+(1024-3*xEyes)*rand(1);
    y = -512+2*xEyes+(1024-4*xEyes)*rand(1);
    f = 1024/2/tand(fovList(j)/2);
    X = d*x/f;
    Y = d*y/f;
    ry = -15+30*rand(1);
    rx = -10+20*rand(1);
    rtList(j,:) = [ry, rx, X, Y, d];
    content{9} = ['LookAt 0 0 ', num2str(d), '  0 0 0  0 1 0'];
    content{11} = ['Camera "perspective" "float fov" [', num2str(fovList(j)), ']'];
    content{22} = ['Translate ', num2str(X), ' ', num2str(Y), ' 0'];
    content{24} = ['Rotate ', num2str(ry), ' 0 1 0'];
    content{26} = ['Rotate ', num2str(rx), ' 1 0 0'];
    writeFile('dataset.pbrt', content);
    % execute the altered pbrt file and create tiff images.
    [~, ~] = system(['pbrt dataset.pbrt']);
    s = ['exrtotiff ./', saveP, '/001.exr ./', saveP, '/', idx, '.tiff'];
    [~, ~] = system(s);
end
save(fullfile(saveP, 'rtList.mat'), 'rtList');
cd(baseP);
end