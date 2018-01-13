function mypresskeycallback_hAnnot(obj,evd)
global landmark;
if (strcmp(evd.Key, '1'))   % left side ear top (photo's left)
    [x,y] = ginput(1);
    landmark(20,:) = [x,y];
    plot(x, y, 'r.');
elseif (strcmp(evd.Key, '2'))   % left side ear bottom
    [x,y] = ginput(1);
    landmark(21,:) = [x,y];
    plot(x, y, 'r.');
elseif (strcmp(evd.Key, '3'))   % right side ear top
    [x,y] = ginput(1);
    landmark(22,:) = [x,y];
    plot(x, y, 'r.');
elseif (strcmp(evd.Key, '4'))   % right side ear bottom
    [x,y] = ginput(1);
    landmark(23,:) = [x,y];
    plot(x, y, 'r.');
elseif (strcmp(evd.Key, '5'))   % chin
    [x,y] = ginput(1);
    landmark(24,:) = [x,y];
    plot(x, y, 'r.');
end
end