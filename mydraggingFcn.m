function mydraggingFcn(varargin)
    pt = get(get(gcbf, 'CurrentAxes'), 'CurrentPoint');
    pt = pt(1, 1:2);    % get CurrentPoint returns a 2 x 3 matrix
    global contour;
    contour = [contour; pt];
    plot(pt(end,1), pt(end,2), 'b.', 'ButtonDownFcn',@mystartDragFcn);
    hold on;
    drawnow ;
end