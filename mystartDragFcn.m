function mystartDragFcn(varargin)
    set(gcbf, 'WindowButtonMotionFcn', @mydraggingFcn );
    pt = get(get(gcbf, 'CurrentAxes'), 'CurrentPoint');
    pt = pt(1, 1:2);    % get CurrentPoint returns a 2 x 3 matrix
    global countour;
    countour = [countour; pt];
end