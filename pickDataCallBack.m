function output_txt = pickDataCallBack(obj,event_obj)
% Display nothing. Use this func instead of the default display to save
% space when running pickLandmarks.m
pos = get(event_obj,'Position');
output_txt = {['']};
