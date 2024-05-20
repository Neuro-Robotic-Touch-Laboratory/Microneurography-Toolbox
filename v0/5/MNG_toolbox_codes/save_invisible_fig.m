function save_invisible_fig(h, filename)
%SAVE_INVISIBLE_FIG Summary of this function goes here
%   Detailed explanation goes here
set(h, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 

savefig(h, 'filename')
end

