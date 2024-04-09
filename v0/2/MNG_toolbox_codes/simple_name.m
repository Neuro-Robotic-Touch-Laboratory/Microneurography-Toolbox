function output_str = simple_name(input_str)
%SIMPLE_NAME removes special caracters from the template interval name 
%   Detailed explanation goes here

tmp_idx = strfind(input_str, '<$§');
if ~isempty(tmp_idx)
    input_str(tmp_idx:end) = [];
end
tmp_idx = strfind(input_str,'§$>');
if ~isempty(tmp_idx)
    input_str(1:3) = [];
end
tmp_idx = strfind(input_str, '§$');
if ~isempty(tmp_idx)
    input_str(tmp_idx:tmp_idx+3) = [];
end
tmp_idx = strfind(input_str, '$§');
if ~isempty(tmp_idx)
    input_str(tmp_idx-2:tmp_idx+1) = [];
end
output_str = input_str;
end