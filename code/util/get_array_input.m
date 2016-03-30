function array_out = get_array_input(array_in,dflt)
% out_array = get_array_input(array_in,dflt)
%
% Utility function to return a cell array given an input which can be a
% cell array or a .grp file
%
% Input: 
%       array_in: either a .grp or cell array of inputs
%       dflt: a default value to return if array_in is empty
%
% Output:
%       array_out: a cell array

if isempty(array_in)
    array_out = dflt;
elseif iscell(array_in)
    array_out = array_in;
elseif exist(array_in,'file')
    array_out = parse_grp(array_in);
else
    error('Unknown input type')
end

end

