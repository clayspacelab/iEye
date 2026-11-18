function [rsvp_params] = rsvp_prf_loadparams(params_file)
%RSVP_PRF_LOADPARAMS Read task parameters from an RSVP pRF rsvp_params.txt file.
%
%   rsvp_params = rsvp_prf_loadparams(params_file)
%
%   Each parameter in the file is a comment line naming it, followed by its
%   value on the next line, then a separator line of dashes ending in '#':
%       # tr = TR length. Format: float
%       1.3
%       ---------#
%   The name is the text between '#' and the first '='. The number, order,
%   and names of the parameters are taken from the file. A parameter whose
%   value line is missing or blank gets []. Values are converted to
%   numeric/logical/etc. by rsvp_parse_value.
%
%   ppd (pixels per degree) is computed from screen_res, screen_height, and
%   view_dist and added to the output.
%
%   params_file defaults to rsvp_params.txt in the current folder.
%
%   See also rsvp_parse_value, rsvp_read_params.


%look for rsvp_params.txt in cwd by default
if ~exist('params_file','var')
    params_file = 'rsvp_params.txt';
end

% Load parameters from RSVP pRF task
lines = strip(readlines(params_file));
is_name = startsWith(lines, "#") & contains(lines, "=");
is_sep  = matches(lines, asManyOfPattern("-", 1) + "#");

% name is between '#' and the first '='; value is the line after the name,
% unless that line is a separator or another parameter (i.e., value missing)
name_idx = find(is_name);
names    = strip(extractBefore(extractAfter(lines(name_idx), "#"), "="));
val_idx  = name_idx + 1;
has_val  = val_idx <= numel(lines);
has_val(has_val) = ~is_sep(val_idx(has_val)) & ~is_name(val_idx(has_val));
values   = repmat("", size(names));
values(has_val) = lines(val_idx(has_val));

% convert values to numeric/logical/etc. (see rsvp_parse_value)
names  = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));
values = arrayfun(@rsvp_parse_value, values, 'UniformOutput', false);
rsvp_params = cell2struct(values, cellstr(names), 1);

%fix size is specified as a percentage in rsvp_params and output as a
%decimal proportion in summary csv, so harmonizing here
rsvp_params.fix_size = rsvp_params.fix_size/100;

% Get screen params, convert to DVA
if all(isfield(rsvp_params, {'screen_res','screen_height','view_dist'}))
    rsvp_params.screen_width = rsvp_params.screen_height * rsvp_params.screen_res(1) / rsvp_params.screen_res(2);  
    rsvp_params.screen_h_deg = (180 / pi) * 2 * atan(rsvp_params.screen_height / (2 * rsvp_params.view_dist));
    rsvp_params.screen_w_deg = (180 / pi) * 2 * atan(rsvp_params.screen_width / (2 * rsvp_params.view_dist));
    rsvp_params.ppd  = rsvp_params.screen_res(2) / rsvp_params.screen_h_deg;
else
    warning('Screen parameters missing from %s; ppd not computed.', params_file)
end

rsvp_params.task = 'rsvp';
rsvp_params.filename = params_file;


end
