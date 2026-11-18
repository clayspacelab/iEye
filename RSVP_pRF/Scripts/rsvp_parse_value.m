function v = rsvp_parse_value(s)
%RSVP_PARSE_VALUE Convert one parameter value written by the RSVP task.
%
%   v = rsvp_parse_value(s)
%
%   s is a char or string value as written by the (Python) task, e.g. in the
%   summary CSV or rsvp_params.txt. It is converted with jsondecode after
%   mapping Python literals (True/False/None, single quotes) to JSON. A list
%   without brackets, as in rsvp_params.txt, is bracketed first:
%       '1.3'               -> 1.3
%       'True' / 'False'    -> true / false
%       '[1920, 1080]'      -> [1920 1080]   (row vector)
%       '1920,1080'         -> [1920 1080]
%       '[[1, 2], [3, 4]]'  -> [1 2; 3 4]
%       'None' or empty     -> []
%       anything else       -> char, e.g. '192.168.1.5', 'pix'
%
%   See also rsvp_read_params, rsvp_prf_loadparams, jsondecode.

    s = strip(string(s));
    if s == ""
        v = [];
        return
    end
    json = regexprep(s, {'\<True\>', '\<False\>', '\<None\>', ''''}, ...
                        {'true',     'false',     'null',     '"'});
    if contains(json, ",") && ~startsWith(json, "[")
        json = "[" + json + "]";     % bare list, e.g. 1920,1080
    end
    try
        v = jsondecode(json);
        if (isnumeric(v) || islogical(v)) && isvector(v)
            v = reshape(v, 1, []);   % jsondecode returns columns
        end
    catch
        v = char(s);                 % not a literal: keep the text
    end
end
