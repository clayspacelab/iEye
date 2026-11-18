function params = rsvp_read_params(filename)
%RSVP_READ_PARAMS Read the name/value parameters below the table in a summary CSV.
%
%   params = rsvp_read_params(filename)
%
%   Every non-blank line after the first blank line that follows the table
%   is read as "name,value". The value is everything after the first comma,
%   so unquoted lists like "background_color,[0, 0, 0]" stay intact. Lines
%   of only commas/whitespace count as blank.
%
%   NOTE: screen_h_deg, screen_w_deg, and ppd params are recalculated even
%   if present due to the existence of a known bug in the output of past
%   versions of the RSVP task.
%
%   Values are converted to numeric, logical, etc. by rsvp_parse_value.
%
%   Requires R2020b or later (readlines).
%
%   See also rsvp_read_summary_table, rsvp_parse_value.

    lines   = strip(readlines(filename));
    isBlank = strlength(regexprep(lines, '[,\s]', '')) == 0;

    params = struct();
    hdrIdx = find(~isBlank, 1);
    gapIdx = hdrIdx + find(isBlank(hdrIdx+1:end), 1);   % first blank after the table
    if isempty(gapIdx)
        return
    end
    pLines = lines(gapIdx:end);
    pLines = pLines(~isBlank(gapIdx:end));

    names  = strip(extractBefore(pLines + ",", ","));
    values = strip(extractAfter(pLines, ","));
    values(ismissing(values)) = "";                      % name with no value

    names  = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));
    values = arrayfun(@rsvp_parse_value, values, 'UniformOutput', false);
    params = cell2struct(values, cellstr(names), 1);


    %recalculate screen deg/ppd
    %there is a bug in RSVP code that means we can't trust these values as
    %output from the task, and they are recoverable from the other
    %parameters anyway (that is how the task computes them)

    if all(ismember({'screen_res','screen_height','view_dist'},fieldnames(params)))
        %screen_width = params.screen_height * params.screen_res(1) / params.screen_res(2);
        params.screen_h_deg = (180 / pi) * 2 * atan(params.screen_height / (2 * params.view_dist));
        params.screen_w_deg = (180 / pi) * 2 * atan(params.screen_width / (2 * params.view_dist));
        params.ppd  = params.screen_res(2) / params.screen_h_deg;
    else
        warning('Screen parameters missing from summary paramters file. Check preprocessing carefully for correct values (e.g., ppd)')
    end
    
    params.task = 'rsvp';
    params.filename = filename;
end
