function T = rsvp_read_summary_table(filename)
%RSVP_READ_SUMMARY_TABLE Read the trial table at the top of a summary CSV.
%
%   T = rsvp_read_summary_table(filename)
%
%   The header is line 1 and the table runs to the first blank line. Lines
%   of only commas/whitespace count as blank. Any number of trials, blank
%   lines and parameters is fine.
%
%   readtable reads the whole file, then the parameter rows are trimmed off.
%   Because readtable's type detection also sees the parameter rows (e.g.
%   parameter names land in column 1), every column is read as text and the
%   columns whose table entries are all numbers are converted to double
%   afterwards. Other columns stay string.
%
%   See also rsvp_read_params, readtable.

    opts = detectImportOptions(filename, 'Delimiter', ',', 'VariableNamesLine', 1, ...
                               'NumHeaderLines', 0, 'EmptyLineRule', 'read');
    opts = setvartype(opts, 'string');
    T = readtable(filename, opts);

    % Keep rows up to the first blank line, and drop the extra columns that
    % long parameter values (e.g. trial_onset lists) add past the header.
    isBlank = all(ismissing(T), 2);
    T = T(1:find([isBlank; true], 1) - 1, :);
    extra = matches(T.Properties.VariableNames, "Var" + digitsPattern) & all(ismissing(T), 1);
    T(:, extra) = [];

    isNum = varfun(@isNumericText, T, 'OutputFormat', 'uniform');
    T = convertvars(T, isNum, 'double');
end


function tf = isNumericText(x)
    % True if every entry is a number, empty, or "nan" (as Python writes it).
    tf = all(~isnan(str2double(x)) | ismissing(x) | strcmpi(x, "nan"));
end
