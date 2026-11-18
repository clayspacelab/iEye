function [T, params] = rsvp_read_summary_csv(filename)
%RSVP_READ_SUMMARY_CSV Read a trial table and a block of name/value parameters from one CSV.
%
%   [T, params] = rsvp_read_summary_csv(filename)
%
%   Expected file layout:
%     1. A header row followed by one row per trial (any number of rows/columns).
%     2. One or more blank rows. Rows containing only commas/whitespace
%        (as written by Excel) also count as blank.
%     3. Any number of "name,value" rows. Values may contain unquoted commas,
%        e.g. "background_color,[0, 0, 0]", so everything after the FIRST
%        comma on a line is treated as the value.
%
%   Outputs
%     T      - MATLAB table. A column whose entries are all numeric becomes
%              double (empty cells -> NaN); all True/False becomes logical;
%              anything else becomes a string array.
%     params - struct with one field per parameter. Values are converted:
%                '1.3'                -> 1.3
%                'True' / 'False'     -> true / false
%                '[1920, 1080]'       -> [1920 1080]   (row vector)
%                '[[1, 2], [3, 4]]'   -> [1 2; 3 4]
%                '[a, b]'             -> {'a','b'}
%                'None' or empty      -> []
%                anything else        -> char, e.g. '192.168.1.5', 'pix'
%              Names that aren't valid MATLAB identifiers are made valid with
%              matlab.lang.makeValidName.
%
%   Example
%     [T, params] = readSummaryCSV('HL_run1_..._rsvp_sweep_summary.csv');
%     mean(T.Accuracy)
%     params.screen_res        % [1920 1080]

    if isstring(filename)
        filename = char(filename);
    end

    % ---- Read the whole file as text ---------------------------------------
    fid = fopen(filename, 'r', 'n', 'UTF-8');
    if fid < 0
        error('readSummaryCSV:cannotOpen', 'Could not open file: %s', filename);
    end
    txt = fread(fid, [1 Inf], '*char');
    fclose(fid);

    if ~isempty(txt) && double(txt(1)) == 65279      % strip UTF-8 byte-order mark
        txt = txt(2:end);
    end

    lines   = strtrim(regexp(txt, '\r\n|\n|\r', 'split'));
    isBlank = cellfun(@(s) isempty(regexprep(s, '[,\s]', '')), lines);

    % ---- Locate the table: first non-blank line to the next blank line -----
    hdrIdx = find(~isBlank, 1, 'first');
    if isempty(hdrIdx)
        error('readSummaryCSV:emptyFile', 'File contains no data: %s', filename);
    end
    nextBlank = find(isBlank(hdrIdx+1:end), 1, 'first');
    if isempty(nextBlank)
        tblEnd = numel(lines);                      % no parameter section
    else
        tblEnd = hdrIdx + nextBlank - 1;
    end

    T = parseTable(lines{hdrIdx}, lines(hdrIdx+1:tblEnd));

    % ---- Everything non-blank after the table is a parameter ---------------
    rest   = lines(tblEnd+1:end);
    params = parseParams(rest(~isBlank(tblEnd+1:end)));
end


% =========================================================================
function T = parseTable(hdrLine, dataLines)
    names = splitCsvLine(hdrLine);
    nCols = numel(names);
    nRows = numel(dataLines);

    raw = repmat({''}, nRows, nCols);
    for r = 1:nRows
        f = splitCsvLine(dataLines{r});
        if numel(f) > nCols
            if any(~cellfun(@isempty, f(nCols+1:end)))
                warning('readSummaryCSV:extraFields', ...
                    'Table row %d has %d fields but the header has %d; extra fields ignored.', ...
                    r, numel(f), nCols);
            end
            f = f(1:nCols);
        end
        raw(r, 1:numel(f)) = f;
    end

    % Drop columns with no header and no data (trailing commas from Excel).
    dropCol = cellfun(@isempty, names) & all(cellfun(@isempty, raw), 1);
    names(dropCol)  = [];
    raw(:, dropCol) = [];
    nCols = numel(names);

    for j = 1:nCols
        if isempty(names{j})
            names{j} = sprintf('Var%d', j);
        end
    end
    varNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));

    cols = cell(1, nCols);
    for j = 1:nCols
        cols{j} = convertColumn(raw(:, j));
    end
    T = table(cols{:}, 'VariableNames', varNames);

    if ~isequal(varNames, names)                    % keep original headers
        T.Properties.VariableDescriptions = names;
    end
end


% =========================================================================
function col = convertColumn(c)
    c = strtrim(c(:));
    if isempty(c)
        col = zeros(0, 1);
        return
    end
    isEmpty = cellfun(@isempty, c);
    num     = str2double(c);
    numOK   = (~isnan(num) | isEmpty | strcmpi(c, 'nan'));

    if all(isEmpty)
        col = NaN(numel(c), 1);
    elseif all(numOK) && isreal(num)
        col = num;
    elseif ~any(isEmpty) && all(ismember(lower(c), {'true', 'false'}))
        col = strcmpi(c, 'true');
    else
        col = string(c);
    end
end


% =========================================================================
function params = parseParams(lines)
    params = struct();
    for k = 1:numel(lines)
        s = regexprep(lines{k}, ',[,\s]*$', '');    % trailing commas from Excel

        % Split at the first comma that is not inside double quotes.
        c = regexp(s, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'once');
        if isempty(c)
            key = s;
            val = '';
        else
            key = s(1:c-1);
            val = s(c+1:end);
        end
        key = unquote(strtrim(key));
        val = unquote(strtrim(val));
        if isempty(key)
            continue
        end

        fname = matlab.lang.makeValidName(key);
        if isfield(params, fname)
            warning('readSummaryCSV:duplicateParam', ...
                'Duplicate parameter "%s"; keeping the last value.', key);
        end
        params.(fname) = parseValue(val);
    end
end


% =========================================================================
function v = parseValue(s)
    s = strtrim(s);
    if isempty(s) || any(strcmpi(s, {'none', 'null'}))
        v = [];
    elseif strcmpi(s, 'true')
        v = true;
    elseif strcmpi(s, 'false')
        v = false;
    elseif strcmpi(s, 'nan')
        v = NaN;
    elseif numel(s) >= 2 && ((s(1) == '[' && s(end) == ']') || ...
                             (s(1) == '(' && s(end) == ')'))
        v = parseList(s);
    else
        d = str2double(s);
        if ~isnan(d) && isreal(d)
            v = d;
        else
            v = s;                                  % plain text, e.g. an IP address
        end
    end
end


% =========================================================================
function v = parseList(s)
    inner = strtrim(s(2:end-1));
    if isempty(inner)
        v = [];
        return
    end

    % Nested list, e.g. [[1, 2], [3, 4]] -> one matrix row per sub-list.
    if any(inner == '[') || any(inner == '(')
        tok  = regexp(inner, '[\[\(]([^\[\]\(\)]*)[\]\)]', 'tokens');
        rows = cellfun(@(t) parseList(['[' t{1} ']']), tok, 'UniformOutput', false);
        lens = cellfun(@numel, rows);
        isNum = cellfun(@(r) isnumeric(r) || islogical(r), rows);
        if ~isempty(rows) && all(isNum) && all(lens == lens(1))
            v = vertcat(rows{:});
        else
            v = rows;                               % ragged or mixed: cell of rows
        end
        return
    end

    items = strtrim(strsplit(inner, ','));
    items(cellfun(@isempty, items)) = [];           % e.g. Python 1-tuple "(5,)"
    if isempty(items)
        v = [];
        return
    end

    nums = str2double(items);
    if all(~isnan(nums) | strcmpi(items, 'nan')) && isreal(nums)
        v = nums;
    elseif all(ismember(lower(items), {'true', 'false'}))
        v = strcmpi(items, 'true');
    else
        v = cellfun(@unquote, items, 'UniformOutput', false);
    end
end


% =========================================================================
function f = splitCsvLine(s)
    % Split on commas that are not inside double quotes.
    f = strtrim(regexp(s, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split'));
    f = cellfun(@unquote, f, 'UniformOutput', false);
end


function s = unquote(s)
    % Remove one layer of matching surrounding quotes ("..." or '...').
    if numel(s) >= 2 && any(s(1) == '"''') && s(end) == s(1)
        q = s(1);
        s = strrep(s(2:end-1), [q q], q);
    end
end
