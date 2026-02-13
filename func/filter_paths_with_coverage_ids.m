function T_out = filter_paths_with_coverage_ids(T_in, method, param)
%FILTER_PATHS_WITH_COVERAGE_IDS  Filter numeric-encoded paths with start/end coverage.
%
% T_out = filter_paths_with_coverage_ids(T_in, method, param)
%
% Inputs:
%   T_in   : table with variables:
%              - path  : cell array; each cell is a numeric vector of node IDs (start..end)
%              - count : numeric vector
%   method : 'cumshare' -> keep fewest rows whose cumulative count share >= param
%            'topfrac'  -> keep top fraction of rows by count
%   param  : fraction in (0,1]; if >1 && <=100, treated as percent (e.g., 50 -> 0.5)
%
% Coverage rule: Ensure every unique START ID and every unique END ID present in T_in
% appears in T_out at least once (not necessarily every start–end pair).
%
% Output:
%   T_out  : filtered table, sorted by count desc

    % --- basic checks & normalize param ---
    T = T_in;
    if ~iscell(T.path), error('T_in.path must be a cell array of numeric vectors.'); end
    if height(T) ~= numel(T.count), error('T_in.count must have one value per row.'); end
    n = height(T);
    counts = T.count(:);

    if param > 1 && param <= 100, param = param/100; end
    if param <= 0 || param > 1, error('param must be in (0,1] (or percent in (0,100]).'); end

    % --- extract start/end IDs ---
    try
        startID = cellfun(@(v) v(1), T.path);
        endID   = cellfun(@(v) v(end), T.path);
    catch
        error('Each T_in.path cell must be a non-empty numeric vector.');
    end

    % --- initial selection by method ---
    keep = false(n,1);

    switch lower(method)
        case 'cumshare'
            total = sum(counts);
            [c_sorted, ord] = sort(counts, 'descend');
            cum = cumsum(c_sorted);
            k = find(cum >= param * total, 1, 'first');
            if isempty(k), k = n; end
            keep(ord(1:k)) = true;

        case 'topfrac'
            [~, ord] = sort(counts, 'descend');
            k = max(1, ceil(param * n));
            keep(ord(1:k)) = true;

        otherwise
            error('Unknown method "%s". Use "cumshare" or "topfrac".', method);
    end

    % --- enforce coverage of ALL starts present in T_in ---
    allStarts = unique(startID);
    coveredStarts = unique(startID(keep));
    missingStarts = setdiff(allStarts, coveredStarts);
    for s = 1:numel(missingStarts)
        idx = find(startID == missingStarts(s));
        % choose highest-count row for this start
        [~, jrel] = max(counts(idx));
        keep(idx(jrel)) = true;
    end

    % --- enforce coverage of ALL ends present in T_in ---
    allEnds = unique(endID);
    coveredEnds = unique(endID(keep));
    missingEnds = setdiff(allEnds, coveredEnds);
    for e = 1:numel(missingEnds)
        idx = find(endID == missingEnds(e));
        [~, jrel] = max(counts(idx));
        keep(idx(jrel)) = true;
    end

    % --- output sorted by count desc ---
    T_out = T(keep, :);
    T_out = sortrows(T_out, 'count', 'descend');
end
