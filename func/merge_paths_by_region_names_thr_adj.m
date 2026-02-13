function T_out = merge_paths_by_region_names_thr_adj( ...
    T_in, label_threshold, keep_at_least_one, max_labels, ...
    atlasAdj, parc_labels, cover_all_rows)
% Merge by region NAMES (ignore labels for grouping), but:
%  - filter labels by contribution threshold,
%  - enforce adjacency between merged labels (atlasAdj + parc_labels),
%  - and recompute count/percentage from rows that MATCH the kept sets.
% Optionally partition each name-only group to cover all rows.
%
% Inputs
%   T_in : table with vars
%           - path (string/char like '_rh_x_1__rh_y_2_')
%           - count (numeric)
%           - percentage (optional). If missing, computed as count/sum(count).
%   label_threshold : fraction or percent (>1 becomes /100). Default 0.
%   keep_at_least_one : logical. Default true.
%   max_labels : positive integer cap per node. Default 3.
%   atlasAdj : N×N logical/sparse adjacency (parcel IDs 1..N).
%   parc_labels : N×1 cell array; parc_labels{id} is token for that ID.
%   cover_all_rows : logical. If true (default), iteratively partition
%                    each group so all rows are assigned to some merged row.
%
% Output
%   T_out : table with {'path','count','percentage'}, sorted by count desc.

    if nargin < 2 || isempty(label_threshold),  label_threshold = 0; end
    if nargin < 3 || isempty(keep_at_least_one), keep_at_least_one = true; end
    if nargin < 4 || isempty(max_labels),        max_labels = 3; end
    if nargin < 7 || isempty(cover_all_rows),    cover_all_rows = true; end

    if label_threshold > 1, label_threshold = label_threshold/100; end
    if ~isstring(T_in.path), T_in.path = string(T_in.path); end

    % Ensure percentage exists
    if ~ismember('percentage', T_in.Properties.VariableNames)
        total_cnt = sum(T_in.count, 'omitnan');
        T_in.percentage = T_in.count / max(total_cnt, eps);
    end

    % Build token -> atlas IDs (case-insensitive)
    N = numel(parc_labels);
    if size(atlasAdj,1) ~= N || size(atlasAdj,2) ~= N
        error('atlasAdj must be N×N with N == numel(parc_labels).');
    end
    atlasAdj = logical(atlasAdj);
    tokenToIds = containers.Map('KeyType','char','ValueType','any');
    for id = 1:N
        tok = lower(string(parc_labels{id}));
        if isKey(tokenToIds, tok)
            tokenToIds(tok) = [tokenToIds(tok), id];
        else
            tokenToIds(tok) = id;
        end
    end

    % Parse all paths
    n = height(T_in);
    nodeNamesAll  = cell(n,1);
    nodeLabelsAll = cell(n,1);
    for i = 1:n
        nodes = split_path_to_nodes(T_in.path(i));
        L = numel(nodes);
        names_i  = strings(1,L);
        labels_i = nan(1,L);
        for k = 1:L
            [nm, lab] = split_name_label(nodes{k});
            names_i(k)  = string(nm);
            labels_i(k) = lab;
        end
        nodeNamesAll{i}  = names_i;
        nodeLabelsAll{i} = labels_i;
    end

    % Group by NAME-ONLY trajectory
    keys = strings(n,1);
    for i = 1:n
        keys(i) = "_" + strjoin(nodeNamesAll{i}, "__") + "_";
    end
    [G, keyLevels] = findgroups(keys);

    mergedPaths = strings(0,1);
    mergedCount = zeros(0,1);
    mergedPct   = zeros(0,1);

    % Process each group (optionally partition to cover all rows)
    groups = unique(G(:))';
    for g = groups
        idx_all = find(G==g);
        names_g = nodeNamesAll{idx_all(1)};
        L = numel(names_g);

        % Stack labels
        R = numel(idx_all);
        labels_stack = nan(R,L);
        for r = 1:R
            labels_stack(r,:) = nodeLabelsAll{idx_all(r)};
        end

        % Work list of unassigned rows (indices into idx_all)
        unassigned = true(R,1);

        while true
            inds = find(unassigned)';
            if isempty(inds), break; end

            % Contributions on the CURRENT remainder only
            pct_group = sum(T_in.percentage(idx_all(inds)), 'omitnan');

            % If negligible remainder, stop
            if pct_group <= eps, break; end

            % For each node position, compute label shares on remainder
            kept_sets = cell(1,L);
            contrib_map = cell(1,L); % store (uniq_labels, contribs) for fallback

            for k = 1:L
                nm = names_g(k);
                labs_k = labels_stack(inds, k);
                valid = ~isnan(labs_k);
                uniq_labs = unique(labs_k(valid))';
                if isempty(uniq_labs)
                    kept_sets{k} = []; % unlabeled node (e.g., thalamusproper)
                    contrib_map{k} = {[],[]};
                    continue;
                end
                contrib = zeros(size(uniq_labs));
                for j = 1:numel(uniq_labs)
                    lab = uniq_labs(j);
                    rows_with_lab = (labs_k == lab);
                    contrib(j) = sum(T_in.percentage(idx_all(inds(rows_with_lab))), 'omitnan');
                end
                shares = contrib ./ max(pct_group, eps);
                keep_idx = shares >= label_threshold;
                if ~any(keep_idx) && keep_at_least_one
                    [~, imax] = max(contrib);
                    keep_idx(imax) = true;
                end
                % adjacency-constrained greedy selection
                cand_labels   = uniq_labs(keep_idx);
                cand_contribs = contrib(keep_idx);
                [~, ord_desc] = sort(cand_contribs,'descend');
                cand_labels = cand_labels(ord_desc);

                kept = [];
                kept_ids = [];
                for j = 1:numel(cand_labels)
                    lab = cand_labels(j);
                    token = lower(char(nm + "_" + string(lab)));
                    if ~isKey(tokenToIds, token)
                        % If token missing, skip this label
                        continue;
                    end
                    ids = tokenToIds(token);
                    if isempty(kept)
                        kept(end+1) = lab; %#ok<AGROW>
                        kept_ids = [kept_ids, ids]; %#ok<AGROW>
                    else
                        if any_is_adjacent(ids, kept_ids, atlasAdj)
                            kept(end+1) = lab; %#ok<AGROW>
                            kept_ids = [kept_ids, ids]; %#ok<AGROW>
                        end
                    end
                    if numel(kept) >= max_labels
                        break;
                    end
                end
                % Fallback: if adjacency blocked everything, keep top-1 by contrib
                if isempty(kept) && ~isempty(cand_labels)
                    kept = cand_labels(1);
                end
                kept_sets{k} = sort(unique(kept));
                contrib_map{k} = {uniq_labs, contrib};
            end

            % Identify rows (within remainder) that MATCH kept sets
            match_mask = true(numel(inds),1);
            for k = 1:L
                labs_k = labels_stack(inds, k);
                S = kept_sets{k};
                if isempty(S) || all(isnan(labs_k))
                    % unlabeled node (or all NaN): no restriction
                    continue;
                end
                present = false(numel(inds),1);
                for s = 1:numel(S)
                    present = present | (labs_k == S(s));
                end
                match_mask = match_mask & present;
            end

            matched_inds = inds(match_mask);

            % Robust fallback: if no actual row matches the multi-node kept sets,
            % choose a single top label per labeled node and match exactly that combo.
            if isempty(matched_inds)
                mode_labels = nan(1,L);
                for k = 1:L
                    nm = names_g(k);
                    labs_k = labels_stack(inds, k);
                    if all(isnan(labs_k))
                        mode_labels(k) = NaN;
                        continue;
                    end
                    % pick highest-contrib single label on remainder
                    uv  = contrib_map{k}{1};
                    con = contrib_map{k}{2};
                    if isempty(uv)
                        mode_labels(k) = NaN;
                    else
                        [~, imax] = max(con);
                        mode_labels(k) = uv(imax);
                    end
                end
                % exact match for that label combo
                exact_mask = true(numel(inds),1);
                for k = 1:L
                    labs_k = labels_stack(inds, k);
                    if ~isnan(mode_labels(k))
                        exact_mask = exact_mask & (labs_k == mode_labels(k));
                    end
                end
                matched_inds = inds(exact_mask);
                % and overwrite kept_sets with the singletons we picked
                for k = 1:L
                    if ~isnan(mode_labels(k))
                        kept_sets{k} = mode_labels(k);
                    else
                        kept_sets{k} = [];
                    end
                end
            end

            % Build path string from kept_sets
            mergedNodes = strings(1,L);
            for k = 1:L
                nm = names_g(k);
                S  = kept_sets{k};
                if isempty(S)
                    mergedNodes(k) = nm;
                elseif numel(S)==1
                    mergedNodes(k) = nm + "_" + string(S);
                else
                    mergedNodes(k) = nm + "_(" + strjoin(string(S), ",") + ")";
                end
            end
            path_str = "_" + strjoin(mergedNodes, "__") + "_";

            % Sum counts/pcts over matched rows ONLY
            cnt = sum(T_in.count(idx_all(matched_inds)), 'omitnan');
            pct = sum(T_in.percentage(idx_all(matched_inds)), 'omitnan');

            mergedPaths(end+1,1) = path_str; %#ok<AGROW>
            mergedCount(end+1,1) = cnt;      %#ok<AGROW>
            mergedPct(end+1,1)   = pct;      %#ok<AGROW>

            % Mark matched rows as assigned
            if cover_all_rows
                unassigned(matched_inds) = false;
            else
                % If not covering all rows, we do a single pass per name-only group.
                break;
            end
        end
    end

    % If multiple name-only groups produced identical path strings,
    % coalesce them here (sum counts/pcts).
    [uPaths, ~, ic] = unique(mergedPaths);
    aggCount = accumarray(ic, mergedCount);
    aggPct   = accumarray(ic, mergedPct);
    T_out = table(uPaths, aggCount, aggPct, ...
        'VariableNames', {'path','count','percentage'});

    % Sort by descending count
    T_out = sortrows(T_out, 'count', 'descend');
end

% ---------- helpers ----------
function nodes = split_path_to_nodes(pathStr)
    s = strip(string(pathStr), '_');
    if strlength(s)==0, nodes = {}; return; end
    parts = split(s, "__");
    nodes = cellstr(parts);
end

function [nameOut, labelOut] = split_name_label(nodeStr)
    nameOut  = nodeStr;
    labelOut = NaN;
    us = find(nodeStr=='_',1,'last');
    if isempty(us), return; end
    if us < strlength(nodeStr)
        tail = extractAfter(nodeStr, us);
        num = str2double(tail);
        if ~isnan(num) && isfinite(num) && all(isstrprop(tail,'digit'))
            nameOut  = extractBefore(nodeStr, us);
            labelOut = num;
        end
    end
end

function tf = any_is_adjacent(idsA, idsB, A)
    if isempty(idsB), tf = false; return; end
    idsA = unique(idsA(idsA>=1 & idsA<=size(A,1)));
    idsB = unique(idsB(idsB>=1 & idsB<=size(A,1)));
    if isempty(idsA) || isempty(idsB), tf = false; return; end
    tf = any(any(A(idsA, idsB)));
end
