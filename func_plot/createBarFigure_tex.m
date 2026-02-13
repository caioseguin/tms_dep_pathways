function createBarFigure_tex(leftValues, rightValues, leftLabels, rightLabels, tickPositions, font_size, label_font_size, font_name)
% createBarFigure_tex Creates a horizontal bar figure with left and right bars.
% Uses MATLAB's TeX interpreter to render mixed font sizes (single text object).
%
% Inputs:
%   - leftValues, rightValues: vectors (same length)
%   - leftLabels, rightLabels: cellstr arrays (same length as values)
%   - tickPositions (optional): vector of positive tick positions (mirrored)
%   - font_size: base font size for region names
%   - label_font_size: font size for label parts inside parentheses "(...)"
%   - font_name (optional): font family for text (default: 'Helvetica')
%
% Example:
%   createBarFigure_tex(L, R, Llbl, Rlbl, [], 12, 8);                 % Helvetica (default)
%   createBarFigure_tex(L, R, Llbl, Rlbl, [], 12, 8, 'Arial');        % Use Arial

    if nargin < 8 || isempty(font_name)
        font_name = 'Helvetica';
    end

    % ---- Prettify labels (your original replacements) ----
    replacements = {
        'accumbensarea','accumbens';
        'superiorfrontal','sup frontal';
        'caudalanteriorcingulate','cau ant cingulate';
        'thalamusproper','thalamus';
        'temporalpole','temp pole';
        'parstriangularis','pars tri';
        'rostralmiddlefrontal','ros mid fontral';
        'RH','R';
        'LH','L'
    };
    for r = 1:size(replacements,1)
        leftLabels  = replace(leftLabels,  replacements{r,1}, replacements{r,2});
        rightLabels = replace(rightLabels, replacements{r,1}, replacements{r,2});
    end

    % ---- Validate ----
    if length(leftValues) ~= length(rightValues) || ...
       length(leftValues) ~= length(leftLabels)  || ...
       length(rightValues) ~= length(rightLabels)
        error('All inputs must have the same length.');
    end

    % ---- Plot bars ----
    bar_offset = 0.5;
    y_right = 1:length(leftValues);
    y_left  = (1:length(leftValues)) + bar_offset;

    barh(y_right, rightValues, 'FaceColor', [0.3010 0.7450 0.9330], ...
        'EdgeColor', 'k', 'LineWidth', 1, 'BarWidth', 0.4);
    hold on;
    barh(y_left, -leftValues, 'FaceColor', [0.4940 0.1840 0.5560], ...
        'EdgeColor', 'k', 'LineWidth', 1, 'BarWidth', 0.4);

    % ---- Axes ----
    maxValue = max(max(rightValues), max(leftValues));
    xlim([-maxValue*1.12, maxValue*1.12]);
    ylim([0.25, length(y_right) + 2.5*bar_offset]);

    if nargin <= 5 || isempty(tickPositions)
        xticks(linspace(-maxValue, maxValue, 7));
        xticklabels(round(abs(linspace(-maxValue, maxValue, 7))));
    else
        xticks(sort([-tickPositions, 0, tickPositions]));
        xticklabels(abs(sort([-tickPositions, 0, tickPositions])));
    end

    ax = gca;
    ax.XAxisLocation = 'origin';
    xlabel('Contribution (%)');
    box off;
    yticks([]);
    ax.YAxis.Visible = 'off';
    ax.LineWidth = 1;
    line([0,0], ylim(), 'Color', 'black', 'LineWidth', 2);

    % ---- Labels (single text object each; TeX interpreter) ----
    spacing = maxValue*0.025;

    for i = 1:length(y_right)
        txt = toTexMixedSize(rightLabels{i}, font_size, label_font_size, font_name);
        text(rightValues(i) + spacing, y_right(i), txt, ...
            'HorizontalAlignment','left',  'Interpreter','tex');
    end

    for i = 1:length(y_left)
        txt = toTexMixedSize(leftLabels{i}, font_size, label_font_size, font_name);
        text(-leftValues(i) - spacing, y_left(i), txt, ...
            'HorizontalAlignment','right', 'Interpreter','tex');
    end
end

% ===== Helper: build one TeX string with mixed sizes (region text vs "(...)") =====
function out = toTexMixedSize(str_in, baseSize, parenSize, fontName)
    s = escapeTex(str_in);

    % Split into non-parenthesis text vs '(...)' groups (keep delimiters)
    tokens = regexp(s, '(\([^)]+\))', 'split');
    parens = regexp(s, '(\([^)]+\))', 'match');

    chunks = strings(0);
    for k = 1:numel(tokens)
        if ~isempty(tokens{k})
            % Region text: base size
            chunks(end+1) = "\fontname{" + fontName + "}\fontsize{" + baseSize + "}" + tokens{k}; %#ok<AGROW>
        end
        if k <= numel(parens) && ~isempty(parens{k})
            % Parentheses text: smaller size
            chunks(end+1) = "\fontname{" + fontName + "}\fontsize{" + parenSize + "}" + parens{k}; %#ok<AGROW>
        end
    end
    out = strjoin(chunks, '');
end

% ===== Helper: escape MATLAB TeX special characters =====
function s = escapeTex(s)
    % Backslash first
    s = strrep(s, '\', '\\');
    % Braces and others
    s = strrep(s, '{', '\{');
    s = strrep(s, '}', '\}');
    s = strrep(s, '_', '\_');
    s = strrep(s, '^', '\^{}');
    s = strrep(s, '#', '\#');
    s = strrep(s, '%', '\%');
    s = strrep(s, '$', '\$');
    s = strrep(s, '&', '\&');
    s = strrep(s, '~', '\~{}');
end
