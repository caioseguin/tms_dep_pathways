function adjust_histogram_limits(axHandle, paddingFactor)
% adjustHistogramLimits Adjusts axis limits for histograms
%
%   adjustHistogramLimits(axHandle, paddingFactor) adjusts the x-axis limits 
%   of the specified axes handle to include extra space around the data.
%   The y-axis starts from zero without any padding.
%   - axHandle: Handle to the axes containing the histogram (default: gca).
%   - paddingFactor: Fraction of the data range to use as padding 
%                    (default: 0.1, i.e., 10% of the range).

    % Default values if arguments are not provided
    if nargin < 1 || isempty(axHandle)
        axHandle = gca; % Use current axes if not specified
    end
    if nargin < 2 || isempty(paddingFactor)
        paddingFactor = 0.1; % Default padding of 10%
    end

    % Get the histogram object
    histObj = findobj(axHandle, 'Type', 'Histogram');
    if isempty(histObj)
        error('No histogram found in the specified axes.');
    end

    % Get histogram data
    xData = get(histObj, 'BinEdges'); % Bin edges for x-axis
    yData = get(histObj, 'Values'); % Heights for y-axis

    % Calculate the x-axis padding
    xMin = min(xData);
    xMax = max(xData);
    xPadding = paddingFactor * (xMax - xMin);

    % Adjust x-axis limits
    newXLim = [xMin - xPadding, xMax + xPadding];
    set(axHandle, 'XLim', newXLim);

    % Adjust y-axis limits to start from zero
    yMax = max(yData);
    newYLim = [0, yMax + paddingFactor * yMax]; % Optional upper padding
    set(axHandle, 'YLim', newYLim);
end
