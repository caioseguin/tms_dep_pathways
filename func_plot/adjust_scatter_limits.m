function adjust_scatter_limits(axHandle, paddingFactor)
% adjustAxisLimits Adjusts axis limits to include padding if necessary
%
%   adjustAxisLimits(axHandle, paddingFactor) adjusts the x- and y-axis 
%   limits of the specified axes handle to include extra space around the data.
%   - axHandle: Handle to the axes to adjust (default: gca if not provided).
%   - paddingFactor: Fraction of the data range to use as padding 
%                    (default: 0.1, i.e., 10% of the range).

    % Default values if arguments are not provided
    if nargin < 1 || isempty(axHandle)
        axHandle = gca; % Use current axes if not specified
    end
    if nargin < 2 || isempty(paddingFactor)
        paddingFactor = 0.1; % Default padding of 10%
    end

    % Get current axis limits
    xLimits = get(axHandle, 'XLim');
    yLimits = get(axHandle, 'YLim');

    % Get the data range
    xData = get(axHandle.Children, 'XData');
    yData = get(axHandle.Children, 'YData');
    
    % Handle cases where there are multiple children in the axes
    if iscell(xData)
        xData = cell2mat(xData(:));
    end
    if iscell(yData)
        yData = cell2mat(yData(:));
    end

    % Determine the min and max of the data
    xMin = min(xData(:));
    xMax = max(xData(:));
    yMin = min(yData(:));
    yMax = max(yData(:));

    % Calculate required padding
    xPadding = paddingFactor * (xMax - xMin);
    yPadding = paddingFactor * (yMax - yMin);

    % Determine new limits only if padding is insufficient
    newXLim = xLimits;
    newYLim = yLimits;

    % Check if the current limits provide enough padding
    if xLimits(1) > xMin - xPadding
        newXLim(1) = xMin - xPadding;
    end
    if xLimits(2) < xMax + xPadding
        newXLim(2) = xMax + xPadding;
    end
    if yLimits(1) > yMin - yPadding
        newYLim(1) = yMin - yPadding;
    end
    if yLimits(2) < yMax + yPadding
        newYLim(2) = yMax + yPadding;
    end

    % Apply the new limits if they've changed
    if ~isequal(newXLim, xLimits)
        set(axHandle, 'XLim', newXLim);
    end
    if ~isequal(newYLim, yLimits)
        set(axHandle, 'YLim', newYLim);
    end
end
