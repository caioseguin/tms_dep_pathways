function L = log_transform(C, offset)
    
    if nargin == 1
        offset = max(1e-14, min(C(C>0))); % 1e-14 is the minimum offset so that log10(C./(max(C(:)) + offset)) ~= 0
    end

    assert(min(C(:)) >= 0, 'C must not contain negative values');

    % Check if C has values larger than 1
    if ~isempty(find(C >= 1, 1))
        C = C./(max(C(:)) + offset); 
    end

    L = -log10(C);
    L(L == Inf) = 0;

end