function [f] = flat(X)
    f = reshape(X, numel(X), 1);
end