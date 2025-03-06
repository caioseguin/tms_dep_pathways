function x = block_wei_avg(X, B1, B2, W1)

    if nargin == 3
        W1 = ones(size(B1));
    end

    b1 = B1(:);
    b2 = B2(:);
    w1 = W1(:);

    b1 = b1(b1>0);
    b2 = b2(b2>0);
    w1 = w1(b1>0);

    Xb = X(b1,b2);
    Wb = repmat(w1, 1, length(b2));

    xb = Xb(:);
    wb = Wb(:);
    idx = xb ~= Inf & wb ~= Inf;
    x = sum(xb(idx).*wb(idx))./sum(wb(idx));

end