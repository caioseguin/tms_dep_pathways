function [idx, mask] = strsea(cell_string_vec, s, contains)

    if nargin == 2
        contains = 1;
    end

    if iscell(s)
        idx = [];
        mask = zeros(length(cell_string_vec), 1);
        for i = 1:length(s)
            [aux, bux] = strsea_fnc(cell_string_vec, s{i}, contains);
            idx = [idx; aux];
            mask(bux == 1) = 1;
        end
        idx = unique(idx);
    else
        [idx, mask] = strsea_fnc(cell_string_vec, s, contains);
    end

end

function [idx, mask] = strsea_fnc(cell_string_vec, s, contains)

    n = length(cell_string_vec);
    mask = zeros(n,1);

    for i = 1:n
        
        if contains
            if strfind(cell_string_vec{i}, s)
                mask(i) = 1;
            end
        else
            if strcmp(cell_string_vec{i}, s)
                mask(i) = 1;
            end
        end 
    end
    
    idx = find(mask);

end