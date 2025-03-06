function T_new = add_to_path_table(T, path_string, n)

    idx = strcmp(T.path, path_string);

    if any(idx)
        T.count(idx) = T.count(idx) + n;
        T_new = T;
    else
        T_new = [T; {path_string, n}];
    end

end