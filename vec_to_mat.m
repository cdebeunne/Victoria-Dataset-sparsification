function I = vec_to_mat(i)
    i = i(~isnan(i));
    if (length(i) == 6)
        I = [i(1), i(2), i(3);
            i(2), i(4), i(5);
            i(5), i(3), i(6)];
    else
        I = [i(1) i(2);
             i(2) i(3)];
    end
            
end