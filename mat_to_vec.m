function i = mat_to_vec(I)
    if length(I) == 3
        i = [I(1,1), I(1,2), I(1,3), I(2,2), I(2,3), I(3,3)];
    elseif length(I) == 2
        i = [I(1,1), I(1,2), I(2,2)];
    end
end

