function val = uptri_1d(M)
    [r, c] = size(M);
    val = zeros(1, (c-1+c-r)*r/2);
    count = 1;
    for i = 1:r
        for j = i+1:c
            val(count) = M(i,j);
            count = count+1;
        end
    end
end % uptri_1d