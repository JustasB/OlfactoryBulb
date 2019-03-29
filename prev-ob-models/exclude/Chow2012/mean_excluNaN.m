function val = mean_excluNaN(V)
    temp = V(find(isnan(V)==0));
    if length(temp)>0
        val = mean(temp);
    else val = NaN;
    end
end