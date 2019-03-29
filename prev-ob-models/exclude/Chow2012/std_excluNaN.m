function val = std_excluNaN(V)
    temp = V(find(isnan(V)==0));
    if length(temp)>0
        val = std(temp);
    else val = NaN;
    end
end