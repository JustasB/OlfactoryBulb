function update_Iplot(cont_density,t_axis,I,Icorr,Iage,Wmg,N_popu,H)
    if length(find(Iage>=0)) <= 0
        Iage(1) = 0;
    end
    set(H{1}, 'CData', I);
    set(H{2}, 'CData', Icorr);
end % setup_Iplot