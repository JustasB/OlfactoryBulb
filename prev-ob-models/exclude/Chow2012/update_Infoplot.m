function update_Infoplot(VST,Pcorr,H)
    set(H{1}, 'YData', uptri_1d(Pcorr));
    for i = 1:length(VST)
        set(H{i+1}, 'YData', VST{i});
    end
end