function update_Pplot(P,Pcorr,Ccorr,rg,H,Wmg,Wgm,N_t)

    set(H{1}, 'CData', P);
    set(H{2}, 'CData', Ccorr);
    set(H{3}, 'CData', Pcorr);
    if nargin > 5
        conn = Wmg*Wgm;
        set(H{4}, 'CData', conn-diag(NaN*diag(conn)));
    end
end % setup_Pplot