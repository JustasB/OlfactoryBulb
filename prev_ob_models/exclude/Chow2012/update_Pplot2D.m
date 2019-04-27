function update_Pplot2D(P,coord,H,Haxis,Pmin,Pmax)
    [Nc, Ns] = size(P);
    Cmin = min(coord')'; Cmax = max(coord')';
    c = ceil(sqrt(Ns)); r = ceil(Ns/c);
    for ii = 1:r
        for jj = 1:c
            kk = (ii-1)*c+jj;
            if kk <= Ns
                % subplot(r,c,kk);
                subplot(1,r*c,kk);
                temp = NaN*zeros(ceil(Cmax-Cmin)'+1);
                for ll = 1:Nc
                    temp(ceil(coord(1,ll)-Cmin(1))+1, ceil(coord(2,ll)-Cmin(2))+1) = P(ll,kk);
                end
                set(H{kk}, 'CData', temp(:,end:-1:1)');
                if nargin < 6
                    mm = mean_excluNaN(temp(:,end:-1:1));
                    ss = std_excluNaN(temp(:,end:-1:1));
                    Pmin = mm-3*ss;
                    Pmax = mm+3*ss;
                end
                set(Haxis{kk}, 'CLim', [Pmin Pmax]);
            end
        end
    end
end % setup_Pplot2D