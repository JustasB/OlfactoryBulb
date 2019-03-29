function [H P2daxis] = setup_Pplot2D(P,coord,N,name,Pmin,Pmax)
    [Nc Ns] = size(P);
    H = cell(1,Ns);
    P2daxis = cell(1,Ns);
    figure(N); clf; set(gcf,'DoubleBuffer','on');
%     if nargin < 6
%         Pmin = min(min(P)); Pmax = max(max(P));
%     end
    Cmin = min(coord')'; Cmax = max(coord')';
    c = ceil(sqrt(Ns)); r = ceil(Ns/c);
    for ii = 1:r
        for jj = 1:c
            kk = (ii-1)*c+jj;
            if kk <= Ns
                subplot(r,c,kk);
                % subplot(1,r*c,kk);
                temp = NaN*zeros(ceil(Cmax-Cmin)'+1);
                for ll = 1:Nc
                    temp(ceil(coord(1,ll)-Cmin(1))+1, ceil(coord(2,ll)-Cmin(2))+1) = P(ll,kk);
                end
                H{kk} = imagesc(temp(:,end:-1:1)');
                P2daxis{kk} = gca;
                colormap(jet);
                if nargin < 6
                    mm = mean_excluNaN(temp(:,end:-1:1));
                    ss = std_excluNaN(temp(:,end:-1:1));
                    Pmin = mm-3*ss;
                    Pmax = mm+3*ss;
                end
                caxis([Pmin Pmax]);      
                set(gca, 'color', 'black');
                axis off; axis image;
                if nargin > 3
                    title(name{kk});
                end
            end
        end
    end
end % setup_Pplot2D