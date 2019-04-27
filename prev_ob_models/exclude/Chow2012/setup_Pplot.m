function H = setup_Pplot(P,Pcorr,Ccorr,rg,N,Wmg,Wgm,time_axis,N_t)

    H = cell(1,4);
    figure(N); clf; set(gcf,'DoubleBuffer','on');
    subplot(2,2,1);
    H{1} = imagesc(P); colormap(jet);
    xlabel('stimulus'); ylabel('MC'); title('P');  axis square;
    subplot(2,2,2);
    H{2} = imagesc(Ccorr); colormap(jet);
    title('channel corr'); caxis([-1 1]); axis square;
    subplot(2,2,3);
    H{3} = imagesc(Pcorr); colormap(jet);
    title('pattern corr'); caxis([-1 1]); axis square;
    if nargin > 5
        conn = Wmg*Wgm;
        subplot(2,2,4);
        H{4} = imagesc(conn-diag(NaN*diag(conn)));
        title('Connection'); axis square;
    end
end % setup_Pplot