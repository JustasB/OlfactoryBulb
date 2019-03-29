function [H, Iaxis] = setup_Iplot(cont_density,t_axis,I,Icorr,Iage,Wmg,N_popu,N)
    Ns = size(I,2);
    H = cell(1,4);
    figure(N); clf; set(gcf,'DoubleBuffer','on');
    subplot(2,1,1);
    H{1} = imagesc(I); colormap(jet);
    Iaxis = gca;
    xlabel('stimulus'); ylabel('GC'); title('I');  axis square;
    subplot(2,1,2);
    H{2} = imagesc(Icorr); colormap(jet);
    title('pattern corr'); xlim([0.5 Ns+.5]); ylim([0.5 Ns+.5]); caxis([-1 1]); axis square;
end % setup_Iplot