function H = setup_Infoplot(time_axis,VST,line_style,line_name,Scorr,Pcorr,N)
    figure(N); clf;
    H = cell(1,length(VST)+1);
    subplot(1,2,2);
    H{1} = plot(uptri_1d(Scorr), uptri_1d(Pcorr), 'x');
    hold on; line([0 0],[-1 1],'color',[0 0 0]);
    hold on; line([-1 1],[0 0],'color',[0 0 0]);
    hold on; line([-1 1],[-1 1],'color',[0 0 0]);
    hold off;
    xlabel('input'); ylabel('output'); title('in/out corr');  axis square;
    subplot(1,2,1);
    xlim([min(time_axis) max(time_axis)]); ylim([-1 1]);
    xlabel('time'); ylabel('correlation/focality'); title('correlation');  axis square;
    subplot(1,2,1);
    for i = 1:length(VST)
        hold on;
        H{i+1} = plot(time_axis,VST{i},line_style{i});
    end
    hold off;
end