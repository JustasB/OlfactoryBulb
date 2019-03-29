%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots a rasterplot of the neurons
%
% INPUTS:
% Spikes - ncells x tp matrix of spikes
% dt     - time step (ms)
% tsim   - simulation time (ms)
% color  - plot color
% fs     - fontsize
% dl     - draw labels (0 - NO, 1 - YES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RasterPlot(Spikes,dt,tsim,color,fs,dl)


hold on;

for ii = 1:size(Spikes,1)
    J = find(Spikes(ii,:));
    for jj = 1:length(J)
        spkx = [J(jj),J(jj)] .* dt;
        spky = [ii,ii + 0.9];
        line(spkx,spky,'color',color,'LineWidth',1);
    end
end
set(gca,'fontsize',fs)

axis([0,tsim + dt,0,size(Spikes,1) + 2]);
if dl == 1
    xlabel('time (ms)','fontsize',fs);
    ylabel('neuron','fontsize',fs);
end