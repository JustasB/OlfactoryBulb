%
%
%

load('pars.mat');
N   = all_pars(1);

if exist('raster.x', 'file'),   
 raster = load('raster.x'); 
 if (~isempty(raster))

 tmax  = max(raster(:,1));
 delta = 50.;
 edges = 0:delta:tmax;
 [y, x] = histc(raster(:,1), edges); 

 y = 1000 * y / (N * delta);
 B = bar(edges, y, 'histc');
 set(B, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);

 xlabel('time [msec]');
 ylabel('Mean Firing Rate [Hz]');
 xlim([0 max(raster(:,1))]);
 end
end