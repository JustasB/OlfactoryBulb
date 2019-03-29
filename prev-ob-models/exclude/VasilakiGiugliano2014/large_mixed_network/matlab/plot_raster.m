%
%
%

load('pars.mat');
N   = all_pars(1);

if exist('raster.x', 'file'),   
 raster = load('raster.x'); 
 if (~isempty(raster))
 cla;
 for k=1:size(raster,1),
     tsp = raster(k,1);
     ind = raster(k,2);
     if ((N>50) && (ind <50))
      L = line([1 1] * tsp,[0 1] + ind);
      set(L, 'Color', [0 0 0], 'LineWidth', 1);
     elseif (N<50)
      L = line([1 1] * tsp,[0 1] + ind);
      set(L, 'Color', [0 0 0], 'LineWidth', 1);
     end
 end
 

 xlabel('time [msec]');
 ylabel('Network unit (first 50)');
 xlim([0 max(raster(:,1))]);
 end
end