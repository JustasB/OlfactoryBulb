%
%
%

% load('pars.mat');
% N   = all_pars(1);
% 
% if exist('raster.x', 'file'),   
%  raster = load('raster.x'); 
%  if (~isempty(raster))
%   tmax  = max(raster(:,1));
%   rates = zeros(N,1);
%   for k=1:N,
%    rates(k) = length(find(raster(:,2) == (k-1)));
%   end
%   rates = 1000* rates / tmax; 

 rates = load('rates.x');  
 delta = 5.;
 edges = 0:delta:max(rates);
 [y, x] = histc(rates, edges); 

  %y = y / N;
  B = bar(edges, y, 'histc');

%  B = bar(1:N, rates,1);
  set(B, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);

  xlabel('Firing Rate [Hz]');
  ylabel('Fraction of neurons');
  
%  xlim([1 max(rates)]);
  xlim([0 max(rates)]);
  ylim([0 length(rates)/5.]);

