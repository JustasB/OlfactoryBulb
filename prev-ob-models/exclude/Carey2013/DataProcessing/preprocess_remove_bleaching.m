function traces = preprocess_remove_bleaching(rawtraces)
    
	crv_tc = {1:100 801:900 1301:1400 1751:1850 2201:2300}; %spaces in between bouts
	crv_t = cellfun(@mean,crv_tc); %[crv_tc{:}];
	
    [numrows, trlen] = size(rawtraces);
    traces = rawtraces;
	threshes = zeros(1,numrows);
    baselines = zeros(1,numrows);
    
	for row = 1:numrows
		dec = rawtraces(row,:);%deconvolve(dec',100,0.325);
		dec(1) = dec(2); %first point usually shows filtering artifact
        dec(end-25:end) = dec(end-25); %because it's easier to just do it this way
		
		%remove bleaching:
		%calc exp. fit to inter-bout intervals:
		%expcfit = fit(crv_t',dec(crv_t)','exp2'); %going directly doesn't work all the time
		crv_binned = cellfun(@mean,subscell(dec,crv_tc));
		
        [min1hz_y, min1hz_t] = min(dec(250:650));
        min1hz_t = min1hz_t + 250;
        
        fit_t = [crv_t(1) min1hz_t crv_t(2:end)];
        fit_y = [crv_binned(1) min1hz_y crv_binned(2:end)];
        
		opt = optimset('Display','off','DiffMinChange',1.0e-8,'DiffMaxChange',0.1,'MaxFunEvals',600,'MaxIter',400);
	    expcfit = lsqcurvefit(@expcurve,[0.001 0.002 0 0], fit_t, fit_y, [0 0 -Inf -Inf], [Inf Inf Inf Inf], opt);
        
        %figure(row),plot(dec),hold all, plot(expcurve(expcfit,1:trlen))
        
		decf = dec - expcurve(expcfit,1:trlen); %subtract fit
        decf = smooth(decf')';
		
		threshes(row) = mean(decf([crv_tc{:}])) + std(decf([crv_tc{:}])) .* 1;
        baselines(row) = mean(decf([crv_tc{:}]));
        
		%decf(decf < threshes(row)) = baselines(row);
        %decf([crv_tc{:}]) = baselines(row);
        decf = smooth(decf',15)';

        decf = deconvolve(decf,100,0.325);
		decf(1:100) = decf(100);
        decf(decf < 0) = 0;
        
        decf = hpf(decf,100,2,0.01); %this doesn't work?
        
        traces(row,:) = decf;
	end
    
end

function out = subscell(c,inds)
	out = inds;
	for n = 1:length(inds)
		out{n} = c(inds{n});
	end
end

function y = expcurve(A,x)
	y = A(1) .* exp(-A(2).*x + A(3)) + A(4);
end