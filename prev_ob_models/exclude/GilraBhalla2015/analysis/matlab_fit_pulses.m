function matlab_fit_pulses()

    kernels = [];
    fitted_mitrals = [0 1];
    % ========================================================================
    figure;
    hold on;
    for fitted_mitral = fitted_mitrals
        py = load(strcat('mit',int2str(fitted_mitral),'.mat'));

        for odornum = [0 1];
            % NOTE: matlab indices start with 1,
            % so careful with pulsenums and start_i (python-convention)
            if odornum==0
                pulsenums = [3 5];
            else
                pulsenums = [4 6];
            end

            % matlab cannot add int64-s, so convert to int32
            starti = int32(py.start_i)+1;
            ydata = py.firingbinsmeanList(pulsenums-1,starti:end);
            len_data = size(ydata,2);
            xdata = py.pulseList(pulsenums,starti:starti+len_data-1);
            start_kernel = py.start_kernels(odornum+1,:);
            %start_kernel = zeros(size(start_kernel));

            % Copied from Priyanka
            % ========================================================================
            % initialize and define the optimization function
            % ========================================================================
            Eval_max = 1e+6; Iter_max = 1e+6; lb = []; ub = [];
            options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max,...
                    'TolFun',1e-15,'TolX',1e-15);
            [kernel,resnorm,residual,exitflag] = ...
                    lsqcurvefit(@conv_adi,start_kernel,...
                    xdata,ydata,lb,ub,options);
            % normalize chi-sq to the number of dof-s
            num_dof = prod(size(ydata)) - prod(size(kernel));
            resnorm/num_dof
            exitflag
            % extra iteration seems to make no difference in resnorm
            %randnoise = rand(size(kernel))*5;
            %[kernel,resnorm,residual,exitflag] = ...
            %        lsqcurvefit(@conv_adi,kernel.*randnoise,...
            %        xdata,ydata,lb,ub,options);
            %resnorm
            %exitflag

            subplot(3,2,odornum*2+fitted_mitral+1);
            hold on;
            % plot only the 1st pulselist & response for this odor
            plot(xdata(1,:)*20,'g-');
            plot(ydata(1,:),'b-');
            [fit] = conv_adi(kernel,xdata(1,:));
            plot(fit,'r-');
            %[origfit] = conv_adi(start_kernel,xdata(1,:));
            %plot(origfit,'m-');

            kernels = [kernels; kernel];
        end
        
        subplot(3,2,5+fitted_mitral);
        hold on;
        AplusBresponse = py.firingbinsmeanList(6,starti:end);
        Apulse = py.pulseList(7,starti:starti+len_data-1);
        Bpulse = py.pulseList(8,starti:starti+len_data-1);
        AplusBpulse = Apulse+Bpulse;
        plot(AplusBpulse*20,'g-');
        plot(AplusBresponse,'b-');
        % bgnd gets added twice by conv_adi, so subtract once.
        [fit] = conv_adi(kernels(fitted_mitral*2+1,:),Apulse) ...
            + conv_adi(kernels(fitted_mitral*2+2,:),Bpulse) - py.bgnd;
        plot(fit,'r-');

    end
    
    % plot kernels
    % savitsky golay filtering (operates on col-s, hence ')
    kernels = sgolayfilt(kernels',4,11)'
    figure(2);
    subplot(2,1,1);
    plot(kernels(1,:),'r-');
    hold on;
    plot(kernels(3,:),'b-');
    subplot(2,1,2);
    plot(kernels(2,:),'r-');
    hold on;
    plot(kernels(4,:),'b-');

    % partly from Priyanka
    % ========================================================================
    % convolution and fit generation for optimization
    % ========================================================================
    function [fit] = conv_adi(kernel,xdata)
        fit = [];
        % convolve for each pulse data
        for i = 1:size(xdata,1)
            M = [];
            M = conv(kernel,xdata(i,:));
            M = M(1:len_data).*py.pulserebindt;
            M = M + py.bgnd;
            fit(i,1:length(M)) = M;
        end
        fit(fit<0) = 0;
    end

end