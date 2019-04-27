function deconv_trs = deconvolve(trs,samplingrate,tau)
    %DECONVOLVE  Deconvolves an exponential from the given trace
    % Trace should be divided by the resting light (RLI) prior to calling this function.
    % Example usage:
    %    tr = o.trials(1).rois.traces(glom_num,:) ./ o.trials(1).rois.RLIs(glom_num);
    %    deconv_tr = deconvolve(tr,o.trials(1).rois.samplingrate);
    %

    if nargin < 3
        tau = 0.248;
    end
    
    butter_Wn = 7.5; %Hz 7.5?
    
    
    
    %butter_a = 0; butter_b = 0; kernelpad = 0; kernel = 0;

    kernellength = round(2 * tau * samplingrate) - 1;

    kernelpad = zeros(1, kernellength);
    kernel = exp((0:-1:-kernellength) / (tau * samplingrate));
    
    [butter_b, butter_a] = butter(4, butter_Wn / (samplingrate / 2));
    
    num = size(trs,1);
    deconv_trs = trs;
    
    for n = 1:num
        filt_tr = filtfilt(butter_b,butter_a,trs(n,:));
        deconv_trs(n,:) = deconv([filt_tr kernelpad], kernel);
    end
end