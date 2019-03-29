function y = hpf(x,Fs,order,freq)
    [b1,a1] = butter(order,freq/(Fs/2),'high');   %design butterworth filter coefficients
    y = filtfilt(b1,a1,x);