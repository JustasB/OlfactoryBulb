function y = epsc(t,samplingrate)
    tau1 = 0.045 * samplingrate; %convert from seconds to samples
    tau2 = 0.035 * samplingrate;
    y = (exp(-t./tau1) - exp(-t./tau2))./(tau1 - tau2);