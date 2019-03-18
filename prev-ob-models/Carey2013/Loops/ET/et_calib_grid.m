function et_calib_grid(input)

    inputgain = 200; % Not divided by ET cell capacitance
    
    %input is integer between 0 and 4   
    inputtrace = ['orn_inputs_et_calib_' num2str(input+1) 'hz'];
    
    doloop_ET_grid(inputgain, inputtrace)
    
end
