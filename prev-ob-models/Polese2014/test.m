%% Script to test the olfactory bulb model
percCompl_part = 0;
percCompl_total = 0;

I = ones(1,5000);
n_glomeruli = 16;


%% odor patterns
OdorA = [[32.5889474557272;36.2316774830248;5.07947265174024;36.5350342455608;25.2943698490164;3.90161619997638;11.1399287546819;21.8752607681994;38.3002734173719;38.5955414079711;6.30452326710193;38.8237112704246;38.2866779297178;19.4150259489136;32.0112187555520;5.67545354508861;]];

% Stimulation current 
I_ext{1} = OdorA*I;


%% Variables initialization
n_measure = length(I_ext);

Out_M1v    = cell(1,n_measure);
Out_T1v    = cell(1,n_measure);


%% initialization of glomeruli
for j = 1:n_glomeruli
    GLOM(j) = glomerulus(0, n_glomeruli - 1);
    InitGlo(GLOM(j));
    indx{j} = [1:n_glomeruli];
    indx{j}(j) = [];
end


for k  = 1:n_measure
    l_measure = length(I_ext{labindex + (k-1)});
    %% Total completetion Percentange
        if(labindex + (k-1) >= percCompl_total*n_measure/100)
            fprintf('Total completion percentage is: %3.2f %% \r',percCompl_total)
            percCompl_total=percCompl_total+10;
            if(percCompl_total>100)
                percCompl_total = 0;
            end
        end
        
        
    Out_M1v = zeros(n_glomeruli,l_measure);   
    Out_T1v = zeros(n_glomeruli,l_measure); 
    
  
    
    for i = 1:l_measure
        
        
        %% Partial completetion Percentange
        if(i >= percCompl_part*l_measure/100)
            fprintf('Partial completion percentage is: %3.2f %% \r',percCompl_part)
            percCompl_part=percCompl_part+10;
            if(percCompl_part>100)
                percCompl_part = 0;
            end
        end
        %%
        for j = 1:n_glomeruli
            stat_inib(j) = GLOM(j).out_inib;
            flag_inib(j) = GLOM(j).spike_flag_out_inib;
        end
        
        
        for j = 1:n_glomeruli
            State(GLOM(j),[], [], stat_inib(indx{j}), flag_inib(indx{j}),I_ext{labindex + (k-1)}(j,i)+15);
            
            Out_M1v(j,i) = GLOM(j).vout;
            Out_T1v(j,i) = GLOM(j).T1.v;

            Out_P1v(j,i) = GLOM(j).P1.v;
            Out_sSA1v(j,i) = GLOM(j).sSA1.v;	
            
        end
    end
end
save('Results.mat')
