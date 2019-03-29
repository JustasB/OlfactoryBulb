classdef NeuronIzh < handle
    %NEURONIZH 
    %   The model of neurons is modeled with following equation
    %   C*v' = a v^2 + b v + c - u + I
    %   u' = d (e v - u)
    %   where:
    %   v represent membrane potential of neuron
    %   u represent membrane recovery variable
    %   a, b, c, d, e are parameter that change type of neuron
    %   I is input synaptic current
    %%
    %  Written by Davide Polese davide.polese@artov.imm.cnr.it
    %  Consiglio Nazionale delle Ricerche - Istituto per la Microelettronica e Microsistemi
    %  Via del Fosso del Cavaliere, 100 00133 Roma
    
    properties(SetAccess = private)
        v;                      % membrane potential of neuron
        u;                      % membrane recovery variable
        I;                      % input current
        t_step;                 % step time of simulation
        t;                      % current time
        
    %% PARAMETER MODEL
        a;  
        Cap;                    % Membrane capacitance
        d;
        e;
        V_peak;                 % spike peak voltage
        V_th;                   % Threshold voltage of neuron
        v_rest;                 % Resting voltage of neuron
        v_rep;                  % Repolarising Voltage
        v_inib;                 % Inhibition voltage
        
    %% PARAMETER OUTPUT
        out_state;              % State of output current
        Iout;                   % Output current
        tau_m_in;               % time costant decaying of inhibition kernel
        tau_s_in;               % time costant risiing of inhibition kernel
        tau_m;                  % time costant decaying of kernel 
        tau_s;                  % time costant rising of kernel 
        
     %% FLAG
        spike_flag;             % it is 1 if neurons make a spike
        
        n_neuron;               % number ID of neuron
        in_neurons;             % numbers ID of input neurons
        n_spike;                % number of spike
        
      %% weight
        wij                     % Synapse weight
        
    end
    
    properties (Dependent = true, SetAccess = private)
        ker;                    % kernel of current output
        
        %% PARAMETER MODEL
        b;
        c;
    end
   
    methods
        function ker = get.ker(NI)
            ker = kernel([0:NI.t_step:5*NI.tau_m], 0 , NI.tau_m_in, NI.tau_s_in, NI.tau_m, NI.tau_s);
        end
        
        function b = get.b(NI)
            b = -NI.a*(NI.v_rest+NI.V_th);
        end
        
        function c = get.c(NI)
            c =NI.a*NI.v_rest*NI.V_th;
        end
            
    
    
        function NI = NeuronIzh(PAR, IC, T_STEP, T0, V_PEAK, V_TH, V_REST, V_REP, V_INIB, TAU, ID, IN_NEURONS)
            % Inizialization of parameter
            NI.a = PAR(1);
            NI.Cap = PAR(2);
            NI.d = PAR(3);
            NI.e = PAR(4);
            NI.V_peak = V_PEAK;
            NI.V_th = V_TH;            
            NI.v_rest = V_REST;
            NI.v_inib = V_INIB;
            NI.v_rep = V_REP;
            
            % Initial condition
            NI.v = IC(1);
            NI.u = IC(2);
            NI.t_step = T_STEP;
            NI.t = T0;
            NI.I = 0;
            
            % out parameter
            NI.out_state = 0;
            NI.Iout = 0;                   
            NI.tau_m_in = TAU(1);                  
            NI.tau_s_in = TAU(2);                  
            NI.tau_m = TAU(3);                  
            NI.tau_s = TAU(4);                  
            
            % flag
            NI.spike_flag = 0;
            
            NI.n_neuron = ID;
            NI.in_neurons = IN_NEURONS;
            NI.n_spike = 0;
            
        end
        
        function Init_wij(NI,SINAPSI,WIJ)
            NI.wij = [];
            if(~isempty(SINAPSI))
               if (isempty(WIJ))    
                   NI.wij=40*ones(1,length(SINAPSI));
                   NI.wij = NI.wij .* sign(SINAPSI);            
               else
                   NI.wij = SINAPSI.*WIJ;
               end
            end
        end
            
        function PotentialNI(NI, Input, I_ext, GLO)
            %global SpikeMatrix
            %% Conductance of synapse
            g = 0.1;
            
            %% reset flag
            NI.spike_flag = 0;
            
            if (NI.t < NI.t_step)
                NI.out_state = zeros(1, length(NI.ker));
            end
                
            
            %% Membrane potential computation
            
            %% Current Input
            I_syn = Input * NI.wij';
            if(isempty(I_syn))
                I_syn = 0;
            end
            
            NI.I = I_syn + I_ext;
            %% Izhikevich model
            
            function dy = Izh(y)
               dy = zeros(2,1);    % a column vector
               dy(1) = (NI.a*y(1)^2 + NI.b*y(1) + NI.c - y(2) + NI.I)/NI.Cap;
               dy(2) = NI.d*(NI.e*(y(1)- NI.v_rest) - y(2));
            end
            %% Solution of next time
            Y = RK4(@Izh, 0.1, [NI.v NI.u]);
            NI.v = Y(1);
            NI.u = Y(2);
            if(NI.v >= NI.V_peak)
                NI.v = NI.v_rep;
                NI.u = NI.u + NI.v_inib;
                NI.spike_flag = 1;
                NI.n_spike = NI.n_spike + 1;
               set_SpikeMatrix(GLO, NI.n_neuron,NI.n_spike, NI.t);
            end
            
            %% CURRENT OUTPUT
            NI.out_state = [NI.out_state(2:end),0] + NI.spike_flag * NI.ker;
            NI.Iout = NI.out_state(1);
            
            %% Increment time
            NI.t = NI.t + NI.t_step;
        end      
    end
end

