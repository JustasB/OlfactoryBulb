 classdef glomerulus < handle
    %% GLOMERULUS
    % This class implements a glomerulus.
    % The glomerulus have four neurons Mitral Cell, Periglomerular Cell 
    % External Tufted Cell and superficial Short Axons.  
    %%
    %  Written by Davide Polese davide.polese@artov.imm.cnr.it
    %  Consiglio Nazionale delle Ricerche - Istituto per la Microelettronica e Microsistemi
    %  Via del Fosso del Cavaliere, 100 00133 Roma
    
    properties(SetAccess = private)
        
        N_input_spike = 0;
        n_input_ORN = 0;
        n_lateral_input = 0;
        %% Matrix of connection include synaptic type
        % -1 Inhitory synapse
        %  0 No synapse
        % +1 Excitatory synapse
        MatrixConnection = [0, -1, 0, 0, 0;       %MC PG ET sSA GloIn
                            1,  0, 1, 1, 0;       %PG 
                            0,  0, 0, 0, 1;       %ET
                            0,  0, 1, 0, 0];      %sSA
        
        M1;
        P1;
        T1;
        sSA1;
        
        
        out = 0;
        out_inib = 0;
        spike_flag_out = 0;
        spike_flag_out_inib = 0;
        vout = 0;
        
        %% Matrix containing time of spike
        SpikeMatrix 
    end
    
    
    properties (Dependent = true, SetAccess = private)
        sinapsiM1;
        sinapsiP1;
        sinapsiT1;
    end
    
    methods
        
        function set_SpikeMatrix(GLO, row, col, value)
            GLO.SpikeMatrix(row,col) = value;
        end
        
        function sinapsiM1 = get.sinapsiM1(GLO)
                        % ORN, PG, ET
            sinapsiM1 = [ones(1,GLO.n_input_ORN),-1, 1]; 
        end
        
        function sinapsiP1 = get.sinapsiP1(GLO)
                        %sSA MC ET 
            sinapsiP1 = [ones(1,GLO.n_input_ORN), ones(1,GLO.n_lateral_input), +1, +1];  
        end
        
        function sinapsiT1 = get.sinapsiT1(GLO)
            sinapsiT1 = [ones(1,GLO.n_input_ORN), ones(1,GLO.n_lateral_input)];
        end
        
        function GLO = glomerulus(N_INPUT, N_INIBIZIONE)
            if nargin==0
                N_INPUT = 0;
                N_INIBIZIONE = 0;
            end
            GLO.n_input_ORN = N_INPUT;
            GLO.n_lateral_input = N_INIBIZIONE;
        end
        
        function InitGlo(GLO)
            %% Initialization of neurons        
            % N = NeuronIzh([a, d, e ],  [v(0), u(0)], t_step, t0, V_peak, V_th, V_rest, V_inib, [tau1, tau2], N_neuron, input_neuron) 
            GLO.M1 = NeuronIzh([1, 40, 0.4, 2.6], [-55, 0], 0.1, 0, 35, -50, -55, -50, 200,[0 0 10,2],1, [2,3,5:(5+GLO.n_input_ORN - 1)]);
            GLO.P1 = NeuronIzh([0.049, 5.9, 0.0167, -0.94], [-53.1, 0], 0.1, 0, 35, -20, -53.1, -20, 50, [0 0 10, 2],2, [1, 3, 5:(5+GLO.n_input_ORN + GLO.n_lateral_input - 1)]);
            GLO.T1 = NeuronIzh([1, 40, 0.4, 2.6], [-55, 0], 0.1, 0, 35, -50, -55, -50,200,[0 0 10,2],3, [5:(5+GLO.n_input_ORN + GLO.n_lateral_input - 1)]);
            GLO.sSA1 = NeuronIzh([0.061 58 0.049 -0.68], [-67 0], 0.1, 0, 35, -30, -67, -30, 150, [5 1 10 2], 4, [3]);
            %% Initialization of weight
            Init_wij(GLO.M1,GLO.sinapsiM1, [40 0] );
            Init_wij(GLO.P1,GLO.sinapsiP1, [40*ones(1,15), 40 40 ] );
            Init_wij(GLO.T1,GLO.sinapsiT1, [40*ones(1,15)] );
            Init_wij(GLO.sSA1, [+1], [19.5]);
            
            GLO.SpikeMatrix(1:4+GLO.n_input_ORN+GLO.n_lateral_input,1:100)=-1000;
        end
       

        function State(GLO, Input, Input_spike, Input_lateral, Input_lateral_spike, I_ext)
            %% updating Spike matrix
            GLO.N_input_spike = GLO.N_input_spike + [Input_spike, Input_lateral_spike];
            ind = find([Input_spike, Input_lateral_spike]);                           % index of neuron spiking
            GLO.SpikeMatrix([4+ind], GLO.N_input_spike(ind) ) = GLO.M1.t;

            
            %% Calculation of potential
            PotentialNI(GLO.T1,[Input, Input_lateral],I_ext, GLO);
            PotentialNI(GLO.sSA1,[GLO.T1.Iout], 0, GLO);
            PotentialNI(GLO.P1,[Input, Input_lateral, GLO.M1.Iout, GLO.T1.Iout], I_ext, GLO);
            PotentialNI(GLO.M1,[Input,GLO.P1.Iout, GLO.T1.Iout],1*I_ext, GLO);
            
            %% Output of Glomerulo
            GLO.spike_flag_out = GLO.M1.spike_flag;
            GLO.spike_flag_out_inib = GLO.sSA1.spike_flag;
            GLO.out_inib = GLO.sSA1.Iout;
            GLO.out = GLO.M1.Iout;
            GLO.vout = GLO.M1.v;
        end
    end
end
