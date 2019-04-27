%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classGranule includes all proprieties and methods specific for granule
% cells. This class is used to create other spiking neurons, like
% FeedForward cells and Feedback cells in the cortical model.
%
% Licurgo de Almeida
% 04/22/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classGranule < classSpkNeuron
    properties
        % from Mitral cells
        AMPAFf % struct with parameters from the excitatory synapses
        ConnAMPAFf  = 0.4; % percentage of mitral cells connected to each
        % granule
        MAMPAFf % connection matrix with the mitral cells
        WAMPAFf % synaptic weight matrix with the mitral cells
        
        % from other Granule cells
        GABAFb % struct with parameters from the inhibitory synapses
        ConnGABAFb  = 0.08; % percentage of granule cells connected to each
        % granule
        MGABAFb % connection matrix with other granule cells
        WGABAFb % synaptic weight matrix with other granule cells
        
        % NE receptors parameters
        Tminne = [0e-3,0e-3]; % effects of NE modulation in the resting 
        % potential, represents the range THETA min in eq. The order here
        % is [Mod OFF, Mod ON];
        Kne = 3; % ligand concentration producing half occupation (Hill eq.)
        % See the funtcion CalculateAffinity.m
        % for details.
        Cmaxne = 10; % concentration range of the ligant (Hill eq.). See the
        % funtcion CalculateAffinity.m for details.
        bne = 1; % non-linearity of the modulation See the funtcion
        % CalculateAffinity.m for details.
    end
    methods
        function obj = classGranule(tsim,ncells,Ninput)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 15; %ms
            obj.CellName = 'Granule';
            if nargin < 3
                Ninput = obj.ncells; % if the number of cells from the
                % other network is not provided, the program assumes both
                % networks have the same number of neurons
            end
            obj.AMPAFf = struct('E',70e-3,'tau1',1,'tau2',2,'G',0.38);
            obj.GABAFb = struct('E',-15e-3,'tau1',4,'tau2',8,'G',0);
            % where the elements of the struct are:
            % E: reversal potential
            % tau1: rising time of the conductance
            % tau2: falling time of the conductance
            % G: max conductance
            obj.MAMPAFf = obj.SetConnections(Ninput,obj.ConnAMPAFf);
            obj.WAMPAFf = obj.MAMPAFf; % if there' no learning, the
            % synaptic weights between connections are either 0 or 1
            obj.MGABAFb = obj.SetConnections(obj.ncells,obj.ConnGABAFb,'auto');
            obj.WGABAFb = obj.MGABAFb;
        end
    end
end