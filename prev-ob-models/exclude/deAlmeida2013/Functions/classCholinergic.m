%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classCholinergic includes all proprieties and methods specific for
% the ACh cells in the HDB.
%
% Licurgo de Almeida
% 07/22/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef classCholinergic < classSpkNeuron
    properties
        % from Mitral cells
        AMPAFf % struct with parameters from the excitatory synapses
        ConnAMPAFf  = 0; % percentage of ACh cells connected to each
        % mitral
        MAMPAFf % connection matrix with the mitral cells
        WAMPAFf % synaptic weight matrix with the mitral cells
        
        % from Gabaergic interneurons
        GABAFf % struct with parameters from the inhibitory synapses
        ConnGABAFf  = 0.5; % percentage of Ff cells connected to each
        % gabaergic cell from the HDB
        MGABAFf % connection matrix with the Gabaergic interneurons
        WGABAFf % synaptic weight matrix with the Gabaergic interneurons
        
    end
    methods
        function obj = classCholinergic(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 10; %ms
            obj.CellName = 'Cholinergic';
            obj.AMPAFf = struct('E',70e-3,'tau1',1,'tau2',2,'G',0);
            % where the elements of the struct are:
            % E: reversal potential
            % tau1: rising time of the conductance
            % tau2: falling time of the conductance
            % G: max conductance
            obj.MAMPAFf = obj.SetConnections(obj.ncells,obj.ConnAMPAFf);
            obj.WAMPAFf = obj.MAMPAFf;
            obj.GABAFf = struct('E',-15e-3,'tau1',4,'tau2',8,'G',0.38);
            obj.MGABAFf = obj.SetConnections(obj.ncells,obj.ConnGABAFf);
            obj.WGABAFf = obj.MGABAFf;
        end
    end
end