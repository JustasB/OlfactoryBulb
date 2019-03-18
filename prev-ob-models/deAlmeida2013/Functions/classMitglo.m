%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classMitglo includes all proprieties and methods specific for the apical
% dendrite of mitral cells. Mitral cells are divided in two compartments
%
% Licurgo de Almeida
% 04/19/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classMitglo < classNoSpkNeuron
    properties
        % from OSNs
        AMPAFf % struct with parameters from the excitatory synapses
        MAMPAFf % connection matrix with the OSNs
        WAMPAFf % synaptic weight matrix with the OSNs
        
        % from Pglo cells
        GABAFf % struct with parameters from the inhibitory synapses
        MGABAFf % connection matrix with the Pglo cells
        WGABAFf % synaptic weight matrix with the Pglo cells
        
        % from soma dendrite
        Rsom = 1; % core resistance from the soma compartment.
    end
    methods
        function obj = classMitglo(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classNoSpkNeuron(inputsuper{:});
            obj.tau = 5; %ms
            obj.CellName = 'Mitglo';
            obj.AMPAFf = struct('E',70e-3,'G',0.27);
            obj.MAMPAFf = eye(obj.ncells);
            obj.WAMPAFf = obj.MAMPAFf; % if there' no learning, the
            % synaptic weights between connections are either 0 or 1
            obj.GABAFf = struct('E',-15e-3,'G',0.38);
            % where the elements of the struct are:
            % E: reversal potential
            % G: max conductance
            obj.MGABAFf = eye(obj.ncells);
            obj.WGABAFf = obj.MGABAFf;
        end
    end
end