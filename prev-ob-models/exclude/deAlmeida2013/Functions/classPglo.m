%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classPglo includes all proprieties and methods specific for
% periglomerular cells
%
% Licurgo de Almeida
% 04/19/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classPglo < classNoSpkNeuron
    properties
        AMPAFf % struct with parameters from the excitatory synapses
        MAMPAFf % connection matrix with the OSNs
        WAMPAFf % % synaptic weight matrix with the OSNs
    end
    methods
        function obj = classPglo(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classNoSpkNeuron(inputsuper{:});
            obj.tau = 2; %ms
            obj.AMPAFf = struct('E',70e-3,'G',0.17);
            % where the elements of the struct are:
            % E: reversal potential
            % G: max conductance
            obj.MAMPAFf = eye(obj.ncells);
            obj.WAMPAFf = obj.MAMPAFf; % if there' no learning, the
            % synaptic weights between connections are either 0 or 1
            obj.CellName = 'Pglo';
        end
    end
end