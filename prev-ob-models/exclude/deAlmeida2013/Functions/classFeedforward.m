%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classFeedforward includes all proprieties and methods specific for
% feedforward cells.
%
% Licurgo de Almeida
% 04/29/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classFeedforward < classSpkNeuron
    properties
        % from Mitral cells
        AMPAFf % struct with parameters from the excitatory synapses
        ConnAMPAFf  = 0.4; % percentage of mitral cells connected to each
        % granule
        MAMPAFf % connection matrix with the mitral cells
        WAMPAFf % synaptic weight matrix with the mitral cells
    end
    methods
        function obj = classFeedforward(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 15; %ms
            obj.CellName = 'feedforward';
            obj.AMPAFf = struct('E',70e-3,'tau1',1,'tau2',2,'G',0.38);
            % where the elements of the struct are:
            % E: reversal potential
            % tau1: rising time of the conductance
            % tau2: falling time of the conductance
            % G: max conductance
            obj.MAMPAFf = obj.SetConnections(obj.ncells,obj.ConnAMPAFf);
            obj.WAMPAFf = obj.MAMPAFf; % if there' no learning, the
            % synaptic weights between connections are either 0 or 1
        end
    end
end