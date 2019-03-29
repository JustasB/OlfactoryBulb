%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classMitsoma includes all proprieties and methods specific for the soma
% compartment of mitral cells. Mitral cells are divided in two compartments
%
% Licurgo de Almeida
% 04/22/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classMitsoma < classSpkNeuron
    properties
        % from apical dendrite
        Rglo = 1; % core resistance from the apical compartment.
        
        % from Granule cells
        GABAFb % struct with parameters from the inhibitory synapses
        ConnGABAFb = 0.4; % percentage of granule cells connected to each
        % mitral cell
        MGABAFb % connection matrix with the granule cells
        WGABAFb % synaptic weight matrix with the granule cells
    end
    methods
        function obj = classMitsoma(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 20; %ms
            obj.CellName = 'Mitsoma';
            obj.GABAFb = struct('E',-15e-3,'tau1',4,'tau2',8,'G',0.38);
            % where the elements of the struct are:
            % E: reversal potential
            % tau1: rising time of the conductance
            % tau2: falling time of the conductance
            % G: max conductance
            obj.MGABAFb = obj.SetConnections(obj.ncells,obj.ConnGABAFb);
            obj.WGABAFb = obj.MGABAFb; % if there' no learning, the
            % synaptic weights between connections are either 0 or 1
        end
    end
end