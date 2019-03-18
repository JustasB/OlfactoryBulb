%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classFeedback includes all proprieties and methods specific for
% feedback cells.
%
% Licurgo de Almeida
% 04/29/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classFeedback < classSpkNeuron
    properties
        % from Mitral cells
        AMPAFf % struct with parameters from the excitatory synapses
        ConnAMPAFf  = 0.2; % percentage of pyr cells connected to each
        % feedback cell
        MAMPAFf % connection matrix with the mitral cells
        WAMPAFf % synaptic weight matrix with the mitral cells
                % from learning
        Tau11 = 200; %ms
        Tau01 = 200; %ms
        Tau10 = 200; %ms
        Taupost = 2; %ms
        Tauf = 7; %ms
        Taur = 1; %ms
        Tdelay = 1; %ms
    end
    methods
        function obj = classFeedback(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 15; %ms
            obj.CellName = 'feedback';
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