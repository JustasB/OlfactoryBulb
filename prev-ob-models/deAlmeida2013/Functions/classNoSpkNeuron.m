%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classNoSpkNeuron has some general proprieties methods for non-spiking
% neurons
%
% Licurgo de Almeida
% 04/17/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef classNoSpkNeuron < classNetGeneral
    properties
        PFire % output probability
    end
    methods
        function obj = classNoSpkNeuron(tsim,ncells)
            if nargin == 0
                input = {};
            else
                input = {tsim,ncells};
            end
            obj = obj@classNetGeneral(input{:});
            obj.PFire = zeros(obj.ncells,1);
        end
    end
end