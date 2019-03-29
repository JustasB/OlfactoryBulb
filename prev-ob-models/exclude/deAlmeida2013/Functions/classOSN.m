%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classOSN includes all proprieties and methods specific for OSNs
%
% Licurgo de Almeida
% 04/17/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classOSN < classNoSpkNeuron
    properties
        AMPAFf % struct with parameters from the excitatory synapses. We
        % are using the same struct we use in other neurons to help in the
        % implementation
        OdorDistM % mean of the input distribution when the inputs are not
        % randomyzed
        OdorDistSigma % standard deviation of the input distribution
        FlagRandomDist = false; % if true, the input distributions are
        % randomized
        RespirationFreq = 0; % if 0, the respiration frequency is off
        OdorInput % input to OSNs
        InputTimes % beginning and ending of the odor input
        MInput % connection matrix for the inputs and the OSNs
        WInput % weight matrix for the inputs and the OSNs
    end
    properties(Hidden = true)
        DefaultM = 0.5; % by default, the most active neurons is the one
        % right in the middle of the raster plot
        DefaultSigma = 0.1; % default sigma proportion to the size of the
        % network;
    end
    methods
        function obj = classOSN(tsim,ncells,OdorDistM,OdorDistSigma)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classNoSpkNeuron(inputsuper{:});
            if nargin > 2
                obj.OdorDistM = OdorDistM;
                obj.OdorDistSigma = OdorDistSigma;
            else
                obj.OdorDistM = ceil(obj.ncells * obj.DefaultM); %makes the
                % most active neuron the one right in the center of the
                % raster plot
                obj.OdorDistSigma = ceil(obj.ncells * obj.DefaultSigma);
                % arbitrary value for the standard deviation
            end
            obj.Tmax = [10e-3,10e-3];
            obj.InputTimes = [0,obj.tsim];
            obj.OdorInput = SetOdorInput(obj);
            obj.AMPAFf = struct('E',70e-3,'G',0.16);
            obj.MInput = eye(obj.ncells);
            obj.CellName = 'OSN';
            obj.WInput = obj.MInput;
        end
        function OdorInput = SetOdorInput(obj,OdorVector,odorC)
            % This function creates the odor input to OSNs
            if nargin < 3
                odorC = 1;
            end
            R = CreateRespFreq(obj);
            x = 1:obj.ncells;
            InputVector = normpdf(x,round(obj.ncells / 2),obj.OdorDistSigma);
            InputVector = InputVector' / max(InputVector);
            if nargin < 2
                if obj.FlagRandomDist == true
                    InputVector = InputVector(randperm(length(InputVector)));
                else
                    InputVector = SetAffinity(obj,InputVector);
                end
            else %if the odor pattern vector is provided
                InputVector = InputVector(OdorVector);
            end
            OdorInput = InputVector * R * odorC;
        end
    end
    methods(Hidden = true)
        function R = CreateRespFreq(obj)
            % This function creates an artificial respiration to modulate
            % the OSN activity
            time = obj.dt:obj.dt:obj.tsim;
            time = time / 1000; %converting the time to seconds
            R = -cos(2 * pi * obj.RespirationFreq * time);
            if obj.RespirationFreq == 0
                R = abs(R);
            else
                R = R + 1;
                R = R / 2;
            end
            inputtimes = obj.InputTimes / 1000; %converting the time to seconds;
            I = time < inputtimes(1); % begin input
            R(I) = 0;
            I = time > inputtimes(2); % end input
            R(I) = 0;
        end
        function Input = SetAffinity(obj,Input)
            % This function rotates the odor affinity of OSNs
            if obj.OdorDistM ~= round(obj.ncells / 2)
                Input = circshift(Input,obj.OdorDistM - round(obj.ncells / 2));
            end
        end
    end
end
