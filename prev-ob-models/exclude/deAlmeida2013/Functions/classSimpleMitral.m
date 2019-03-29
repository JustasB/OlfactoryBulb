%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classSimpleMitral creates a simplified version of the mitral activity
%
% Licurgo de Almeida
% 04/25/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef classSimpleMitral
    properties
        CellName = 'Mitsimple';
        dt = 0.5; % timestep (in ms)
        tsim = 1000; % simulation time (in ms)
        ncells = 100; % number of neurons in the network
        tinit = 1; % stimulus begin (in ms)
        tfinal = 1000; % stimulus end (in ms)
        InputTimes % beginning and ending of the odor input
        TFire % time since last spike
        noise = 0; % noise level of non-active neurons
        f = 30; % frequency of modulation
        r = 60; % sets the average firing rate in response to the odor
        % stimulation
        delta = 16; % controls the duration of the positive excursion of
        % the oscillation. higher delta = higher synchronization.
        A = 4; % controls the sparseness of the odor representation.
        % higher A = low sparseness and vice-versa.
        odorant %reference number of the odorant to used in the 
        % save file.
        refrac = 20; % Period o inactivity to make the activation patter of
        % Mitral cells from this model similar to the biophysical model. This
        % should not be seeing as something like a refractory period of the
        % integrate and fire neurons.
        O % Output (0 or 1 for spiking neurons); between 0 and 1 for non-
        % spiking neurons
    end
    methods
        function obj = classSimpleMitral(tsim,ncells,odorant)
            if(nargin > 0)
                obj.tsim = tsim;
                obj.ncells = ncells;
                obj.tfinal = tsim;
            end
            if nargin > 2
                obj.odorant = odorant;
            else
                obj.odorant = round(obj.ncells / 2);
            end
            obj.TFire = ones(obj.ncells,1) * 10000000000;
            obj.InputTimes = [0,obj.tsim];
        end
        function obj = CreateOutput(obj)
            inputmat = CreateAffinity(obj);
            time = (obj.dt:obj.dt:obj.tsim) / 1000;
            R = ((1 - cos(2 * pi * (obj.f) * time)) / 2).^obj.delta;
            R = R / (sum(R) / (obj.tsim / 1000));
            R = R * obj.r;
            inputtimes = obj.InputTimes / 1000; %converting the time to seconds;
            I = time < inputtimes(1); % begin input
            R(I) = 0;
            I = time > inputtimes(2); % end input
            R(I) = 0;
            Act = inputmat * R;
            Z = rand(size(Act));
            Miraster = Act >= Z;
            Miraster = Miraster * 1.0;
            time = round(obj.refrac / obj.dt);
            T = size(Miraster,2);
            Out = obj.CreateRefPeriod(Miraster,T,time);
            % Add expontaneous activity to cells
            if obj.noise > 0
                N = rand(size(Act));
                noisec = obj.noise - (sum(Out,2) ./ (obj.tsim / 1000));
                I = noisec < 0;
                noisec(I) = 0;
                for ii = 1:obj.ncells
                    N(ii,:) = N(ii,:) <= noisec(ii) / (1000 / obj.dt); % noise
                end
                
                Out = Out | N;
                Out = Out * 1.0;
                Out = obj.CreateRefPeriod(Out,T,time);
            end
            obj.O = Out;
        end
        function S = Sparseness(obj)
            Vact = sum(obj.O,2);
            N = obj.ncells;
            S = (1 - ((sum(Vact / N)^2) / sum((Vact.^2) / N))) / (1 - 1 / N);
        end
        function C = Coherence(obj)
            % This function calculates the normalized coherence in the 
            % network activity by comparing it to a random network with
            % same firing rate
            twindow = 2; %ms
            sizeblock = twindow / obj.dt;
            NS = obj.O;
            Creal = obj.CalculateCoherence(NS,sizeblock);
            randNS = rand(size(NS));
            Act = sum(NS,2);
            Act = Act / size(NS,2);
            for ii = 1:length(Act)
                randNS(ii,:) = randNS(ii,:) <= Act(ii);
            end
            Crand = obj.CalculateCoherence(randNS,sizeblock);
            C = 1 - (Crand / Creal);
        end
        function F = Frequency(obj)
            % This function extracts the average frequency of a group of
            % cells
            F = sum(obj.O,2) / (obj.tsim / 1000);
        end
        function Raster(obj)
            % This function plots a rasterplot of the neurons
            scrsz = get(0,'ScreenSize');
            figH = figure;
            set(figH,'position',[0,400,scrsz(3)-0.4*scrsz(3),scrsz(4)-0.6*scrsz(4)]);
            cla;
            hold on;
            box on;
            title(obj.CellName,'fontsize',16);
            for ii = 1:size(obj.O,1)
                J = find(obj.O(ii,:));
                for jj = 1:length(J)
                    spkx = [J(jj),J(jj)] .* obj.dt;
                    spky = [ii,ii + 0.9];
                    line(spkx,spky,'color','k','LineWidth',1);
                end
            end
            axis([0,obj.tsim + obj.dt,0,size(obj.O,1) + 2]);
            xlabel('time (ms)','fontsize',14);
            ylabel('neuron','fontsize',14);
        end
    end
    methods(Hidden = true)
        function inputmat = CreateAffinity(obj)
            x = 1:obj.ncells;
            y = normpdf(x,round(obj.ncells / 2),obj.A);
            y = y / max(y);
            inputmat = SetAffinity(obj,y');
        end
        function Input = SetAffinity(obj,Input)
            % This function rotates the odor affinity of cells
            if obj.odorant ~= round(obj.ncells / 2)
                Input = circshift(Input,obj.odorant - round(obj.ncells / 2));
            end
        end
    end
    methods(Hidden = true,Static = true)
        function Miraster = CreateRefPeriod(Miraster,T,time)
            for ii = 1:T
                I = Miraster(:,ii) == 1;
                if ii + time <= T
                    Miraster(I,ii + 1:ii + time) = 0;
                else
                    Miraster(I,ii + 1:end) = 0;
                end
            end
        end
        function C = CalculateCoherence(NS,sizeblock)
            C = 0;
            for ii = 1:size(NS,1) - 1
                for jj  = ii + 1:size(NS,1)
                    sv1 = sum(NS(ii,:));
                    sv2 = sum(NS(jj,:));
                    if sv1 <= sv2
                        mins = find(NS(ii,:));
                        maxs = find(NS(jj,:));
                    else
                        maxs = find(NS(ii,:));
                        mins = find(NS(jj,:));
                    end
                    if ~isempty(mins)
                        sumcoh = 0;
                        for kk = 1:length(mins)
                            sumcoh = sumcoh + sum(abs(maxs - mins(kk)) <= sizeblock);
                        end
                        C = C + sumcoh;
                    end
                end
            end
            C = C / sum(sum(NS));
        end
    end
end
        