%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classSpkNeuron has some general properties methods for spiking neurons
%
% Licurgo de Almeida
% 04/15/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef classSpkNeuron < classNetGeneral
    properties
        SpikeV = 65e-3; % spike voltage for plot figures
        Vhyper = -10e-3; % hyperpolarization potential
        TFire % time since last spike
        Trefrac = 2; % refractory period
        InitalAct = 0; % activity before simulation. This activity is
        % proportional to the frequency at a given time
        TauAct = 100; % ms.
        ActOverT; % vector with the calculated activity over time
        LearnUnlearn = true; % if true, network have nearning and unlearning
        % if false, network has only learning.
    end
    methods
        function obj = classSpkNeuron(tsim,ncells)
            if nargin == 0
                input = {};
            else
                input = {tsim,ncells};
            end
            obj = obj@classNetGeneral(input{:});
            obj.TFire = ones(obj.ncells,1) * 10000000000;
            obj.ActOverT = ones(1,obj.tsim / obj.dt);
        end
        function Raster(obj,PlotActOverTime)
            % This function plots a rasterplot of the neurons
            % The input parameter is:
            % PlotActOverTime: if true, add a second plot to the figure 
            % where we count the number of spikes in a given time window.
            if nargin < 2
                PlotActOverTime = false;
            end
            scrsz = get(0,'ScreenSize');
            figH = figure;
            set(figH,'position',[0,400,scrsz(3)-0.4*scrsz(3),scrsz(4)-0.6*scrsz(4)]);
            if PlotActOverTime == true
                window = 4; %ms
                x = 1:window:obj.tsim;
                y = SumSpikes(obj,window);
                subplot(2,1,2);
                cla;
                hold on;
                box on;
                bar(x,y);
                axis([0,obj.tsim + 1,0,max(y) + 2]);
                xlabel('time (ms)','fontsize',12);
                ylabel('spikes','fontsize',12);
                subplot(2,1,1);
            end
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
            if PlotActOverTime == false
                xlabel('time (ms)','fontsize',14);
            end
            ylabel('neuron','fontsize',14);
        end
        function ActOverTime(obj)
            % This function plots the activity over time.
            v = obj.ActOverT;
            x = (obj.dt:obj.dt:obj.tsim);
            scrsz = get(0,'ScreenSize');
            figH = figure;
            set(figH,'position',[0,400,scrsz(3)-0.4*scrsz(3),scrsz(4)-0.6*scrsz(4)]);
            title(obj.CellName,'fontsize',16);
            plot(x,v,'k');
            axis([0,obj.tsim + obj.dt,0,max(v) + max(v) * 0.1]);
            xlabel('time (ms)','fontsize',14);
            ylabel('activity','fontsize',14);
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
            A = sum(NS,2);
            A = A / size(NS,2);
            for ii = 1:length(A)
                randNS(ii,:) = randNS(ii,:) <= A(ii);
            end
            Crand = obj.CalculateCoherence(randNS,sizeblock);
            C = 1 - (Crand / Creal);
        end
        function STS = STSync(obj)
            % This function calculates the pike train synchrony (STS) index
            % (Brunel and Wang, 2003), suggested in (Li and Cleland, 2013).
            window = 5; %ms
            y = SumSpikes(obj,window);
            [C,lags] = xcov(y);
            I = find(lags == 0,1);
            STS = C(I) / mean(Frequency(obj))^2;
        end
        function F = Frequency(obj)
            % This function extracts the average frequency of a group of
            % cells
            F = sum(obj.O,2) / (obj.tsim / 1000);
        end
        function W = ChangeWeights(obj,M,W,Cpre)
            % This function changes the weight matrix between cells
            % Input parameters are:
            % M: connection matrix between two groups os cells
            % W: tweight matrix between two groups os cells
            % Cpre: presynaptic network object
            [ipost,bglu] = WeightCurves(obj,Cpre);
            Tlearning = size(ipost,2);
            W = CalculateWeights(obj,M,W,ipost,bglu,Tlearning);
        end
    end
    methods(Hidden = true)
        function y = SumSpikes(obj,window)
            Scount = zeros(obj.ncells,ceil(obj.tsim / window));
                for ii = 1:obj.ncells
                    c = (1:round(window / obj.dt));
                    for jj = 1:size(Scount,2)
                        if sum(obj.O(ii,c)) >=1
                            Scount(ii,jj) = 1;
                        end
                        c = c + round(window / obj.dt);
                    end
                end
                y = sum(Scount);
        end
        function [ipost,bglu] = WeightCurves(obj,Cpre)
            % This functions creates the curves for the synaptic weight
            % changes
            
            % Create variables to reduce overhead
            Spost = obj.O;
            Spre = Cpre.O;
            dt = obj.dt;
            Taur = obj.Taur;
            Tauf = obj.Tauf;
            Taupost = obj.Taupost;
            t0post = zeros(size(Spost,1),1);
            t0post = t0post - 10000000000;
            t0pre = t0post;
            
            ipost = zeros(size(Spost));
            bglu = ipost;
            
            for ii = 1:size(Spost,2)
                t = ii * dt;
                I = Spost(:,ii) == 1;
                t0post(I) = t;
                J = Spre(:,ii) == 1;
                t0pre(J) = t;
                ipost(:,ii) = ((t - t0post) ./ Taupost) .* ...
                    exp(1 - (t - t0post) ./ Taupost);
                bglu(:,ii) = exp(-((t - t0pre) ./ Tauf)) .* ...
                    (1 - exp(-((t - t0pre) ./ Taur)));
            end
        end
        function W = CalculateWeights(obj,M,W,ipost,bglu,Tlearning)
            % This function calculates the changes in the synaptic weights
            % between two cells
            dt = obj.dt;
            Mod = obj.Mod;
            delay = round(obj.Tdelay / obj.dt);
            b = circshift(bglu,[0 -delay]);
            M = M';
            Tau11 = obj.Tau11;
            if obj.LearnUnlearn == true
                i = ipost ./ obj.Tau01;
                bb = b ./ obj.Tau10;
                for ii = 1:size(M,1)
                    auxW = W(ii,:)';
                    for jj = 1:Tlearning - delay
                        auxW = auxW + Mod(jj) * dt * (((ipost(ii,jj) .* b(:,jj)...
                            .* M(:,ii)) ./ Tau11) .* (1 - auxW) + (i(ii,jj)...
                            + bb(:,jj) .* M(:,ii)) .* (0 - auxW));

                    end
                    W(ii,:) = auxW';
                end
            else
                for ii = 1:size(M,1)
                    auxW = W(ii,:)';
                    for jj = 1:Tlearning - delay
                        auxW = auxW + Mod(jj) * dt * (((ipost(ii,jj) .* b(:,jj)...
                            .* M(:,ii)) ./ Tau11) .* (1 - auxW));
                        %                         auxW = auxW + dt * (((ipost(ii,jj) .* b(:,jj)...
                        %                             .* M(:,ii)) ./ Tau11) .* (1 - auxW)); % For NE simulations
                    end
                    W(ii,:) = auxW';
                end
            end
        end
    end
    methods(Hidden = true,Static = true)
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