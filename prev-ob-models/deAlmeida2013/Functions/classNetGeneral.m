%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Superclass classNetGeneral has some general proprieties common to all 
% networks
%
% Licurgo de Almeida
% 04/15/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef classNetGeneral
    properties
        CellName = 'Neuron';
        dt = 0.5; % timestep (in ms)
        tsim = 1000; % simulation time (in ms)
        ncells = 100; % number of neurons in the network
        tinit = 1; % stimulus begin (in ms)
        tfinal = 1000; % stimulus end (in ms)
        tau = 5; % charging time constant of the neuron (ms).
        Tmin = [0e-3,0e-3]; % effects of modulation in the resting 
        % potential, represents the range THETA min in eq. The order here
        % is [Mod OFF, Mod ON];
        Tmax = [15e-3,15e-3]; % effects of modulation in the firing 
        % threshold, represents THETA max in eq. the order here is
        % [Mod OFF, Mod ON];
        O % Output (0 or 1 for spiking neurons); between 0 and 1 for non-
        % spiking neurons
        V % Voltage at the time step t
        ModValue = 0; % Value of modulation (varies between 0 and 1), if set to -1 the
        % changes in Mod are dynamical
        Mod % Modulation over time
        Beta = 1; % non-linearity of output
        K = 3; % ligand concentration producing half occupation (Hill eq.)
        % This parameter is very sensitive to the overall activation of the
        % overall activation of the HDB. See the funtcion CalculateAffinity.m
        % for details.
        Cmax = 10; % concentration range of the ligant (Hill eq.). See the
        % funtcion CalculateAffinity.m for details.
        b = 1; % non-linearity of the modulation See the funtcion
        % CalculateAffinity.m for details.
    end
    methods
         function obj = classNetGeneral(tsim,ncells)
            if(nargin > 0)
                obj.tsim = tsim;
                obj.ncells = ncells;
                obj.tfinal = tsim;
                obj.O = zeros(ncells,tsim / obj.dt);
                obj.V = zeros(ncells,tsim / obj.dt);
                obj.Mod = zeros(1,tsim / obj.dt);
            else
                obj.O = zeros(obj.ncells,obj.tsim / obj.dt);
                obj.V = zeros(obj.ncells,obj.tsim / obj.dt);
                obj.Mod = zeros(1,tsim / obj.dt);
            end
         end
        function M = SetConnections(obj,Ninput,Cchance,Conntype)
            % This function creates the connection matrix between two
            % groups of neurons.
            % Ninput the number of cells in the input network
            % Ninput the number of cells in the input network
            % Cchance is the chance of connections between networks
            % Conntype is the type of connection. This can be 'normal' for
            % connections between two different networks or 'auto' for
            % autoassociative connections. When Conntype == 'auto', neurons
            % cannot be connected to each other.
            if nargin == 3
                Conntype = 'normal';
            end
            M = zeros(obj.ncells,Ninput);
            switch lower(Conntype)
                case 'normal'
                    for ii = 1:size(M,1)
                        r = randperm(Ninput);
                        r = r(1:round(Ninput * Cchance));
                        M(ii,r) = 1;
                    end
                case 'auto' % if Conntype is auto, one extra step is added
                    % to make sure that neurons are not connected to itself
                    for ii = 1:size(M,1)
                        r = randperm(Ninput);
                        r = r(1:round(Ninput * Cchance) + 1);
                        M(ii,r) = 1;
                        if M(ii,ii) == 1 % test if neuron is connected to
                            % itself
                            M(ii,ii) = 0;
                        else
                            M(ii,r(end)) = 0;
                        end
                    end
                case 'grouped' % if Conntype is grouped, the connections is
                    % made from neurons close to each other. This method should
                    % only be used in cases where
                    for ii = 1:size(M,1)
                        n = round((Ninput * Cchance) / 2); % divide the number
                        % of connections in two groups
                        i = ii * round(Ninput / size(M,1));
                        nhigh = (i:i + n);
                        nlow = (i:-1:i - n);
                        I = nhigh > Ninput;
                        nhigh(I) = nhigh(I) - Ninput;
                        I = nlow < 1;
                        nlow(I) = nlow(I) + Ninput;
                        M(ii,[nlow,nhigh]) = 1;
                    end
            end
        end
        function S = Sparseness(obj)
            Vact = sum(obj.O,2);
            N = obj.ncells;
            S = (1 - ((sum(Vact / N)^2) / sum((Vact.^2) / N))) / (1 - 1 / N);
        end
        function PlotVoltage(obj,neuron)
            for ii = 1:length(neuron)
                % This function plots a rasterplot of the neurons
                scrsz = get(0,'ScreenSize');
                figH = figure;
                set(figH,'position',[0,400,scrsz(3)-0.4*scrsz(3),scrsz(4)-0.6*scrsz(4)]);
                cla;
                hold on;
                box on;
                title([obj.CellName ', neuron ' num2str(neuron(ii))],'fontsize',16);
                x = (obj.dt:obj.dt:obj.tsim); %time
                plot(x,obj.V(neuron(ii),:),'k');
                ylabel('Voltage (V)','FontSize',12);
                xlabel('time (ms)','FontSize',12);
                axis([0,x(end) + obj.dt,-20e-3,70e-3]);
            end
        end
    end
end