function [Mitral GraProximal GraDistal param InputCurrent] = OB_network_GCE(inputFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boleslaw Osinski 2015
%
% This is a modification of code originally developed by Licurgo de
% Almeida 2013
%
% Publication describing the model...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
% inputFile - txt file with parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS:
% Mitral, GraProximal, GraDistal - structures containing simulated neurons
% param      -  model parameters
% InputCurrent - structure containing all the input currents
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting program

tic;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Reading parameters from input file and creating neurons
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param.inputFile = inputFile;
    
    % Open the input file for reading
    fid1 = fopen(param.inputFile,'r');
    if fid1 == -1
        msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
        return;
    end
    
    str = fgetl(fid1);
    
    % Get network parameters.
    while ~strcmpi(str,'neurons')
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % network parameters
            otherwise
                param = SetNetworkParameters(param,str);
        end  
        str = fgetl(fid1);
    end
    
    
    % Create cells
    [param,Mitral,GraProximal,GraDistal] = CreateCells(param);

    str = fgetl(fid1);
    
    % Get cell parameters.
    while ~strcmpi(str,'end')
        
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % neuronal parameters
            case 'mitral'
                celltype = 'mitral';
            case 'graproximal'
                celltype = 'graproximal';
            case 'gradistal'
                celltype = 'gradistal';
            otherwise
                switch celltype
                    case 'mitral'
                        Mitral = SetNeuronParameters(Mitral,param.nMitral,str);
                    case 'graproximal'
                        GraProximal = SetNeuronParameters(GraProximal,param.nGraprox,str);
                    case 'gradistal'
                        GraDistal = SetNeuronParameters(GraDistal,param.nGradist,str);    
                end
                
        end
        
        str = fgetl(fid1);
    end

    
    fclose(fid1); % Close input file
    fname = inputFile(1:end - 3);
    fname = strcat(fname,'mat');
    save(fname,'Mitral','GraProximal','GraDistal','param');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neuronal activity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Mitral GraProximal GraDistal param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,param);

toc;
end

function param = SetNetworkParameters(param,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----------OB--------------OB----------------OB--------------------------
%
% This function sets the parameters for the different neurons
%
% Modified by Boleslaw Osinski on 06/14/2013
%
% Licurgo de Almeida
% 12/20/2010
%
% Information not related with the parameters of different neuros.
% Path: path where we save input and output file
% dt: timestep (in ms)
% tsim: simulation time (in ms)
% tinit: stimulus begin (in ms)
% tfinal: stimulus end (in ms)
% nMitral: number of mitral cells
% nGradist: number of granule distal synapses
% nGraprox: number of granule cell soma
% GraGracon = if true, granule cells connect to each other
% DistalON = if true graded inhibitory distal Granule dendrites are present
% ProximalON = if true spiking inhibitory proximal Granule soma are present
% BulbH = Bulb height (in distance units)
% BulbW = Bulb width (in distance units)
% NoiseMit = Mitral cell noise std
% NoiseGraprox = Granule proximal dendrite noise std
% NoiseGradist = Granule distal dendrite noise std
% hCaflag = if true, h (VDCC innactiavtion) depends on [Ca]
% rhoCa = proportionality constant between [Ca] and ICa
% ExFrac = fraction of excited dGCs
% Respiration = if true, the input is modulated by an oscillation
% representing the respiration
% RespFreq = Respiratory frequency
% Inoise = noise fraction of ORN input
% Wmin = minimum of MC input weight
% SpikeV = Spike voltage
% CChanceGraMit = Chance of connection between Gra and Mit
% CChanceGraGra = Chance of connection between Gra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    
    % OB parameters
    case 'path'
        param.outputPath = ParValue; %path
    case 'dt'
        param.dt = str2double(ParValue); %step
    case 'tsim'
        param.tsim = str2double(ParValue); %simulation time
    case 'tinit'
        param.tinit = str2double(ParValue); %stimulus begin
    case 'tfinal'
        param.tfinal = str2double(ParValue); %stimulus end
    case 'nmitral'
        param.nMitral = str2double(ParValue); %number of Mitral cells
    case 'ngradist'
        param.nGradist = str2double(ParValue); %number of granule distal synapses
    case 'ngraprox'
        param.nGraprox = str2double(ParValue); %number of Granule cell soma
    case 'gragracon'
        param.GraGracon = str2double(ParValue); %Gra-Gra connections
    case 'distalon'
        param.DistalON = str2num(ParValue); %flag distal gra dendrites
    case 'proximalon'
        param.ProximalON = str2num(ParValue); %flag proximal gra dendrites
    case 'bulbh'
        param.BulbH = str2double(ParValue); %Bulb height
    case 'bulbw'
        param.BulbW = str2double(ParValue); %Bulb width
    case 'noisemit'
        param.noisemit = str2double(ParValue); %noise Mitral
    case 'noisegraprox'
        param.noisegraprox = str2double(ParValue); %noise granule prox
    case 'noisegradist'
        param.noisegradist = str2double(ParValue); %noise granule dist
    case 'preset'
        param.flagpreset = str2num(ParValue); %flag preset positions
    case 'hcaflag'
        param.hCaflag = str2num(ParValue); %flag [Ca] dependence of h   
    case 'rhoca'
        param.rhoCa = str2num(ParValue); %proportionality factor between ICa and [Ca]    
    case 'exfrac'
        param.ExFrac = str2num(ParValue); %fraction of excited dGCs     
    case 'respiration'
        param.flagRespiration = str2num(ParValue); %flag respiratory modulation
    case 'respfreq'
        param.RespFreq = str2double(ParValue); %respiratory frequency
    case 'inoise'
        param.Inoise = str2double(ParValue); %noise fraction of ORN input
    case 'wmin'
        param.Wmin = str2double(ParValue); %minimum of MC EXT input weight
    case 'spikev'
        param.SpikeV = str2double(ParValue); %Spike voltage
    case 'cchancegramit'
        param.CChanceGraMit = str2double(ParValue); %chance of connection between Granule and Mitral cells
    case 'cchancegragra'
        param.CChanceGraGra = str2double(ParValue); %chance of connection between Granule cells
            
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end


function [param,Mitral,GraProximal,GraDistal] = CreateCells(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function initiallizes the structure for each neuron type
%
% Modified by Boleslaw Osinski on 06/16/2013 and 08/13/2014
%
% Licurgo de Almeida
% 12/21/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    Mitral = cell(param.nMitral,1);
    for ii = 1:param.nMitral
        Mitral{ii}.input = []; %no input for now
        Mitral{ii}.label = 'Mitral';
    end
    
    GraProximal = cell(param.nGraprox,1);
    GraDistal = cell(param.nGradist,1);
    for ii = 1:param.nGraprox
        GraProximal{ii}.input = [];
        GraProximal{ii}.label = 'GraProximal';
    end
    for ii = 1:param.nGradist
        GraDistal{ii}.input = [];
        GraDistal{ii}.label = 'GraDistal';
    end


end


function N = SetNeuronParameters(N,ncells,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
% Modified by Boleszek Osinski on 06/18/2013
%
%
% * tau = charging time constant of the neuron (ms). ms (not s) is the basic time unit in this program.
% * Fthresh = Firing threshold in V
% * CCaTh = Threshold [Ca] for maximum GABA release
% * Vrest = resting potential
% * Vhyper = hyperpolarization potential
% * EAMPA = AMPA's Nernst potential.
% * ECa = Ca Nernst potential.
% * EGABA = GABA's Nernst potential.
% * tauAMPA1 = AMPA's rising tau.
% * tauAMPA2 = AMPA's falling tau.
% * tauGABA1 = GABA's rising tau.
% * tauGABA2 = GABA's falling tau.
%
%
% * gmaxAMPA = AMPA's max conductance
% * gmaxGABA = GABA's max conductance
% * gmaxGABAP = Periglomerula GABA's max conductance (for Mitral cells
% only)
% * IACh = Addition current when ACh is ON in non-spiking cells.
% * nGracon = number of granule cells connected to mitral or granule cells
% * tauPROX1 = proximal Mit-Gra synapse rising tau (if distal dendrites removed)
% * tauPROX2 = proximal Mit-Gra synapse falling tau (if distal dendrites removed)
% * tauAMPA1 = distal Mit-Gra AMPA receptor rising tau (Gra cell only)
% * tauAMPA2 = distal Mit-Gra AMPA receptor falling tau (Gra cell only)
% * tauNMDA1 = distal Mit-Gra NMDA receptor rising tau (Gra cell only)
% * tauNMDA2 = distal Mit-Gra NMDA receptor falling tau (Gra cell only)
% * tauVDCC = time constant of VDCC activation variable
% * wAMPAMI = excitatory synaptic weight from Mitral cell to Granule cell AMPA
% * wNMDAMI = excitatory synaptic weight from Mitral cell to Granule cell NMDA
% * wVDCCMI = synaptic conductance of Granule cell VDCC
% * wGABAGR = inhibitory synaptic weight from Granule cell to Gra or Mit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    % OB and PC
    case 'tau'
        for ii = 1:ncells
            N{ii}.tau = str2double(ParValue); %time in ms
        end
    case 'fthresh'
        for ii = 1:ncells
            N{ii}.FThresh = str2double(ParValue); %threshold in V
        end        
    case 'ccath'
        for ii = 1:ncells
            N{ii}.CCaTh = str2double(ParValue); %[Ca] in uM
        end
    case 'vrest'
        for ii = 1:ncells
            N{ii}.Vrest = str2double(ParValue); %potential in volts
        end
    case 'vhyper'
        for ii = 1:ncells
            N{ii}.Vhyper = str2double(ParValue); %potential in volts
        end
    case 'noise'
        for ii = 1:ncells
            N{ii}.Noise = str2double(ParValue);
        end
    case 'eampa'
        for ii = 1:ncells
            N{ii}.EAMPA = str2double(ParValue); %AMPA (Na) reversal potential in volts
        end
    case 'eca'
        for ii = 1:ncells
            N{ii}.ECa = str2double(ParValue); %Ca reversal potential in volts
        end
    case 'egaba'
        for ii = 1:ncells
            N{ii}.EGABA = str2double(ParValue); %potential in volts
        end
    case 'tauampa1'
        for ii = 1:ncells
            N{ii}.tauAMPA1 = str2double(ParValue); %time in ms
        end
    case 'tauampa2'
        for ii = 1:ncells
            N{ii}.tauAMPA2 = str2double(ParValue); %time in ms
        end
    case 'taugaba1'
        for ii = 1:ncells
            N{ii}.tauGABA1 = str2double(ParValue); %time in ms
        end
    case 'taugaba2'
        for ii = 1:ncells
            N{ii}.tauGABA2 = str2double(ParValue); %time in ms
        end
        
        % OB
    case 'gmaxampa'
        for ii = 1:ncells
            N{ii}.gmaxAMPA = str2double(ParValue); %AMPA channel
            % conductance in siemens
        end
    case 'gmaxgaba'
        for ii = 1:ncells
            N{ii}.gmaxGABA = str2double(ParValue); %GABA channel
            % conductance in siemens
        end
    case 'ngracon'
        for ii = 1:ncells
            N{ii}.NGraCon = str2double(ParValue); %number of granule cells
            % connected to a mitral or granule cell
        end
    case 'tauprox1'
        for ii = 1:ncells
            N{ii}.tauPROX1 = str2double(ParValue); %time in ms
        end
    case 'tauprox2'
        for ii = 1:ncells
            N{ii}.tauPROX2 = str2double(ParValue); %time in ms
        end
    case 'taudist1'
        for ii = 1:ncells
            N{ii}.tauAMPA1 = str2double(ParValue); %time in ms
        end
    case 'taudist2'
        for ii = 1:ncells
            N{ii}.tauAMPA2 = str2double(ParValue); %time in ms
        end
    case 'taunmda1'
        for ii = 1:ncells
            N{ii}.tauNMDA1 = str2double(ParValue); %time in ms
        end
    case 'taunmda2'
        for ii = 1:ncells
            N{ii}.tauNMDA2 = str2double(ParValue); %time in ms
        end
    case 'tauvdcc'
        for ii = 1:ncells
            N{ii}.tauVDCC = str2double(ParValue); %time in ms
        end        
    case 'wampami'
        for ii = 1:ncells
            N{ii}.wAMPAMI = str2double(ParValue); % excitatory synaptic weight from Mitral to Granule AMPA
        end
    case 'wnmdami'
        for ii = 1:ncells
            N{ii}.wNMDAMI = str2double(ParValue); % excitatory synaptic weight from Mitral to Granule NMDA
        end
     case 'wvdccmi'
        for ii = 1:ncells
            N{ii}.wVDCCMI = str2double(ParValue); % VDCC synaptic conductance
        end
    case 'wgabagr'
        for ii = 1:ncells
            N{ii}.wGABAGR = str2double(ParValue); % inhibitory synaptic weight from Granule to Gra or Mit
        end
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mitral GraProximal GraDistal param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates the neuronal activity
%
% Licurgo de Almeida
% 11/03/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create initial setup -- OB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if param.flagRespiration == true
    Respiration = CreateRespFreq(param);
else
    Respiration = ones(1,round(param.tsim / param.dt));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%   OB   %%%%%

nds = param.nGradist/param.nGraprox; % number of distal synapses per granule cell


    % Each MC connected to random subset CChanceGraMit*nGradist dGCs
    % In this scheme each MC is connected to exactly CChanceGraMit*nGradist
    % dGCs, but each dGC can be connected to a range of MCs
        MatGradistMit = SetConnections(param.nMitral,param.nGradist,param.CChanceGraMit);
        
        MatMitGradist = MatGradistMit'; % reciprocal synapses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set weights matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%   OB   %%%%%

    wGABAMit = zeros(param.nMitral,1);
    for ii = 1:param.nMitral
        wGABAMit(ii) = Mitral{ii}.wGABAGR;
    end
    wGradistMit = SetWeights(MatGradistMit,wGABAMit);
    
    wAMPAGradist = zeros(param.nGradist,1);
    wNMDAGradist = zeros(param.nGradist,1);
    wVDCCGradist = zeros(param.nGradist,1);
    wGABAGraProx = zeros(param.nGraprox,1);
    for ii = 1:param.nGradist
        wAMPAGradist(ii) = GraDistal{ii}.wAMPAMI;
        wNMDAGradist(ii) = GraDistal{ii}.wNMDAMI;
        wVDCCGradist(ii) = GraDistal{ii}.wVDCCMI;
    end
    for ii = 1:param.nGraprox
        wGABAGraProx(ii) = GraProximal{ii}.wGABAGR;
    end
    wMitGradistAMPA = SetWeights(MatMitGradist,wAMPAGradist);
    wMitGradistNMDA = SetWeights(MatMitGradist,wNMDAGradist);
    wMitGradistVDCC = SetWeights(MatMitGradist,wVDCCGradist);
%     wProxProx = SetWeights(MatProxProx,wGABAGraProx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Mitral cell at a given time
Vmit = zeros(param.nMitral,round(param.tsim / param.dt));
Smit = Vmit; % Binary matrix recording only the spikes
Iext_matrix = Vmit; % Initialize extrenal inputs to MCs

refracmit = 3; % Refractory period after firing (ms)

Restmit = zeros(param.nMitral,1); % resting potential
Threshmit = Restmit; % current firing threshold
Hypermit = Restmit; % hyperpolarization potential
taumit = Restmit; % tau neuron
countrefracmit = Restmit; % Refractory period counter. This is assumed to
% be the same for all Mitral cells so we don't need to add it to the for
% below.
gmaxAMPAmit = Restmit; % Max AMPA conductance
tauAMPA1mit = Restmit; % AMPA's rising tau.
tauAMPA2mit = Restmit; % AMPA's falling tau.
EAMPAmit = Restmit; % AMPA's Nernst potential.
gmaxGABAmit = Restmit; % Max GABA conductance
tauGABA1mit = Restmit; % GABA's rising tau.
tauGABA2mit = Restmit; % GABA's falling tau.
EGABAmit = Restmit; % GABA's Nernst potential.


for ii = 1:param.nMitral
    Restmit(ii) = Mitral{ii}.Vrest;
    Hypermit(ii) = Mitral{ii}.Vhyper;
    taumit(ii) = Mitral{ii}.tau;
    gmaxAMPAmit(ii) = Mitral{ii}.gmaxAMPA;
    tauAMPA1mit(ii) = Mitral{ii}.tauAMPA1;
    tauAMPA2mit(ii) = Mitral{ii}.tauAMPA2;
    EAMPAmit(ii) = Mitral{ii}.EAMPA;
    gmaxGABAmit(ii) = Mitral{ii}.gmaxGABA;
    tauGABA1mit(ii) = Mitral{ii}.tauGABA1;
    tauGABA2mit(ii) = Mitral{ii}.tauGABA2;
    EGABAmit(ii) = Mitral{ii}.EGABA;

    Threshmit(ii) = Mitral{ii}.FThresh; 
    
    Mitral{ii}.Connections = MatGradistMit(ii,:);
    Mitral{ii}.ConWeights = wGradistMit(ii,:);
end

% Initialize Mitral cells potentials
Vmit(:,1) = Restmit;
Vmit_nospike = Vmit;


% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Only used if paramal.ProxON is
% true and param.DistalON is false (i.e. spiking inhibition)
tGABA0mit = zeros(param.nGradist,round(param.tsim / param.dt)) - 10000000;

Igradistmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igradistmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

% only used when distal gra dendrites are removed
Igraproxmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igraproxmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

maxgGABAmit = getmaxg(param.dt,tauGABA1mit,tauGABA2mit); % Get max conductance
% amplitude


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set Proximal Granule cell parameters and variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if param.ProximalON == true
%     
% % Stores the voltage of each Granule cell at a given time
% Vgraprox = zeros(param.nGraprox,round(param.tsim / param.dt));
% % Binary matrix recording only the spikes
% Sgraprox = Vgraprox;
% % Probability of each cell firing
% Pfiregraprox = Vgraprox;
% 
% refracgraprox = 5; % Refractory period after firing (ms)
% 
% Restgraprox = zeros(param.nGraprox,1); % resting potential
% Threshgraprox = Restgraprox; % current firing threshold
% Hypergraprox = Restgraprox; % hyperpolarization potential
% taugraprox = Restgraprox; % tau neuron
% countrefracgraprox = Restgraprox; % Refractory period counter. This is assumed to
% % be the same for all Granule cells so we don't need to add it to the for below.
% 
% gmaxAMPAgraprox = Restgraprox; % Max AMPA conductance
% 
% tauPROX1 = Restgraprox; % AMPA's rising tau.
% tauPROX2 = Restgraprox; % AMPA's falling tau.
% 
% EAMPAgraprox = Restgraprox; % AMPA's Nernst potential.
% gmaxGABAgraprox = Restgraprox; % Max GABA conductance
% tauGABA1graprox = Restgraprox; % GABA's rising tau.
% tauGABA2graprox = Restgraprox; % GABA's falling tau.
% EGABAgraprox = Restgraprox; % GABA's Nernst potential.
% 
% 
% for ii = 1:param.nGraprox
%     Restgraprox(ii) = GraProximal{ii}.Vrest;
%     Hypergraprox(ii) = GraProximal{ii}.Vhyper;
%     taugraprox(ii) = GraProximal{ii}.tau;
%     gmaxAMPAgraprox(ii) = GraProximal{ii}.gmaxAMPA;
%     tauPROX1(ii) = GraProximal{ii}.tauPROX1;
%     tauPROX2(ii) = GraProximal{ii}.tauPROX2;
%     EAMPAgraprox(ii) = GraProximal{ii}.EAMPA;
%     gmaxGABAgraprox(ii) = GraProximal{ii}.gmaxGABA;
%     tauGABA1graprox(ii) = GraProximal{ii}.tauGABA1;
%     tauGABA2graprox(ii) = GraProximal{ii}.tauGABA2;
%     EGABAgraprox(ii) = GraProximal{ii}.EGABA;
%     
%     Threshgraprox(ii) = GraProximal{ii}.CCaTh;
%     
% end
% 
% % Initialize Granule cells potentials
% Vgraprox(:,1) = Restgraprox;
% Vgraprox_nospike = Vgraprox;
% 
% % bulbar input current to graprox
% Idistprox = zeros(param.nGraprox,1); % Input coming from dist to prox
% Idistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% Imitgraprox = zeros(param.nGraprox,1);
% Imitgraprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% 
% % input voltage (for direct distprox summing)
% Vdistprox = zeros(param.nGraprox,1); % V input coming from dist to prox
% Vdistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% 
% 
% % AMPA time counter. This variable starts with a very negative value just
% % to make sure that the currents will be = 0
% 
% % tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
% tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
% 
% 
% maxgAMPAgraprox = getmaxg(param.dt,tauPROX1,tauPROX2); % Get max conductance amplitude
% 
% % maxgAMPAgraprox = 1; % Set gmax = 1 because it will be normalized anyways
% 
% % GABA time counter. This variable starts with a very negative value just
% % to make sure that the currents will be = 0. Each Granule cell can be
% % connected to a varied number of other granule cells
% tGABA0graprox = zeros(param.nGraprox,round(param.tsim / param.dt)) - 10000000;
% Igragra = zeros(param.nGraprox,1); % Input coming from other granule cells
% Igragra_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% 
% maxgGABAgraprox = getmaxg(param.dt,tauGABA1graprox,tauGABA2graprox); % Get max conductance amplitude
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % End of Proximal Granule cells parameters and variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Distal Granule cell parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.DistalON == true
% Stores the voltage of each Granule cell at a given time
Vgradist = zeros(param.nGradist,round(param.tsim / param.dt));
CCa = zeros(param.nGradist,round(param.tsim / param.dt));% Ca concentration
Sgradist= Vgradist;% Binary matrix recording only the spikes


Prelease_matrix = Vgradist; % Prelease_matrix is the probability of GABA release 
                             % when using graded inhibition from distal gra dendrites

Restgradist = zeros(param.nGradist,1); % resting potential
Threshgradist = Restgradist; % current firing threshold
Hypergradist = Restgradist; % hyperpolarization potential
taugradist = Restgradist; % tau neuron
gmaxAMPAgradist = Restgradist; % Max AMPA conductance

tauAMPA1 = Restgradist; % AMPA's rising tau.
tauAMPA2 = Restgradist; % AMPA's falling tau.
tauNMDA1 = Restgradist; % NMDA's rising tau.
tauNMDA2 = Restgradist; % NMDA's falling tau.
tauVDCC = Restgradist; % VDCC activation tau.

% L-type Ca activation parameters
mbar = Restgradist;
taum = Restgradist;
m = Vgradist;
h = Vgradist;

EAMPAgradist = Restgradist; % AMPA's Nernst potential.
ENMDAgradist = Restgradist; % NMDA's Nernst potential.
EVDCCgradist = Restgradist; % Ca Nernst potential.

for ii = 1:param.nGradist
    Restgradist(ii) = GraDistal{ii}.Vrest;
    Hypergradist(ii) = GraDistal{ii}.Vhyper;
    taugradist(ii) = GraDistal{ii}.tau;
    gmaxAMPAgradist(ii) = GraDistal{ii}.gmaxAMPA;
    tauAMPA1(ii) = GraDistal{ii}.tauAMPA1;
    tauAMPA2(ii) = GraDistal{ii}.tauAMPA2;
    tauNMDA1(ii) = GraDistal{ii}.tauNMDA1;
    tauNMDA2(ii) = GraDistal{ii}.tauNMDA2;
    tauVDCC(ii) = GraDistal{ii}.tauVDCC;
    EAMPAgradist(ii) = GraDistal{ii}.EAMPA;
    ENMDAgradist(ii) = GraDistal{ii}.EAMPA; % nmda reversal potential is the same as AMPA
    EVDCCgradist(ii) = GraDistal{ii}.ECa; % Ca reversal potential is ~100mV higher than AMPA
    % EVDCCgradist(ii) = GraDistal{ii}.EAMPA;
    
    Threshgradist(ii) = GraDistal{ii}.CCaTh;
    
    GraDistal{ii}.Connections = MatMitGradist(ii,:);
end
% 
% 1-ExFrac subpopulation of GCs in unexcited state
if param.ExFrac < 1
    Restgradist((round(param.ExFrac*param.nGradist)+1):end) = -75e-3;
end

% Initialize Granule cells potentials
Vgradist(:,1) = Restgradist;
CCa(:,1) = 0.1; % [uM]

% proximal granule input currents to gradist
Iproxdist = zeros(param.nGradist,1); % Input coming from prox to dist
Iproxdist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% mitral input to gradist AMPA and NMDA receptors
ImitgradistAMPA = zeros(param.nGradist,1);
ImitgradistAMPA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

ImitgradistNMDA = zeros(param.nGradist,1);
ImitgradistNMDA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% L-type Ca current
ImitgradistVDCC = zeros(param.nGradist,1);
ImitgradistVDCC_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
% Wmitgradist = ones(param.nGranule,param.nMitral); % Synaptic weights

maxgAMPAgradist = getmaxg(param.dt,tauAMPA1,tauAMPA2); % Get max conductance amplitude for AMPA
maxgNMDAgradist = max(exp(-([0:.01:200])./tauNMDA2(1)) - exp(-([0:.01:200])./tauNMDA1(1))); % Get max conductance amplitude for NMDA


% Calculate baseline [Ca] from steady state solution of [Ca]' = 0
mbarrest = 1./(1 + exp(-(1e3*Restgradist+45)./7)); % Restgradist must be in units of mV
RT = 300*8.31; % J/mole
z= 2 ; % Ca ion valence
F = 96485; % Faraday constant Coul/mole
Cout = 1500;
Cin = 0.002:0.002:2.5;
EVDCCgradist = (RT/(z*F))*log(Cout./Cin); % RHS of equation
LHS = zeros(param.nGradist,length(Cin));
for ii = 1:length(Cin)
    for nn = 1:param.nGradist
        LHS(nn,ii) = Restgradist(nn) + (Cin(ii)^2)./(1e-4 .* param.rhoCa .* sum(wMitGradistVDCC(nn,:)) .* mbarrest(nn));
    end
end

% plot LHS and RHS
% subplot(2,1,1)
% plot(Cin,LHS(1,:),'k.',Cin,EVDCCgradist,'k--')
% set(gca,'fontsize',16)
% legend('LHS','RHS','location','best')
% % plot |LHS - RHS|
% subplot(2,1,2)
% plot(Cin,abs(LHS(1,:) - EVDCCgradist),'k')
% set(gca,'fontsize',16)
% title('abs(LHS - RHS)')
% xlabel('[Ca] (\muM)')

% old solution (without explicit [Ca] dpendence in Eca)
% CCaBase = sqrt(1e-4 * nmit_per_gra .* param.rhoCa .* wVDCCGradist .* mbarrest .* (EVDCCgradist - Restgradist));
% Note: the 1e-4 comes from the definition of h for N-type current:
%    h = 1e-4/(1e-4 + [Ca])

% new solution (with explicit [Ca] dpendence in Eca)
CCaBase = zeros(param.nGradist,1);
if wVDCCGradist > 0
    for ii = 1:param.nGradist
    % mind = find(abs(LHS(ii,:) - EVDCCgradist) == min(abs(LHS(ii,:) - EVDCCgradist)));
    CCaBase(ii) = Cin(abs(LHS(ii,:) - EVDCCgradist) == min(abs(LHS(ii,:) - EVDCCgradist)));
    end
end

if sum(wVDCCGradist) > 0
    EVDCCgradistBase = (RT/(z*F))*log(Cout./CCaBase);
    WMGVDCCsum = zeros(param.nGradist,1);
    for ii = 1:param.nGradist
        WMGVDCCsum(ii) = sum(wMitGradistVDCC(ii,:));
    end
    IVDCCBase = 1e-4*WMGVDCCsum.*mbarrest.*(EVDCCgradistBase - Restgradist)./(1e-4 + CCaBase);
else
    IVDCCBase = zeros(param.nGradist,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Distal Granule cell synapse parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin neuron simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iext_max = SetExtInput(param); % set level of external input
% sc = [rand(1,1000).*tanh(linspace(0,pi,1000)) ones(1,round(param.tsim / param.dt)-1000)];
% scaling function designed to soften initial impact

% Start loop
for tt = 2:round(param.tsim / param.dt)
    t = tt * param.dt; % current time
    

    % Mitral Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Get Glo input to Mitral cells. This current is applied directly in
    % Vmit, no channel conductance.

    % External inupt to Mit cells is modulated by respiration if Resp is on
    Iext = Respiration(tt) .* Iext_max(:,tt);
%     Iext = sc(tt) * Respiration(tt) .* Iext_max(:,tt);
    Iext_matrix(:,tt) = Iext(:);
    
            
    % Get Granule graded inputs to Mitral cells
    Igradistmit(:) = SetInoSpike_GraMit(gmaxGABAmit,Prelease_matrix(:,tt-1),EGABAmit,Vmit_nospike(:,tt - 1),...
        wGradistMit);
    Igradistmit_matrix(:,tt) = Igradistmit(:);
    
    Vnoisemit = param.noisemit .* randn(param.nMitral,1);
    
    % Mitral cell potential
    
    % Forwards Euler
    Vmit(:,tt) = Vmit(:,tt - 1) + (param.dt ./ taumit(:)) .* ...
        ((Iext(:) + Igradistmit(:)) - Vmit(:,tt - 1) + Restmit(:) + Vnoisemit);
    
    % Backwards Euler
    
%     g_gramit = (((tauGABA1mit(1) * tauGABA2mit(1)) / (tauGABA1mit(1) - tauGABA2mit(1))) *...
%     (exp(-(t - tGABA0mit(:,tt)) / tauGABA1mit(1)) - exp(-(t - tGABA0mit(:,tt)) / tauGABA2mit(1)))) / maxgGABAmit(1);
% 
%     Vmit(:,tt) = (Vmit(:,tt - 1) + (param.dt ./ taumit(1)) .* ((Iext(:) + (wGraMit_scaled * g_gramit) .* EGABAmit(1) + Restmit(:))))...
%         ./ (1 + (param.dt ./ taumit(1)) .* ((wGraMit_scaled * g_gramit) + 1));
   
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vmit(:,tt - 1) == param.SpikeV;
    Vmit(I,tt) = Hypermit(I);
    Vmit_nospike(:,tt) = Vmit(:,tt);
    
    % I is a vector of 1s or 0s of length ncells
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Granule variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This variable is used for input onto distal gra dendrites
        Mit_spike_time = t;
        tAMPA0mit_gra(I,ceil(Mit_spike_time/param.dt):end) = Mit_spike_time;
        tNMDA0mit_gra = tAMPA0mit_gra;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracmit(I) = refracmit / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracmit > 0;
    countrefracmit(I) = countrefracmit(I) - 1;
    
    I = find(countrefracmit == 0); % if countrefracmit = 0 the neuron can fire
    % again

% spike when V > thresh    
    if ~isempty(I)
    J = find(Vmit(I,tt) >= Threshmit(I));
    if ~isempty(J)
            Vmit(I(J),tt) = param.SpikeV; % Action potential
            Smit(I(J),tt) = 1; % Record spike time
    end
    end
    
    
    % Granule Distal Dendrites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.DistalON == true
        
    % NOTE: In the paper the weights were not included in the definition of the 
    % currents. But we include them in the calculation of the currents because
    % multiplication by weight matrix allows convenient summation of all
    % inputs onto cell. The result is the same.
    
    % Get Mitral input to Distal Granule AMPA receptors
    ImitgradistAMPA(:) = SetI(tauAMPA1(1),tauAMPA2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgradist(1),...
        wMitGradistAMPA,EAMPAgradist(1),Vgradist(:,tt - 1));
    ImitgradistAMPA_matrix(:,tt) = ImitgradistAMPA(:);
    
    % Mitral input to Distal Granule NMDA receptors
    ImitgradistNMDA(:) = SetI_NMDA(tauNMDA1,tauNMDA2,t,tNMDA0mit_gra(:,tt),maxgNMDAgradist,...
        wMitGradistNMDA,ENMDAgradist(1),Vgradist(:,tt - 1));
    ImitgradistNMDA_matrix(:,tt) = ImitgradistNMDA(:);
    
    % Distal Granule VDCC Current
    
            % Calcium Concentration variable
    CCa(:,tt) = CCa(:,tt - 1) + (param.dt ./ taugradist(:)) .* (param.rhoCa .*...
          (ImitgradistNMDA(:) + ImitgradistVDCC(:) )...
          - CCa(:,tt - 1) );
      
%     % L-type
%     % % m' = (mbar - m)/taum;
%    mbar(:) = 1./(1 + exp(-(1e3.*Vgradist(:,tt - 1)+50)./3));
%    taum(:) = 18 * exp(-((1e3.*Vgradist(:,tt - 1)+45)./20).^2) + 1.5;
%    h(:,tt) = 0.00045./(0.00045+CCa(:,tt));

    % Ntype
    mbar(:) = 1./(1 + exp(-(1e3.*Vgradist(:,tt - 1)+45)./7));
    taum(:) = tauVDCC .* exp(-((1e3.*Vgradist(:,tt - 1)+70)./25).^2) + 0.3;
    if param.hCaflag == true
     h(:,tt) = 0.0001./(0.0001+CCa(:,tt));
    else
        h(:,tt) = 0.0001./(0.0001+CCaBase);
        % h(:,tt) = 0.0001./(0.0001+WMGVDCCsum*2.0714e-04);
        % scaling factor to obtain mean [Ca] for WVDCC = 350, WGM = 0.025
        % is 2.0714e-04
    end
     
    mtemp = m(:,tt - 1) + (param.dt ./ taum(:)).*(mbar(:) - m(:,tt - 1));
    mtemp(mtemp < 0) = 0;
    m(:,tt) = mtemp;
    
%     RT = 300*8.31; % J/mole
%     z= 2 ; % Ca ion valence
%     F = 96485; % Faraday constant Coul/mole
%     Cout = 1500;
    Cin = CCa(:,tt);

    EVDCCgradist = (RT/(z*F))*log(Cout./Cin);

    ImitgradistVDCC(:) = SetI_VDCC(wMitGradistVDCC,EVDCCgradist,Vgradist(:,tt - 1),m(:,tt),h(:,tt));
%     ImitgradistVDCC(ImitgradistVDCC<0) = 0;
    ImitgradistVDCC_matrix(:,tt) = ImitgradistVDCC(:);
    
    
    % NOTE!!! All voltages are in V, not mV, so parameters have
    % to be scaled accordingly!!!

    

%     % Get Proximal Granule input to Distal Granule Dendrites
%     if param.ProximalON == true
%         % Prox-dist current is applied directly in Vgradist, no ion channel conductance.
%         % Sgraprox(:,tt-1) is [nGraprox x 1]
%         % MatProxDist is [nGradist x nGraprox]
%         % gmaxAMPAgradist is [nGradist x 1]
%         % we want Iproxdist to be [nGradist x 1] so we need to transpose all 3 terms!
%         Iproxdist(:) = Sgraprox(:,tt-1)'*MatProxDist' .* gmaxAMPAgradist';
%         Iproxdist_matrix(:,tt) = Iproxdist(:);
%     end
%     
    
    Vnoisegradist = param.noisegradist .* randn(param.nGradist,1);
    
    % Distal Granule cell potential
    % Forward Euler
    Vgradist(:,tt) = Vgradist(:,tt - 1) + (param.dt ./ taugradist(:)) .* ...
        ((ImitgradistAMPA(:) + ImitgradistNMDA(:) + (ImitgradistVDCC(:) - IVDCCBase) + Iproxdist(:))...
        - Vgradist(:,tt - 1) + Restgradist(:) + Vnoisegradist);
    
    
    % GABA release probability for Graded inhibition from distal gra dendrites
    Prelease_matrix(:,tt) = Prelease(CCa(:,tt),CCaBase,Threshgradist);
    
end
    
%     % Granule Proximal Soma
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if param.ProximalON == true
%     % Get bulbar input to granule soma
%     
%     
%     if param.DistalON == true
%         % Distal Granule dendrite input to Proximal Granule. This current is applied directly in
%         % Vgraprox, no channel conductance (Nernst potetial not used).
%         % Prelease_matrix is [nGradist x 1]
%         % MatProxDist is [nGradist x nGraprox]
%         % gmaxAMPAgraprox is [nGraprox x 1]
%         % we want Idistprox to be [nGraprox x 1] so we need to transpose only 1st and 3rd terms!
%         Idistprox(:) = Prelease_matrix(:,tt)'*MatProxDist .* gmaxAMPAgraprox';
% %         Vdistprox(:) = (Vgradist(:,tt)'*MatProxDist)./nds; % Avg dendritic membrane potentials per soma
% %         Idistprox(:) = Idistprox(:)./Rgraprox(:);
%         Idistprox_matrix(:,tt) = Idistprox(:);
%     else
%         % if distal dendrites are removed, mit synapse directly onto gra
%         Imitgraprox(:) = SetI(tauPROX1(1),tauPROX2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgraprox(1),...
%         wMitGradistAMPA,EAMPAgraprox(1),Vgraprox_nospike(:,tt - 1));
%         Imitgraprox_matrix(:,tt) = Imitgraprox(:);
%     end
%     
%     
%     
%     % Get Granule input to Granule cells
%     if param.GraGracon == true
%         Igragra(:) = SetI(tauGABA1graprox(1),tauGABA2graprox(1),t,tGABA0graprox(:,tt),maxgGABAgraprox(1),...
%             wProxProx,EGABAgraprox(1),Vgraprox_nospike(:,tt - 1));
%         Igragra_matrix(:,tt) = Igragra(:);
%     end
%     
%     
%     Vnoisegraprox = param.noisegraprox .* randn(param.nGraprox,1);
%     
%     % Proximal Granule cell potential
%     % Forward Euler
%     Vgraprox(:,tt) = Vgraprox(:,tt - 1) + (param.dt ./ taugraprox(:)) .* ...
%         ((Igragra(:) + Idistprox(:) + Imitgraprox(:))...
%         - Vgraprox(:,tt - 1) + Restgraprox(:) + Vnoisegraprox);
% 
%     % Backwards Euler
% %     g_mitgra = (((tauPROX1(1) * tauPROX2(1)) / (tauPROX1(1) - tauPROX2(1))) *...
% %     (exp(-(t - tAMPA0mit_gra(:,tt)) / tauPROX1(1)) - exp(-(t - tAMPA0mit_gra(:,tt)) / tauPROX2(1)))) / maxgAMPAgraprox(1);
% % 
% % 
% %     Vgraprox(:,tt) = (Vgraprox(:,tt - 1) + (param.dt ./ taugraprox(1)) .* (((wMitGradist * g_mitgra) .* EAMPAgraprox(1) + Restgraprox(:))))...
% %         ./ (1 + (param.dt ./ taugraprox(1)) .* ((wMitGradist * g_mitgra) + 1));
% 
%     
%     
%     % If the neuron fired last cycle, neuron potential hyperpotentializes
%     I = Vgraprox(:,tt - 1) == param.SpikeV;
%     Vgraprox(I,tt) = Hypergraprox(I);
%     Vgraprox_nospike(:,tt) = Vgraprox(:,tt);
%     
%     
%     % Proximal Granule Dendrite variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GraGra_spike_time = t; % new t0 for the Distal Granule cell GABA current
%     tGABA0graprox(I,ceil(GraGra_spike_time/param.dt):end) = GraGra_spike_time;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Mitral Cell Variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % NOTE: tGABA0mit has dimensions
%     % [nGradist x dt*tsim]
%         Idist = zeros(param.nGradist,1);
%         for ii = 1:param.nGraprox
%         Idist((nds*(ii-1)+1):nds*ii) = I(ii);
%         end
%         % new t0 for the Mitral cell GABA current
%         Gra_spike_time = t;
%         tGABA0mit(logical(Idist),ceil(Gra_spike_time/param.dt):end) = Gra_spike_time;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     countrefracgraprox(I) = refracgraprox / param.dt; % neurons that just fired get
%     % into the refractory period
%     I = countrefracgraprox > 0;
%     countrefracgraprox(I) = countrefracgraprox(I) - 1;
%     
%     I = find(countrefracgraprox == 0); % if countrefracgra = 0 the neuron can fire
%     % again
%     

% % spike when V > thresh
%     if ~isempty(I)
%         J = find(Vgraprox(I,tt) >= Threshgraprox(I));
%         
%         if ~isempty(J)
%             Vgraprox(I(J),tt) = param.SpikeV; % Action potential
%             Sgraprox(I(J),tt) = 1; % Record spike time
%         end
%     end
%     
% end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ORN and Mitral cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nMitral
    Mitral{ii}.V = Vmit(ii,:); % Save neuronal activity
    Mitral{ii}.S = Smit(ii,:); % Save spike time
end

% Granule cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ii = 1:param.nGraprox
%     if param.ProximalON == true
%         GraProximal{ii}.V = Vgraprox(ii,:); % Save neuronal activity
%         GraProximal{ii}.S = Sgraprox(ii,:); % Save spike time
%     else
%         GraProximal{ii}.V = zeros(1,round(param.tsim / param.dt)); % 0 vector
%         GraProximal{ii}.S = zeros(1,round(param.tsim / param.dt)); % 0 vector
%     end
% end
for ii = 1:param.nGradist
    if param.DistalON == true
        GraDistal{ii}.V = Vgradist(ii,:); % Save neuronal activity
        GraDistal{ii}.S = Sgradist(ii,:); % Save spike time
    else
        GraDistal{ii}.V = zeros(1,round(param.tsim / param.dt)); % 0 vector
        GraDistal{ii}.S = zeros(1,round(param.tsim / param.dt)); % 0 vector
    end
end


% Save input currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to proximal granule
if param.ProximalON == true
    if param.DistalON == true
        InputCurrent.Idistprox = Idistprox_matrix;
    else
        InputCurrent.Imitgraprox = Imitgraprox_matrix;
    end
    InputCurrent.Igragra = Igragra_matrix;
end

% to distal granule
if param.DistalON == true
    if param.ProximalON == true
        InputCurrent.Iproxdist = Iproxdist_matrix;
    end
    InputCurrent.ImitgradistAMPA = ImitgradistAMPA_matrix;
    InputCurrent.ImitgradistNMDA = ImitgradistNMDA_matrix;
    InputCurrent.ImitgradistVDCC = ImitgradistVDCC_matrix;
end

% to mitral
if param.DistalON == true
    InputCurrent.Igradistmit = Igradistmit_matrix;
else
    InputCurrent.Igraproxmit = Igraproxmit_matrix;
end
InputCurrent.Iext = Iext_matrix;


% Calcium currents
InputCurrent.Prelease = Prelease_matrix;
% InputCurrent.Pfiregraprox = Pfiregraprox;
InputCurrent.CCa = CCa;
InputCurrent.CCaBase = CCaBase;
InputCurrent.IVDCCBase = IVDCCBase;
InputCurrent.mgradist = m;
InputCurrent.hgradist = h;

end



function Iext = SetExtInput(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function sets the external input to MCs
% The main function of this program is NeuroActivity.m
%
% Boleszek Osinski
% 03/05/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iext = zeros(param.nMitral,round(param.tsim / param.dt));

Wext = 0.001*(sort(rand(1,param.nMitral),'descend'))+param.Wmin; % hard-coded weights matrix

for ii = 1:param.nMitral
    Iext(ii,round(param.tinit / param.dt) : round(param.tfinal / param.dt)) = ...
        Wext(ii) * (1 + param.Inoise*randn);
end

end

function R = CreateRespFreq(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates an artificial respiration to modulate the bulb
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 04/06/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = param.dt:param.dt:param.tsim;
time = time / 1000; %converting the time to seconds
R = -cos(2 * pi * param.RespFreq * time);
R = R + 1;
R = R / 2;
end

function maxg = getmaxg(dt,tau1,tau2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the max amplitude for g, so we can use this value to
% normalize the curve
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/10/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 100;

x = 0:dt:tmax;
maxg = zeros(length(tau1),1);

for ii = 1: length(tau1)
    if ii == 1
        y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) * (exp(-(x)...
            / tau1(ii)) - exp(-(x) / tau2(ii)));
        maxg(ii) = max(y);
    else
        if tau1(ii) == tau1(ii -1) && tau2(ii) == tau2(ii -1)
            maxg(ii) = maxg(ii - 1);
        else
            y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii)))...
                * (exp(-(x) / tau1(ii)) - exp(-(x) / tau2(ii)));
            maxg(ii) = max(y);
        end
    end
end
end


function Ic = SetInoSpike_GraMit(gmax,P,E,V,W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the distal granule summed graded input to mitral cells
% (weight matrix W is not necessarily square)
%
%
% Modifed by Boleszek Osinski on 07/16/2013 to allow for different numbers
% of neurons (particularly for nonspiking input from Gradist to Mitral)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: if W  = wGraMit then it has dimensions nMitral x nGranule

%     g = P .* gmax;
Ic = zeros(size(W,1),1);
    for ii = 1:size(W,1)
    Ic(ii) = sum(W(ii,:)' .* P * gmax(ii) * (E(ii) - V(ii)));
    end
    
end

function Ic = SetI(tau1,tau2,t,t0,normg,W,E,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the currents at each time step.
% This function is called from NeuroActivity.m
%
% Licurgo de Almeida
% 11/18/2010
%
% Ic: channel current
% gmax: maximum conductance
% tau1: channel's rising time
% tau2: channel's falling time
% t: step time
% t0: last time the pre synaptic neuron fired
% normg: normalizes the conductance curve between 0 and 1
% W: synaptic weights
% E: Nernst potential
% V: neuron's potential from last timestep
% Mcon: connection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = (((tau1 * tau2) / (tau1 - tau2)) *...
    (exp(-(t - t0) / tau1) - exp(-(t - t0) / tau2))) / normg;

Ic = (W * g) .* (E - V);

end


function INMDA = SetI_NMDA(tau1,tau2,t,t0,gmax,W,E,V)
% NMDA Synapse
g_norm = 1;
Mg_conc = 1;
gamma = 0.016;
eta = 0.28;
E_Mg = 0; % assuming that Vrest = -70e-3

% % Code for visualizing all the variables of this current
% t1 = 2;
% t2 = 75;
% 
% t = 0:1:200;
% 
% v = (E_Mg-0.1):0.001:(E_Mg+0.1);
% 
% g = zeros(length(t), length(v));
% for i = 1:length(t)
%     for j = 1:length(v)
%         g(i,j) = g_norm * (exp(-t(i)/t2) - exp(-t(i)/t1))/(1+eta*Mg_conc*exp(-(v(j)-E_Mg)/gamma));
%     end
% end
% figure(4321)
% surf(v,t,g, 'FaceColor', 'interp', 'edgecolor', 'none', 'FaceLighting', 'phong')
% colorbar;
% camlight left;
% xlabel('V_M (V)');ylabel('t (ms)');zlabel('Conductance')
% 

% % variables for debugging
% tau1 = tauNMDA1(1);
% tau2 = tauNMDA2(1);
% t0 = tNMDA0mit_gra(:,tt);
% W = wMitGradist;
% E = ENMDAgradist(1);
% V = Vgradist(:,tt - 1);

% NOTE: If W is Wmitgradist then it has dimensions ngradist x nmit

% g = g_norm * (exp(-(t - t0)./tau1) - exp(-(t - t0)./tau2)); % dim nmit x 1
% Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma)); % dim ngradist x 1
% INMDA = (W * g) .* (V-E) .* Mg_block;

% Modified to allow for distributions over tau1 and tau2
Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma)); % dim ngradist x 1

INMDA = zeros(length(V),1);
for ii = 1:length(V)
    g = (exp(-(t - t0)./tau2(ii)) - exp(-(t - t0)./tau1(ii)))./gmax; % dim nmit x 1
    INMDA(ii) = sum(W(ii,:).* g') * (E-V(ii)) * Mg_block(ii);
end

% NOTE!!! All voltages are in V, not mV, so parameters such as gamma have
% to be scaled accordingly
end


function IVDCC = SetI_VDCC(W,E,V,m,h)
% L-type HVA Synapse
% defined as
% I_L = g*m*h*(V - E_L);

% % Visualize all variables of the current
% E_L = 0;
% V = 1e3*((E_L-0.1):0.001:(E_L+0.1)); % in mV
% 
% Compare L and N-type Ca models
% % Look at gating variable m
% % m' = (mbar - m)/taum;
% mbarL = 1./(1 + exp(-(V+50)./3));
% mbarN = 1./(1 + exp(-(V+45)./7));
% taumL = 18 * exp(-((V+45)./20).^2) + 1.5;
% taumN = 18 * exp(-((V+70)./25).^2) + 0.3;
% subplot(2,1,1)
% plot(V,mbarL,V,mbarN);xlim([-70 0])
% set(gca,'fontsize',14)
% legend('L-Type','N-Type','location','best')
% title('m')
% subplot(2,1,2)
% plot(V,taumL,V,taumN);xlim([-70 0])
% set(gca,'fontsize',14)
% title('\tau_m (ms)')
% xlabel('Vm_{GC} (mV)')
% 
% % plot N-type and NMDA together
% v = -0.1:0.001:0.1; % in V
% Mg_conc = 1;
% E_Mg = 0;
% gamma = 0.016;
% eta = 0.28;
% Mg_block_original = 1./(1+eta*Mg_conc*exp(-(v-E_Mg)/gamma));
% 
% scrsz = get(0,'ScreenSize');
% figH=figure;
% set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
% plot(V,mbarN,'-k',V,Mg_block_original,'.k');xlim([-80 0])
% set(gca,'fontsize',17)
% % legend('Bo''s model','NMDA')
% legend('N-Type','NMDA','location','southeast')
% legend boxoff
% xlabel('V_{rest,GC} (mV)');ylabel('Activation')


% 
% subplot(2,1,1)
% plot(V,mbar);xlim([-70 0])
% set(gca,'fontsize',17)
% xlabel('V (mV)');ylabel('Steady state m')
% subplot(2,1,2)
% plot(V,taum);xlim([-70 0])
% set(gca,'fontsize',17)
% xlabel('V (mV)');ylabel('\tau_{m} (ms)')
% 
% % Look at product of gating variables m*h
% Ccytmin = 0.1;
% Ccytmax = 0.5; % cytosolic [Ca] varries between 0.1 and 0.5 uM durring oscillations
% hmax = 0.00045/(0.00045+Ccytmax);
% hmin = 0.00045/(0.00045+Ccytmin);
% 
% plot(V,hmin*mbar,V,hmax*mbar);xlim([-70 0])
% set(gca,'fontsize',17)
% legend('m*h_{min}','m*h_{max}')
% legend boxoff
% xlabel('V (mV)');ylabel('Steady state m*h')

% Nernst potential for Ca
% RT = 300*8.31; % J/mole
% z= 2 ; % Ca ion valence
% F = 96485; % Faraday constant Coul/mole
% Cout = 1500;
% Cin = 0.5;
% 
% ECa = (RT/(z*F))*log(Cout/Cin)
% % ECa = 100 - 130 mV!


% I_L = g*m*h*(V - E_L)*s(t);
% NOTE: W has dimensions ngradist x nmit
% m has dimensions ngradist x 1

% Cin = 0.25; % cytosolic Ca concentration (uM)
% Cin = Cin(1);
% h = 0.00045/(0.00045+Cin);

IVDCC = zeros(length(V),1);
for ii = 1:length(V)
    IVDCC(ii) = sum(W(ii,:)) * m(ii) * h(ii) * (E(ii)-V(ii));
 %    IVDCC(ii) = max(W(ii,:)) * m(ii) * h(ii) * (E(ii)-V(ii));
end

% NOTE!!! All voltages are in V, not mV, so parameters have
% to be scaled accordingly
end

function P = Prelease(C,B,T)
% GABA release probablity
% 
% Boleslaw Osinski
% March 2015


% matrix of probabilities forced to be within range 0 - 1
% C - [Ca]
% B - baseline [Ca]
% T - threshold for maximum GABA release

P = (C - B) ./ (T - B);
% P = C ./ T;
J = P <= 0;
P(J) = 0;
J = P > 1;
P(J) = 1;

end




function Mat = SetConnections(cell1,cell2,cchance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the connection between different neurons
% The main function of this program is piriformmain.m
%
% Modified by Boleszek Osinski on 07/15/2013
%
% Licurgo de Almeida (original)
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: this function only sets fixed connections
    
        Mat = zeros(cell1,cell2);
            for ii = 1:cell1
                convector = randperm(cell2);
                convector = convector(1:round(cell2 * cchance));
                Mat(ii,convector) = 1;
            end

end

function w = SetWeights(Mat,weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the synaptic weights between different neurons
%
% % Modified by Boleszek Osinski on 07/15/2013
%
% Licurgo de Almeida (original)
% 03/02/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = zeros(size(Mat));
for ii = 1:length(weights)
    w(ii,:) = Mat(ii,:) * weights(ii);
end


end



