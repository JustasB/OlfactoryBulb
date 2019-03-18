%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% classPyramidal includes all proprieties and methods specific for
% pyramidal cells.
%
% Licurgo de Almeida
% 04/29/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef classPyramidal < classSpkNeuron
    properties
        % from Mitral cells
        AMPAFf % struct with parameters from the excitatory synapses
        ConnAMPAFf  = 0.2; % percentage of pyr cells connected to each
        % mitral
        MAMPAFf % connection matrix with the mitral cells
        WAMPAFf % synaptic weight matrix with the mitral cells
        
        % from Feed forward interneurons
        GABAFf % struct with parameters from the inhibitory synapses
        ConnGABAFf  = 0.3; % percentage of Ff cells connected to each
        % pyramidal cell
        MGABAFf % connection matrix with the Feedforward interneurons
        WGABAFf % synaptic weight matrix with the Feedforward interneurons
        
        % from Pyramidal cells
        AMPAFb % struct with parameters from the excitatory synapses
        ConnAMPAFb  = 0.2; % percentage of pyr cells connected to each
        % other
        MAMPAFb % connection matrix with the pyramidal cells
        WAMPAFb % synaptic weight matrix with the pyramidal cells
        
        % from Feedback interneurons
        GABAFb % struct with parameters from the inhibitory synapses
        ConnGABAFb  = 0.4; % percentage of feedback cells connected to each
        % pyramidal cell
        MGABAFb % connection matrix with the Feedback interneurons
        WGABAFb % synaptic weight matrix with the Feedback interneurons
        
        % from AHP
        EAHP = -15e-3; % reversal potential for after hyperpolarization
        % potential
        AAHP = 16; % amplitude of the hyperpolarization potential
        tauAHP = 20; % falling time of the hyperpolarization potential
        
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
        function obj = classPyramidal(tsim,ncells)
            if nargin == 0
                inputsuper = {};
            else
                inputsuper = {tsim,ncells};
            end
            obj = obj@classSpkNeuron(inputsuper{:});
            obj.tau = 10; %ms
            obj.CellName = 'Pyramidal';
            obj.AMPAFf = struct('E',70e-3,'tau1',1,'tau2',2,'G',0.38);
            % where the elements of the struct are:
            % E: reversal potential
            % tau1: rising time of the conductance
            % tau2: falling time of the conductance
            % G: max conductance
            obj.MAMPAFf = obj.SetConnections(obj.ncells,obj.ConnAMPAFf,'normal');
            obj.WAMPAFf = obj.MAMPAFf;
            obj.GABAFf = struct('E',-15e-3,'tau1',4,'tau2',8,'G',0.38);
            obj.MGABAFf = obj.SetConnections(obj.ncells,obj.ConnGABAFf);
            obj.WGABAFf = obj.MGABAFf;
            obj.AMPAFb = struct('E',70e-3,'tau1',1,'tau2',2,'G',[4.5,1.8]);
            % AMPAFb.G changes with ACh modulation
            obj.MAMPAFb = obj.SetConnections(obj.ncells,obj.ConnAMPAFb,'auto');
            obj.WAMPAFb = rand(size(obj.MAMPAFb)) .* obj.MAMPAFb * 0.01; % autoassociative
            % connection weights start with 0
            obj.GABAFb = struct('E',-15e-3,'tau1',4,'tau2',8,'G',[200e-3,70e-3]);
            obj.MGABAFb = obj.SetConnections(obj.ncells,obj.ConnGABAFb);
            obj.WGABAFb = obj.MGABAFb;
        end
    end
end