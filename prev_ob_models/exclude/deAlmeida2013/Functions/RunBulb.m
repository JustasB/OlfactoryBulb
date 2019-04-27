function [OSN,Pg,M1,M2,Gr] = RunBulb(OSN,Pg,M1,M2,Gr)
% [OSN,Pg,M1,M2,Gr] = RunBulb(OSN,Pg,M1,M2,Gr)
% This function calculates cell activity over time. The cell groups are:
% OSN: Olfactory sensory neurons
% Pg: Periglomerular cells
% M1: Apical compartment of Mitral cells
% M2: Soma compartment of Mitral cells
% Gr: Granule cells

% Create OSNs variables. These variables are created to reduce overhead and
% improve performance
Woi = OSN.WInput * OSN.AMPAFf.G;
Vo = OSN.V(:,end);
Ioi = OSN.OdorInput;
Eoi = OSN.AMPAFf.E;
Oo = zeros(size(Ioi));
Voo = zeros(size(Ioi));
Modo = CalculateAffinity(OSN.K,OSN.Cmax,OSN.ModValue,OSN.b);
OSN.Mod(:) = Modo;
Vresto = SetModulation(OSN.Tmin,Modo);
FThresho = SetModulation(OSN.Tmax,Modo);
Betao = OSN.Beta;
stepo = OSN.dt / OSN.tau; %timestep

% Create Pglo cells variables
Wpgff = Pg.WAMPAFf * Pg.AMPAFf.G;
Vpg = Pg.V(:,end);
Epgff = Pg.AMPAFf.E;
Opg = zeros(size(Ioi));
Vopg = zeros(size(Ioi)); % pg cells fully active.
Modpg = CalculateAffinity(Pg.K,Pg.Cmax,Pg.ModValue,Pg.b);
Pg.Mod(:) = Modpg;
Vrestpg = SetModulation(Pg.Tmin,Modpg);
FThreshpg = SetModulation(Pg.Tmax,Modpg);
Betapg = Pg.Beta;
steppg = Pg.dt / Pg.tau; % timestep

% Create apical compartment of mitral cells variables
Wm1eff = M1.WAMPAFf * M1.AMPAFf.G;
Vm1 = M1.V(:,end);
Em1eff = M1.AMPAFf.E;
Vom1 = zeros(size(Ioi));
stepm1 = M1.dt / M1.tau; % timestep
Wm1iff = M1.WGABAFf * M1.GABAFf.G;
Em1iff = M1.GABAFf.E;
Rm1eff = M1.Rsom;

% Create soma compartment of mitral cells
Rm2eff = M2.Rglo;
Vm2 = M2.V(:,end);
Om2 = zeros(size(Ioi));
Vom2 = zeros(size(Ioi));
stepm2 = M2.dt / M2.tau; % timestep
Wm2ifb = M2.WGABAFb * M2.GABAFb.G;
Em2ifb = M2.GABAFb.E;
tau1m2ifb = M2.GABAFb.tau1;
tau2m2ifb = M2.GABAFb.tau2;
Tm2 = M2.TFire;
Modm2 = CalculateAffinity(M2.K,M2.Cmax,M2.ModValue,M2.b);
M2.Mod(:) = Modm2;
Vrestm2 = SetModulation(M2.Tmin,Modm2);
FThreshm2 = SetModulation(M2.Tmax,Modm2);
Outm2line = {M2.Trefrac,M2.Beta,M2.dt,M2.SpikeV,M2.Vhyper};

% Create Granule cells
Vg = Gr.V(:,end);
Og = zeros(size(Ioi));
Vog = zeros(size(Ioi));
stepg = Gr.dt / Gr.tau; % timestep
Wgeff = Gr.WAMPAFf * Gr.AMPAFf.G;
Egeff = Gr.AMPAFf.E;
tau1geff = Gr.AMPAFf.tau1;
tau2geff = Gr.AMPAFf.tau2;
Wgifb = Gr.WGABAFb * Gr.GABAFb.G;
Egifb = Gr.GABAFb.E;
tau1gifb = Gr.GABAFb.tau1;
tau2gifb = Gr.GABAFb.tau2;
Tg = Gr.TFire * 0; % Gr Tfire starts at 0
Modg = CalculateAffinity(Gr.K,Gr.Cmax,Gr.ModValue,Gr.b);
Modneg = CalculateAffinity(Gr.Kne,Gr.Cmaxne,Gr.ModValue,Gr.bne);
Gr.Mod(:) = Modg;
Vrestg = SetModulation(Gr.Tmin,Modg) + SetModulation(Gr.Tminne,Modneg);
FThreshg = SetModulation(Gr.Tmax,Modg);
Outgline = {Gr.Trefrac,Gr.Beta,Gr.dt,Gr.SpikeV,Gr.Vhyper};

% Part that matters
for tt = 1:size(OSN.OdorInput,2)
    %OSN
    Coi = SetINoSpike(Vo,Ioi(:,tt),Woi,Eoi);
    Vo = Vo + stepo * (Coi - Vo);
    Po = NoSpikeOutput(Vo,Vresto,FThresho,Betao);
    Oo(:,tt) = Po;
    Voo(:,tt) = Vo;
    
    %Pglo
    Cpgff = SetINoSpike(Vpg,Po,Wpgff,Epgff);
    Vpg = Vpg + steppg * (Cpgff - Vpg);
    Ppg = NoSpikeOutput(Vpg,Vrestpg,FThreshpg,Betapg);
    Opg(:,tt) = Ppg;
    Vopg(:,tt) = Vo;
    
    %Mitral apical compartment
    Cm1eff = SetINoSpike(Vm1,Po,Wm1eff,Em1eff);
    Cm1iff = SetINoSpike(Vm1,Ppg,Wm1iff,Em1iff);
    Vm1 = Vm1 + stepm1 * ((Vm2 - Vm1) / ((Rm2eff + Rm1eff) / 2) + Cm1eff + Cm1iff - Vm1);
    Vom1(:,tt) = Vm1;
    
    %Mitral soma compartment
    Cm2ifb = SetISpike(Vm2,tau1m2ifb,tau2m2ifb,Tg,Wm2ifb,Em2ifb);
    Vm2 = Vm2 + stepm2 * ((Vm1 - Vm2) / ((Rm1eff + Rm2eff) / 2) + Cm2ifb - Vm2);
    [Vm2,Tm2] = SpikeOutput(Vm2,Tm2,Vrestm2,FThreshm2,Outm2line{:});
    Om2(:,tt) = Tm2 == 0;
    Vom2(:,tt) = Vm2;
    
    %Granule
    Cgeff = SetISpike(Vg,tau1geff,tau2geff,Tm2,Wgeff,Egeff);
    Cgifb = SetISpike(Vg,tau1gifb,tau2gifb,Tg,Wgifb,Egifb);
    Vg = Vg + stepg * (Cgeff + Cgifb - Vg);
    [Vg,Tg] = SpikeOutput(Vg,Tg,Vrestg,FThreshg,Outgline{:});
    Og(:,tt) = Tg == 0;
    Vog(:,tt) = Vg;
end
OSN.O = Oo;
OSN.V = Voo;
Pg.O = Opg;
Pg.V = Vopg;
% M1 doesn't have an output, only the Voltage
M1.V = Vom1;
M2.O = Om2;
M2.O(:,1:10) = 0; % this to avoid a intial burst of almost all neurons in the
% network
M2.V = Vom2;
Gr.O = Og;
Gr.V = Vog;