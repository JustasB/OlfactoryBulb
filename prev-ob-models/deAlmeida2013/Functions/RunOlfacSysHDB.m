function [OSN,Pg,M1,M2,Gr,Ff,Py,Fb,Gb,Ac] = RunOlfacSysHDB(OSN,Pg,M1,M2,Gr,Ff,Py,Fb,Gb,Ac)
% [OSN,Pg,M1,M2,Gr,Ff,Py,Fb] = RunOlfacSysHDB(OSN,Pg,M1,M2,Gr,Ff,Py,Fb)
% This function calculates cell activity over time. The cell groups are:
% OSN: Olfactory sensory neurons
% Pg: Periglomerular cells (interneuron)
% M1: Apical compartment of Mitral cells
% M2: Soma compartment of Mitral cells
% Gr: Granule cells (interneuron)
% Ff: Feedforward cells (interneuron)
% Py: Pyramidal cells
% Fb: Feedback cells (interneuron)
% Gb: HDB Gabaergic (interneuron)
% Ac: HDB Cholinergic

% Create OSNs variables. These variables are created to reduce overhead and
% improve performance
Woi = OSN.WInput * OSN.AMPAFf.G;
Vo = OSN.V(:,end);
Ioi = OSN.OdorInput;
Eoi = OSN.AMPAFf.E;
Oo = zeros(size(Ioi));
Voo = zeros(size(Ioi));
Tmino = OSN.Tmin;
Tmaxo = OSN.Tmax;
if OSN.ModValue == -1 % if mod is -1, the level of modulation changes at
    % each time step
    dynmod = true;
    Ko = OSN.K;
    bo = OSN.b;
else
    dynmod = false;
    Modo = ModulationActivation(OSN.K,OSN.ModValue,OSN.b);
    OSN.Mod(:) = Modo;
    Vresto = SetModulation(Tmino,Modo);
    FThresho = SetModulation(Tmaxo,Modo);
end
Betao = OSN.Beta;
stepo = OSN.dt / OSN.tau; %timestep

% Create Pglo cells variables
Wpgff = Pg.WAMPAFf * Pg.AMPAFf.G;
Vpg = Pg.V(:,end);
Epgff = Pg.AMPAFf.E;
Opg = zeros(size(Ioi));
Vopg = zeros(size(Ioi)); % pg cells fully active.
Tminpg = Pg.Tmin;
Tmaxpg = Pg.Tmax;
if dynmod == true
    Kpg = Pg.K;
    bpg = Pg.b;
else
    Modpg = ModulationActivation(Pg.K,OSN.ModValue,Pg.b);
    Pg.Mod(:) = Modpg;
    Vrestpg = SetModulation(Tminpg,Modpg);
    FThreshpg = SetModulation(Tmaxpg,Modpg);
end
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
Tminm2 = M2.Tmin;
Tmaxm2 = M2.Tmax;
if dynmod == true
    Km2 = M2.K;
    bm2 = M2.b;
else
    Modm2 = ModulationActivation(M2.K,OSN.ModValue,M2.b);
    M2.Mod(:) = Modm2;
    Vrestm2 = SetModulation(Tminm2,Modm2);
    FThreshm2 = SetModulation(Tmaxm2,Modm2);
end
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
Tming = Gr.Tmin;
Tmaxg = Gr.Tmax;
if dynmod == true
    Kg = Gr.K;
    bg = Gr.b;
else
    Modg = ModulationActivation(Gr.K,OSN.ModValue,Gr.b);
    Gr.Mod(:) = Modg;
    Vrestg = SetModulation(Tming,Modg);
    FThreshg = SetModulation(Tmaxg,Modg);
end
Outgline = {Gr.Trefrac,Gr.Beta,Gr.dt,Gr.SpikeV,Gr.Vhyper};

% Create Feedforward variables. These variables are created to reduce
% overhead and improve performance
Vf = Ff.V(:,end);
Of = Ff.O;
Vof = Ff.V;
stepf = Ff.dt / Ff.tau; % timestep

Wfeff = Ff.WAMPAFf * Ff.AMPAFf.G;
Efeff = Ff.AMPAFf.E;
tau1feff = Ff.AMPAFf.tau1;
tau2feff = Ff.AMPAFf.tau2;
Tf = Ff.TFire;
Tminf = Ff.Tmin;
Tmaxf = Ff.Tmax;
if dynmod == true
    Kf = Ff.K;
    bf = Ff.b;
else
    Modf = ModulationActivation(Ff.K,OSN.ModValue,Ff.b);
    Ff.Mod(:) = Modf;
    Vrestf = SetModulation(Tminf,Modf);
    FThreshf = SetModulation(Tmaxf,Modf);
end
Outfline = {Ff.Trefrac,Ff.Beta,Ff.dt,Ff.SpikeV,Ff.Vhyper};

% Create Pyramidal variables
Vp = Py.V(:,end);
Op = Py.O;
Vop = Py.V;
stepp = Py.dt / Py.tau; % timestep

Wpeff = Py.WAMPAFf * Py.AMPAFf.G;
Epeff = Py.AMPAFf.E;
tau1peff = Py.AMPAFf.tau1;
tau2peff = Py.AMPAFf.tau2;

Wpiff = Py.WGABAFf * Py.GABAFf.G;
Epiff = Py.GABAFf.E;
tau1piff = Py.GABAFf.tau1;
tau2piff = Py.GABAFf.tau2;

if Py.LearnUnlearn == true
    Wpefb = Py.WAMPAFb;
else
    Wpefb = Py.WAMPAFb / sum(sum(Py.WAMPAFb)); % if LearnUnlearn == false
    % Ws are normalized
end
AMPAFbGp = Py.AMPAFb.G; % AMPAFb.G changes with ACh modulation
Epefb = Py.AMPAFb.E;
tau1pefb = Py.AMPAFb.tau1;
tau2pefb = Py.AMPAFb.tau2;

Wpifb = Py.WGABAFb;
GABAFbGp = Py.GABAFb.G; % GABAFb.G changes with ACh modulation
Epifb = Py.GABAFb.E;
tau1pifb = Py.GABAFb.tau1;
tau2pifb = Py.GABAFb.tau2;

Epahp = Py.EAHP;
AAHPp = Py.AAHP; % Py.AAHP changes with ACh modulation
Vpahp = zeros(size(Vp)); % AHP
steppahp = Py.dt / Py.tauAHP;

%Apahp = SetModulation(Py.AAHP,Py.Mod);


Tp = Py.TFire;
Tminp = Py.Tmin;
Tmaxp = Py.Tmax;
if dynmod == true
    Kp = Py.K;
    bp = Py.b;
else
    Modp = ModulationActivation(Py.K,OSN.ModValue,Py.b);
    Py.Mod(:) = Modp;
    Vrestp = SetModulation(Tminp,Modp);
    FThreshp = SetModulation(Tmaxp,Modp);
    Gpefb = SetModulation(AMPAFbGp,Modp);
    Gpifb = SetModulation(GABAFbGp,Modp);
    Apahp = SetModulation(AAHPp,Modp);
end
Outpline = {Py.Trefrac,Py.Beta,Py.dt,Py.SpikeV,Py.Vhyper};

% Create Feedback variables
Vb = Fb.V(:,end);
Ob = Fb.O;
Vob = Fb.V;
stepb = Fb.dt / Fb.tau; % timestep

Wbeff = Fb.WAMPAFf * Fb.AMPAFf.G;
Ebeff = Fb.AMPAFf.E;
tau1beff = Fb.AMPAFf.tau1;
tau2beff = Fb.AMPAFf.tau2;

Tb = Fb.TFire;
Tminb = Fb.Tmin;
Tmaxb = Fb.Tmax;
if dynmod == true
    Kb = Fb.K;
    bb = Fb.b;
else
    Modb = ModulationActivation(Fb.K,OSN.ModValue,Fb.b);
    Fb.Mod(:) = Modb;
    Vrestb = SetModulation(Tminb,Modb);
    FThreshb = SetModulation(Tmaxb,Modb);
end
Outbline = {Fb.Trefrac,Fb.Beta,Fb.dt,...
    Fb.SpikeV,Fb.Vhyper};

% Create HDB Gabaergic variables. These variables are created to reduce
% overhead and improve performance
Vgb = Gb.V(:,end);
Ogb = Gb.O;
Vogb = Gb.V;
stepgb = Gb.dt / Gb.tau; % timestep

Wgbeff = Gb.WAMPAFf * Gb.AMPAFf.G;
Egbeff = Gb.AMPAFf.E;
tau1gbeff = Gb.AMPAFf.tau1;
tau2gbeff = Gb.AMPAFf.tau2;

Tgb = Gb.TFire;
Tmingb = Gb.Tmin;
Tmaxgb = Gb.Tmax;
if dynmod == true
    Kgb = Gb.K;
    bgb = Gb.b;
else
    Modgb = ModulationActivation(Gb.K,OSN.ModValue,Gb.b);
    Gb.Mod(:) = Modgb;
    Vrestgb = SetModulation(Tmingb,Modgb);
    FThreshgb = SetModulation(Tmaxgb,Modgb);
end
Outgbline = {Gb.Trefrac,Gb.Beta,Gb.dt,Gb.SpikeV,Gb.Vhyper};

% Create HDB Cholinergic variables. These variables are created to reduce
% overhead and improve performance
Va = Ac.V(:,end);
Oa = Ac.O;
Voa = Ac.V;
stepa = Ac.dt / Ac.tau; % timestep

stepActa = Ac.dt / Ac.TauAct;
Aa = Ac.InitalAct;
Acta = Ac.ActOverT;

Waeff = Ac.WAMPAFf * Ac.AMPAFf.G;
Eaeff = Ac.AMPAFf.E;
tau1aeff = Ac.AMPAFf.tau1;
tau2aeff = Ac.AMPAFf.tau2;

Waiff = Ac.WGABAFf * Ac.GABAFf.G;
Eaiff = Ac.GABAFf.E;
tau1aiff = Ac.GABAFf.tau1;
tau2aiff = Ac.GABAFf.tau2;

Ta = Ac.TFire;
Tmina = Ac.Tmin;
Tmaxa = Ac.Tmax;
if dynmod == true
    Ka = Ac.K;
    ba = Ac.b;
else
    Moda = ModulationActivation(Ac.K,OSN.ModValue,Ac.b);
    Ac.Mod(:) = Moda;
    Vresta = SetModulation(Tmina,Moda);
    FThresha = SetModulation(Tmaxa,Moda);
end
Outaline = {Ac.Trefrac,Ac.Beta,Ac.dt,Ac.SpikeV,Ac.Vhyper};
% Part that matters
for tt = 1:size(OSN.OdorInput,2)
    %Set dynamic modulation of all networks
    if dynmod == true
        %OSN
        Modo = ModulationActivation(Ko,Aa,bo);
        Vresto = SetModulation(Tmino,Modo);
        FThresho = SetModulation(Tmaxo,Modo);
        %Pglo
        Modpg = ModulationActivation(Kpg,Aa,bpg);
        Vrestpg = SetModulation(Tminpg,Modpg);
        FThreshpg = SetModulation(Tmaxpg,Modpg);
        %Mitral soma compartment
        Modm2 = ModulationActivation(Km2,Aa,bm2);
        Vrestm2 = SetModulation(Tminm2,Modm2);
        FThreshm2 = SetModulation(Tmaxm2,Modm2);
        %Granule
        Modg = ModulationActivation(Kg,Aa,bg);
        Vrestg = SetModulation(Tming,Modg);
        FThreshg = SetModulation(Tmaxg,Modg);
        % Feedforward cell
        Modf = ModulationActivation(Kf,Aa,bf);
        Vrestf = SetModulation(Tminf,Modf);
        FThreshf = SetModulation(Tmaxf,Modf);
        % Pyramidal cell
        Modp = ModulationActivation(Kp,Aa,bp);
        Vrestp = SetModulation(Tminp,Modp);
        FThreshp = SetModulation(Tmaxp,Modp);
        Gpefb = SetModulation(AMPAFbGp,Modp);
        Gpifb = SetModulation(GABAFbGp,Modp);
        Apahp = SetModulation(AAHPp,Modp);
        % Feedback cell
        Modb = ModulationActivation(Kb,Aa,bb);
        Vrestb = SetModulation(Tminb,Modb);
        FThreshb = SetModulation(Tmaxb,Modb);
        % HDB Gabaergic cell
        Modgb = ModulationActivation(Kgb,Aa,bgb);
        Vrestgb = SetModulation(Tmingb,Modgb);
        FThreshgb = SetModulation(Tmaxgb,Modgb);
        % HDB Cholinergic cell
        Moda = ModulationActivation(Ka,Aa,ba);
        Vresta = SetModulation(Tmina,Moda);
        FThresha = SetModulation(Tmaxa,Moda);
    end
    
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
    
    % Feedforward cell
    Cfeff = SetISpike(Vf,tau1feff,tau2feff,Tm2,Wfeff,Efeff);
    Vf = Vf + stepf * (Cfeff - Vf);
    [Vf,Tf] = SpikeOutput(Vf,Tf,Vrestf,FThreshf,Outfline{:});
    Of(:,tt) = Tf == 0;
    Vof(:,tt) = Vf;
    
    % Pyramidal cell
    Cpeff = SetISpike(Vp,tau1peff,tau2peff,Tm2,Wpeff,Epeff);
    Cpiff = SetISpike(Vp,tau1piff,tau2piff,Tf,Wpiff,Epiff);
    Cpefb = SetISpike(Vp,tau1pefb,tau2pefb,Tp,Wpefb * Gpefb,Epefb);
    Cpifb = SetISpike(Vp,tau1pifb,tau2pifb,Tb,Wpifb * Gpifb,Epifb);
    Vpahp = Vpahp + steppahp * ((Tp == 0) * Apahp - Vpahp); % AHP
    %XX(:,tt) = Vpahp .* (Epahp - Vp);
    %YY(:,tt) = Cpeff;
    Vp = Vp + stepp * (Cpeff + Cpiff + Cpefb + Cpifb - Vp +...
        Vpahp .* (Epahp - Vp));
    [Vp,Tp] = SpikeOutput(Vp,Tp,Vrestp,FThreshp,Outpline{:});
    Op(:,tt) = Tp == 0;
    Vop(:,tt) = Vp;
    
    % Feedback cell
    Cbeff = SetISpike(Vb,tau1beff,tau2beff,Tp,Wbeff,Ebeff);
    Vb = Vb + stepb * (Cbeff - Vb);
    [Vb,Tb] = SpikeOutput(Vb,Tb,Vrestb,FThreshb,Outbline{:});
    Ob(:,tt) = Tb == 0;
    Vob(:,tt) = Vb;
    
    % HDB Gabaergic cell
    Cgbeff = SetISpike(Vgb,tau1gbeff,tau2gbeff,Tp,Wgbeff,Egbeff);
    Vgb = Vgb + stepgb * (Cgbeff - Vgb);
    [Vgb,Tgb] = SpikeOutput(Vgb,Tgb,Vrestgb,FThreshgb,Outgbline{:});
    Ogb(:,tt) = Tgb == 0;
    Vogb(:,tt) = Vgb;
    
    % HDB Cholinergic cell
    Caeff = SetISpike(Va,tau1aeff,tau2aeff,Tm2,Waeff,Eaeff);
    Caiff = SetISpike(Va,tau1aiff,tau2aiff,Tgb,Waiff,Eaiff);
    Va = Va + stepa * (Caeff + Caiff - Va);
    [Va,Ta] = SpikeOutput(Va,Ta,Vresta,FThresha,Outaline{:});
    Oa(:,tt) = Ta == 0;
    Voa(:,tt) = Va;
    Aa = Aa + stepActa * (sum(Oa(:,tt)) - Aa);
    Acta(tt) = Aa;
        
    
end
OSN.O = Oo;
OSN.V = Voo;
Pg.O = Opg;
Pg.V = Vopg;
% M1 doesn't have an output, only the Voltage
M1.V = Vom1;
M2.O = Om2;
M2.V = Vom2;
Gr.O = Og;
Gr.V = Vog;
Ff.O = Of;
Ff.V = Vof;
Py.O = Op;
Py.V = Vop;
Fb.O = Ob;
Fb.V = Vob;
Gb.O = Ogb;
Gb.V = Vogb;
Ac.O = Oa;
Ac.V = Voa;
Ac.ActOverT = Acta;
if dynmod == true
    OSN.Mod = ModulationActivation(Ko,Acta,bo);
    Pg.Mod = ModulationActivation(Kpg,Acta,bpg);
    M2.Mod = ModulationActivation(Km2,Acta,bm2);
    Gr.Mod = ModulationActivation(Kg,Acta,bg);
    Ff.Mod = ModulationActivation(Kf,Acta,bf);
    Py.Mod = ModulationActivation(Kp,Acta,bp);
    Fb.Mod = ModulationActivation(Kb,Acta,bb);
    Gb.Mod = ModulationActivation(Kgb,Acta,bgb);
    Ac.Mod = ModulationActivation(Ka,Acta,ba);
end