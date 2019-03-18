function [OSN,Pg,M1,M2,Gr,Ff,Py,Fb] = RunOlfacSys(OSN,Pg,M1,M2,Gr,Ff,Py,Fb)
% [OSN,Pg,M1,M2,Gr,Ff,Py,Fb] = RunOlfacSys(OSN,Pg,M1,M2,Gr,Ff,Py,Fb)
% This function calculates cell activity over time. The cell groups are:
% OSN: Olfactory sensory neurons
% Pg: Periglomerular cells (interneuron)
% M1: Apical compartment of Mitral cells
% M2: Soma compartment of Mitral cells
% Gr: Granule cells (interneuron)
% Ff: Feedforward cells (interneuron)
% Py: Pyramidal cells
% Fb: Feedback cells (interneuron)








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
Vgr = Gr.V(:,end);
Ogr = zeros(size(Ioi));
Vogr = zeros(size(Ioi));
stepgr = Gr.dt / Gr.tau; % timestep
Wgreff = Gr.WAMPAFf * Gr.AMPAFf.G;
Egreff = Gr.AMPAFf.E;
tau1greff = Gr.AMPAFf.tau1;
tau2greff = Gr.AMPAFf.tau2;
Wgrifb = Gr.WGABAFb * Gr.GABAFb.G;
Egrifb = Gr.GABAFb.E;
tau1grifb = Gr.GABAFb.tau1;
tau2grifb = Gr.GABAFb.tau2;
Tgr = Gr.TFire * 0; % Gr Tfire starts at 0
Modgr = CalculateAffinity(Gr.K,Gr.Cmax,Gr.ModValue,Gr.b);
Modnegr = CalculateAffinity(Gr.Kne,Gr.Cmaxne,Gr.ModValue,Gr.bne);
Gr.Mod(:) = Modgr;
Vrestgr = SetModulation(Gr.Tmin,Modgr) + SetModulation(Gr.Tminne,Modnegr);
FThreshgr = SetModulation(Gr.Tmax,Modgr);
Outgrline = {Gr.Trefrac,Gr.Beta,Gr.dt,Gr.SpikeV,Gr.Vhyper};

% Create Feedforward variables.
Vf = Ff.V(:,end);
Of = Ff.O;
Vof = Ff.V;
stepf = Ff.dt / Ff.tau; % timestep
Wfeff = Ff.WAMPAFf * Ff.AMPAFf.G;
Efeff = Ff.AMPAFf.E;
tau1feff = Ff.AMPAFf.tau1;
tau2feff = Ff.AMPAFf.tau2;
Tf = Ff.TFire;
Modf = CalculateAffinity(Ff.K,Ff.Cmax,Ff.ModValue,Ff.b);
Ff.Mod(:) = Modf;
Vrestf = SetModulation(Ff.Tmin,Modf);
FThreshf = SetModulation(Ff.Tmax,Modf);
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

Tp = Py.TFire;
Tminp = Py.Tmin;
Tmaxp = Py.Tmax;

Modp = CalculateAffinity(Py.K,Py.Cmax,Py.ModValue,Py.b);

Py.Mod(:) = Modp;
Vrestp = SetModulation(Tminp,Modp);
FThreshp = SetModulation(Tmaxp,Modp);
Gpefb = SetModulation(AMPAFbGp,Modp);
Gpifb = SetModulation(GABAFbGp,Modp);
Apahp = SetModulation(AAHPp,Modp);

Outpline = {Py.Trefrac,Py.Beta,Py.dt,Py.SpikeV,Py.Vhyper};

% Create Feedback variables
Vb = Fb.V(:,end);
Ob = Fb.O;
Vob = Fb.V;
stepb = Fb.dt / Fb.tau; % timestep

Wbeff = Fb.WAMPAFf;
AMPAFfGb = Fb.AMPAFf.G; % AMPAFf.G changes with ACh modulation
Ebeff = Fb.AMPAFf.E;
tau1beff = Fb.AMPAFf.tau1;
tau2beff = Fb.AMPAFf.tau2;

Tb = Fb.TFire;
Tminb = Fb.Tmin;
Tmaxb = Fb.Tmax;

Modb = CalculateAffinity(Fb.K,Fb.Cmax,Fb.ModValue,Fb.b);
Fb.Mod(:) = Modb;
Vrestb = SetModulation(Tminb,Modb);
FThreshb = SetModulation(Tmaxb,Modb);
Gbeff = SetModulation(AMPAFfGb,Modb);

Outbline = {Fb.Trefrac,Fb.Beta,Fb.dt,Fb.SpikeV,Fb.Vhyper};

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
    Cm2ifb = SetISpike(Vm2,tau1m2ifb,tau2m2ifb,Tgr,Wm2ifb,Em2ifb);
    Vm2 = Vm2 + stepm2 * ((Vm1 - Vm2) / ((Rm1eff + Rm2eff) / 2) + Cm2ifb - Vm2);
    [Vm2,Tm2] = SpikeOutput(Vm2,Tm2,Vrestm2,FThreshm2,Outm2line{:});
    Om2(:,tt) = Tm2 == 0;
    Vom2(:,tt) = Vm2;
    
    %Granule
    Cgreff = SetISpike(Vgr,tau1greff,tau2greff,Tm2,Wgreff,Egreff);
    Cgrifb = SetISpike(Vgr,tau1grifb,tau2grifb,Tgr,Wgrifb,Egrifb);
    Vgr = Vgr + stepgr * (Cgreff + Cgrifb - Vgr);
    [Vgr,Tgr] = SpikeOutput(Vgr,Tgr,Vrestgr,FThreshgr,Outgrline{:});
    Ogr(:,tt) = Tgr == 0;
    Vogr(:,tt) = Vgr;

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
    Vp = Vp + stepp * (Cpeff + Cpiff + Cpefb + Cpifb - Vp +...
        Vpahp .* (Epahp - Vp));
    [Vp,Tp] = SpikeOutput(Vp,Tp,Vrestp,FThreshp,Outpline{:});
    Op(:,tt) = Tp == 0;
    Vop(:,tt) = Vp;

    % Feedback cell
    Cbeff = SetISpike(Vb,tau1beff,tau2beff,Tp,Wbeff * Gbeff,Ebeff);
    Vb = Vb + stepb * (Cbeff - Vb);
    [Vb,Tb] = SpikeOutput(Vb,Tb,Vrestb,FThreshb,Outbline{:});
    Ob(:,tt) = Tb == 0;
    Vob(:,tt) = Vb;

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
Gr.O = Ogr;
Gr.V = Vogr;
Ff.O = Of;
Ff.V = Vof;
Py.O = Op;
Py.V = Vop;
Py.O(:,1:10) = 0; % this to avoid a intial burst of almost all neurons in the
% network
Fb.O = Ob;
Fb.V = Vob;