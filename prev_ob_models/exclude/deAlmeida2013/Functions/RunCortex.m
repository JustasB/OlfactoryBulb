function [Ff,Py,Fb] = RunCortex(Mi,Ff,Py,Fb)
% [Ff,Py,Fb] = RunCortex(Mi,Ff,Py,Fb)
% This function calculates cell activity over time. The cell groups are:
% Mi: Mitral cell (soma compartment) coming from the bulb simulations
% Ff: Feedforward cells (interneuron)
% Py: Pyramidal cells
% Fb: Feedback cells (interneuron)

% Get Mitral output data
Om = Mi.O;
Tm = Mi.TFire;
dt = Mi.dt;

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

for tt = 1:size(Om,2)
    % Mitral input
    Out = Om(:,tt) == 1;
    Tm(Out) = 0;
    Tm(~Out) = Tm(~Out) + dt;
    
    % Feedforward cell
    Cfeff = SetISpike(Vf,tau1feff,tau2feff,Tm,Wfeff,Efeff);
    Vf = Vf + stepf * (Cfeff - Vf);
    [Vf,Tf] = SpikeOutput(Vf,Tf,Vrestf,FThreshf,Outfline{:});
    Of(:,tt) = Tf == 0;
    Vof(:,tt) = Vf;

    % Pyramidal cell
    Cpeff = SetISpike(Vp,tau1peff,tau2peff,Tm,Wpeff,Epeff);
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
Ff.O = Of;
Ff.V = Vof;
Py.O = Op;
Py.V = Vop;
Fb.O = Ob;
Fb.V = Vob;