: Sodium ion accumulation with radial and longitudinal diffusion, buffering and pumping


NEURON {
    SUFFIX nadp
    USEION na READ nao, nai, ina WRITE nai, ina, ena
    NONSPECIFIC_CURRENT ik_pump
    RANGE ina_pmp, TotalPump, ik_ratio, na, pump, pumpna
    GLOBAL vrat, DNa, k1, k2, k3, k4
                                                           
}

DEFINE Nannuli 4

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    (mV) = (millivolt)
    FARADAY = (faraday) (10000 coulomb)
    R = (k-mole) (joule/degC)    
    PI = (pi) (1)
    (mol) = (1)

}

PARAMETER {

    DNa = 0.15 (um2/ms)

    k1 = 1.0 (/mM3-ms)

    k2 = 0.001 (/ms)
    
    k3 = 0.3 (/ms)
                            : to eliminate pump, set TotalPump to 0 in hoc
    TotalPump = 1e-14 (mol/cm2)
    
    ik_ratio = -0.66666666 (1)

}

ASSIGNED {
    diam (um)
    L (um)
    ina (mA/cm2)
    nai (mM)
    vrat[Nannuli] : numeric value of vrat[i] equals the volume
                 : of annulus i of a 1um diameter cylinder
                 : multiply by diam^2 to get volume per um length
    
    k4          (/mM3-ms)

    nao (mM)
    ena (mV)

    ina_pmp (mA/cm2)
    parea (um)

    ik_pump (mA/cm2)



}

CONSTANT { volo = 1e10 (um2) }

STATE {
    : na[0] is equivalent to nai
    na[Nannuli] (mM) <1e-3>

    pump (mol/cm2)
    pumpna (mol/cm2)

}


BREAKPOINT {


    SOLVE state METHOD sparse

    ina = ina_pmp
    ik_pump = ik_ratio*ina_pmp

}

LOCAL factors_done

INITIAL {
    nai=15 :by Sara
    k4=(((nai/nao)^3)*k1*k3)/k2    :Set the equilibrium at nai0_na_ion
    parea = PI*diam
    pump = TotalPump/(1 + (nai*k1/k2))
    pumpna = TotalPump - pump
    if (factors_done == 0) {    : flag becomes 1 in the first segment
        factors_done = 1        : all subsequent segments will have
        factors()               : vrat = 0 unless vrat is GLOBAL
    }

    FROM i=0 TO Nannuli-1 {
        na[i] = nai

    }

}

LOCAL frat[Nannuli]     : scales the rate constants for model geometry

PROCEDURE factors() {
    LOCAL r, dr2
    r = 1/2                 : starts at edge (half diam)
    dr2 = r/(Nannuli-1)/2   : full thickness of outermost annulus,
                            : half thickness of all other annuli
    vrat[0] = 0
    frat[0] = 2*r
    FROM i=0 TO Nannuli-2 {
        vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
        r = r - dr2
        frat[i+1] = 2*PI*r/(2*dr2)              : outer radius of annulus
                                                : div by distance between centers
        r = r - dr2
        vrat[i+1] = PI*(r+dr2/2)*2*dr2 : outer half of annulus
    }
}


KINETIC state {
    COMPARTMENT i, diam*diam*vrat[i] {na}
    COMPARTMENT (1e10)*parea {pump pumpna}
    COMPARTMENT volo {nao}

    LONGITUDINAL_DIFFUSION i, DNa*diam*diam*vrat[i] {na}

    :pump
    ~ 3 na[0] + pump <-> pumpna (k1*parea*(1e10), k2*parea*(1e10))
    ~ pumpna <-> pump + 3 nao (k3*parea*(1e10), k4*parea*(1e10))

    CONSERVE pump + pumpna = TotalPump * parea * (1e10)
    ina_pmp = FARADAY*(f_flux - b_flux)/parea

    : all currents except pump
    ~ na[0] << (-(ina-ina_pmp)*PI*diam/(FARADAY)) : ina is Na efflux

    FROM i=0 TO Nannuli-2 {
        ~ na[i] <-> na[i+1] (DNa*frat[i+1], DNa*frat[i+1])
    }
    nai = na[0]
    ena = ((R*(273.15+celsius))/(FARADAY*10))*log(nao/nai)
}

