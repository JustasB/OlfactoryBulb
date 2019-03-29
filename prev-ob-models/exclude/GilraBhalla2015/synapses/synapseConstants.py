from synapseConstantsMinimal import *

############# ORN to mitral!
############# Set these such that 100 ORNs at approx 50Hz make the mitral cell fire in the middle of its linear range.
#SYN_EXC_G = 1 * 8.6516e-9 # Siemens
#SYN_INH_G = 1 * 2.2126e-9 # Siemens

GRANULE_INH_GRADED = False#True

RECEPTOR_SATURATION = 1.0#0.5 # Needed for single synapse KinSynChan and also in the baseline SynChan to correct for usage of synaptic weights across synchan and kinsynchan.
#RECEPTOR_SATN_CORRECTN_NMDA = 0.8785 # above 6mV EPSP #0.22 #for CM=0.04 #0.25 #for CM=0.01 #0.275 # See my onenote notes dt03/07/2010 for derivation of this value.
#RECEPTOR_SATN_CORRECTN_AMPA = 1.0818 # above 6mV EPSP #0.35 #0.19 # See my onenote notes dt03/07/2010 for derivation of this value.
RECEPTOR_SATN_CORRECTN_NMDA = 0.892 #4mV EPSP #0.22 #for CM=0.04 #0.25 #for CM=0.01 #0.275 # See my onenote notes dt03/07/2010 for derivation of this value.
RECEPTOR_SATN_CORRECTN_AMPA = 1.069 #4mV EPSP #0.35 #0.19 # See my onenote notes dt03/07/2010 for derivation of this value.
#RECEPTOR_SATN_CORRECTN_NMDA = 0.9397 #7mV EPSP with CM=0.02 #0.22 #for CM=0.04 #0.25 #for CM=0.01 #0.275 # See my onenote notes dt03/07/2010 for derivation of this value.
#RECEPTOR_SATN_CORRECTN_AMPA = 1.0679 #7mV EPSP with CM=0.02 #0.35 #0.19 # See my onenote notes dt03/07/2010 for derivation of this value.

## Arevian et al's activity dep inhibition work is under physiological conditions 
## i.e. 1mM Mg. as also Urban's previous work.
MG_CONC = 1.0 #0.001 #1.3 #mM SI units mol/m^3 = mmol/liter = mMolar (mM) # value of 1.3 mM From Isaacson 2001 PNAS; [Mg++] should be non-zero, hence 0.001 for 0.0
## Giridhar et al use no Mg2+, hence this setting for testing asymmetrical lateral inhibition
#MG_CONC = 0.1 #0.001 #1.3 #mM SI units mol/m^3 = mmol/liter = mMolar (mM) # value of 1.3 mM From Isaacson 2001 PNAS; [Mg++] should be non-zero, hence 0.001 for 0.0

mitral_granule_AMPA_Ek = 0.0 # Volts
## below is imported from synapseConstantsMinimal.py
#mitral_granule_AMPA_Gbar = 0.35e-9#0.35e-9 # Siemens ## This has been set so as to get roughly 8mV near the 12mV EPSP of Trombley & Shepherd 1992 JNeurosci
mitral_granule_saturatingAMPA_Gbar = 0.35e-9/RECEPTOR_SATN_CORRECTN_AMPA# Siemens ## This has been set so as to get roughly 8mV near the 12mV EPSP of Trombley & Shepherd 1992 JNeurosci
# From Cang and Isaacson 2003, in-vivo whole cell sEPSP data: 1ms rise and 4ms decay.
mitral_granule_AMPA_tau1 = 1.0e-3#2.0e-3 # seconds # Davison etal 2003 assume instantaneous rise but write that 4 ms is experimental, but is this for AMPA or GABA!
mitral_granule_AMPA_tau2 = 4.0e-3#5.5e-3 # seconds # decay time - from Migliore and Shepherd 2008.
mitral_granule_saturatingAMPA_pulseWidth = mitral_granule_AMPA_tau1 # this is the time to peak for KinSynChan
mitral_granule_saturatingAMPA_tau1 = mitral_granule_AMPA_tau2 # saturating syn decay time - roughly from fig 1B of Migliore and Shepherd 2008.
mitral_granule_saturatingAMPA_rInf = RECEPTOR_SATURATION # make some receptors saturate - set it so that it is best for lateral inhibition.

mitral_granule_NMDA_Ek = 0.0 # Volts
mitral_granule_NMDA_Gbar = 0.26*mitral_granule_AMPA_Gbar #0.1e-9 # Siemens ## This has been set so as to get roughly 8mV near the 12mV EPSP of Trombley & Shepherd 1992 JNeurosci
mitral_granule_saturatingNMDA_Gbar = 0.26*mitral_granule_AMPA_Gbar/RECEPTOR_SATN_CORRECTN_NMDA # Siemens ## This has been set so as to get roughly 8mV near the 12mV EPSP of Trombley & Shepherd 1992 JNeurosci
##### For SynChan
mitral_granule_NMDA_tau1 = 25e-3#20e-3 # rise time
mitral_granule_NMDA_tau2 = 200e-3#50e-3 # decay time - roughly from Migliore and Shepherd 2008.
##### For KinSynChan - rinf - fraction open due to a transmitter release and tau1 - decay time
mitral_granule_saturatingNMDA_pulseWidth = mitral_granule_NMDA_tau1 # this is the time to peak for KinSynChan
mitral_granule_saturatingNMDA_tau1 = mitral_granule_NMDA_tau2 # decay time - roughly from Migliore and Shepherd 2008.
mitral_granule_saturatingNMDA_rInf = RECEPTOR_SATURATION # make some receptors saturate - first pass set it so that it is best for lateral inhibition.
mitral_granule_NMDA_KMg_A = 1.0/0.1 #1.0/0.33 # mM
mitral_granule_NMDA_KMg_B = 1.0/73.0 #1.0/60.0 # V
mitral_granule_NMDA_MgConc = MG_CONC

## Choose between using a short-term plastic synapse or a non-plastic synapse.
GABA_plastic = False
GABA_depression_factor = 0.8 # should be 0.5 from Murthy's 2005 paper but aggregated synapses (?)
GABA_recovery_time = 6.0 # seconds # From Venki Murthy's 2005 paper.

granule_mitral_GABA_Ek = -0.078 # Volts
#### averaged inhibitory synapse:
#granule_mitral_GABA_Gbar = 0.6e-9 # Siemens
#granule_mitral_GABA_tau1 = 50e-3 # averaged IPSP from Urban and Sakmann 2002 #1e-3 # seconds # Davison etal 2003 assume instantaneous rise but write that 4 ms is experimental, but is this for AMPA or GABA!
#granule_mitral_GABA_tau2 = 75e-3 # averaged IPSP from Urban and Saknann 2002 #10e-3 # seconds # roughly from fig 2A of Isaacson and Strowbridge 1998.
#### unitary inhibitory synapse:
## below are imported from synapseConstantsMinimal.py
#granule_mitral_GABA_Gbar = 10e-9#15e-9#2e-9 # Siemens
#self_mitral_GABA_Gbar = 50e-12 # Siemens
granule_mitral_GABA_tau1 = 1e-3#3e-3 # roughly from fig 1C of Schoppa et al 1998
granule_mitral_GABA_tau2 = 20e-3#1e-3#20e-3 # from text and fig 1C of Schoppa et al 1998 
