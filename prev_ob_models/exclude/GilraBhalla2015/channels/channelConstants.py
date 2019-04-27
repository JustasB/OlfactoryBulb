GAS_CONST = 8.314 # J/K/mol
FARADAY = 9.6485e4 # C/mol

VMIN = -0.1 # Volts
VMAX = 0.1 # Volts
NDIVS = 200 # number
dv = ( VMAX - VMIN ) / NDIVS # Volts

# Arevian et al's activity dep inhibition work is under physiological conditions i.e. 1mM Mg. as also Urban's previous work.
MG_CONC = 1.0 #0.001 #1.3 #mM SI units mol/m^3 = mmol/liter = mMolar (mM) # value of 1.3 mM From Isaacson 2001 PNAS; [Mg++] should be non-zero, hence 0.001 for 0.0

CELSIUS = 35.0 # deg Celsius

##### Not a global constant but put in to set default values of conductances in the channels
sarea = 5e-9 # m^2 typical surface area of a mitral cell soma
