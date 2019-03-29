RM = 12.0 # ohm m^2 from granule.tem
CM = 0.01 # F/m^2 from Bhalla and Bower 1993 ##### cannot find CM being set in granule.tem!!! though this is perhaps the NEURON default
RA = 0.5 # ohm m from Bhalla and Bower 1993 #### Cannot find RA being set in granule.tem!!!!! I don't think this is the NEURON default
#### actually RA is not used. psuedo compartments are used to model the inter-compartmental resistances.

EREST = -0.065 # Volts
sarea = 5e-9 # m^2 default surface area of soma

SURFACE_AREA_TOTAL	= 8353.0e-12	# m^2
LENGTH		= 50.0e-6		# m
SOMA_SAREA_FRAC	= 0.0136
PERIPH_SAREA_FRAC = 0.308
DEEP_SAREA_FRAC = 1 - SOMA_SAREA_FRAC - PERIPH_SAREA_FRAC
G_SOMA_PERI 	= 3.08e-6	# S.m^-2 from granule.tem
G_SOMA_DEEP 	= 4.34e-6    	# S.m^-2 from granule.tem

SOMASTIM = True
PERISTIM = False

IFULL = 0.0125e-9 # nA ### without electrode leak as here, need to halve the current to match Bhalla and Bower as the input resistance is higher.

