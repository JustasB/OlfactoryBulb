//	PARAMETER FILE FOR NEURON 'block' : mitral cell model under
//		conditions of TEA and TTX block
//	Author : Upi Bhalla 
//	Mar 29 1991
//	Highly detailed model of mit cell with experimental averages for
//		cell geometry.


//	Format of file :
// x,y,z,dia are in microns, all other units are SI (Meter Kilogram Second Amp)
// In polar mode 'r' is in microns, theta and phi in degrees 
// Control line options start with a '*'
// The format for each compartment parameter line is :
//name	parent	r	theta	phi	d	ch	dens ...
//in polar mode, and in cartesian mode :
//name	parent	x	y	z	d	ch	dens ...

//		Coordinate mode
*cartesian
*relative

//		Specifying constants
*set_global	EREST_ACT	-0.065

///// Aditya translated Davison's Bhalla and Bower adapted 7 compartment model from Neuron into Genesis.
///// Davison's model has pseudo compartments interleaved within proper compartments soma, prim, glom and dend.
///// pseudo compartments have only Ra and negligible Cm and Rm.
//        #print self.Rs2d # 5154639 Ohms
//        #print self.Rs2p # 18281535.6 Ohms
//        #print self.Rp2g # 17064846 Ohms
// Ra = RA l / (pi r^2) = RA / (pi r^2), Rm = RM / (2 pi r l), Cm = CM * (2 pi r l) 
///// proper compartments have only Rm and Cm, negligible Ra.
///// But cannot set Ra, Rm, Cm values to less than 1e-15, hence CM = 1e-3. However, 1e-15 leads to numerical error accumulation so set reasonable values
///// This is done to ensure no difference between symmetric and asymmetric compartments.

*set_global	RA	1e-4
*set_global	CM	0.01
*set_global	RM	10
//soma	none	100	0	0	16	LCa3_mit_usb	40	K_mit_usb	28	KA_bsg_yka	58.7	Ca_mit_conc	5.2e-6	Kca_mit_usb	142	Na_mit_usb	1532	K2_mit_usb	1956
//LCa3_mit_usb causes a 15ms current transient with voltage clamp in mitral cell (Ca ionic channels might cause part of observed 'IPSC'.) - fig 2A left of Isaacson & Westbrook 1998.
//To remove double spikes, the fast first one of which does not propagate into the soma, can increase KA to 587 or LCa3 to 400.
soma	none	100	0	0	16	LCa3_mit_usb	40	K_mit_usb	28	KA_bsg_yka	58.7	Ca_mit_conc	5.2e-6	Kca_mit_usb	142	Na_mit_usb	1532	K2_mit_usb	1956

*set_global	RM	1e5
*set_global CM  1e-3
// set such that Rs2p # 18281535.6 Ohms
*set_global RA  57.433089
s2p	    soma	1	0	0	2

*set_global	RA	1e-4
*set_global	CM	0.01
*set_global	RM	10
//prim    s2p 	100	0	0	104	LCa3_mit_usb	22	K_mit_usb	17.4	K2_mit_usb	12.3	Na_mit_usb	13.4
// fast and full propagation along primary dendrite:
prim    s2p 	100	0	0	104	LCa3_mit_usb	10	K_mit_usb	17.4	K2_mit_usb	123	    Na_mit_usb	134

*set_global	RM	1e5
*set_global CM  1e-3
// set such that Rs2p # 17064846 Ohms
*set_global RA  53.6107495
p2g     prim    1   0   0   2

*set_global	RA	1e-4
*set_global	CM	0.01
*set_global	RM	10
// original bhalla and bower model had K_mit_usb as 28, but Davison used 200 as he says glom is too excitable. I also use 200.
//glom	p2g 	100	0	0	26.7	LCa3_mit_usb	95	K_mit_usb	28
glom	p2g 	100	0	0	26.7	LCa3_mit_usb	95	K_mit_usb	200
// I put below values to make mitral fire more and so that glom is less excitable - however it makes mitral double spike methinks.
//glom	p2g 	100	0	0	26.7	LCa3_mit_usb	20	K2_mit_usb	123 Na_mit_usb	134

*set_global	RM	1e5
*set_global CM  1e-3
// set such that Rs2p # 5154639 Ohms - use this Ra and the dia and len of dend below to get RA. Don't use this pseudo segment's dia and len
//*set_global RA  16.193762
// RA was too large, the lateral dendrite was not responding fast enough - too long time constant. See Margie et al 2000.
// Break lateral dendrite into two compartments
*set_global RA  16
s2d     soma    1   0   0   2

// after splitting into two, the electrotonic lengths have changed now.

*set_global	RA	1e-4
*set_global	CM	0.01
*set_global	RM	10
// Using Ra above and l and dia for dend below to get RA + using RM and dia of dend, lambda = 600microns. So electrotonic L of dend = 1/6.
// Compare with lambda of Upi's cell which is 1600 microns, and typical dendrite length is 850microns, hence L=0.5.
// Unlike BBmit1993, Davison doesn't have these LCa3_mit_usb and K_mit_usb in the dendrite!
//dend	s2d		100	0	0	170.9	LCa3_mit_usb	4	K_mit_usb	8.5	K2_mit_usb	226	Na_mit_usb	330
dend	s2d		100	0	0	170.9	LCa3_mit_usb	4	K_mit_usb	8.5	K2_mit_usb	226	Na_mit_usb	330

*set_global	RM	1e5
*set_global CM  1e-3
// set such that Rs2p # 5154639 Ohms - use this Ra and the dia and len of dend below to get RA. Don't use this pseudo segment's dia and len
//*set_global RA  16.193762
// RA was too large, the lateral dendrite was not responding fast enough - too long time constant. See Margie et al 2000.
// Break lateral dendrite into two compartments
*set_global RA  8
//s2d2     dend    1   0   0   2

*set_global	RA	1e-4
*set_global	CM	0.01
*set_global	RM	10
// Unlike BBmit1993, Davison doesn't have these LCa3_mit_usb and K_mit_usb in the dendrite!
//dend2	s2d2	50	0	0	170.9	LCa3_mit_usb	4	K_mit_usb	8.5	K2_mit_usb	226	Na_mit_usb	330
