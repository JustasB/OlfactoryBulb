// Granule cell model of Migliore and Shepherd 
// keep only primary dendrite
// Translated by Aditya NCBS 2010.
// fields are ..
// name
// parent
// x,y,z  // coords of endpoint. 
// dia // needed for memb props. All lengths are in microns.
// ch name density 
// ch name density 
// .....
// Control lines start with '*'. Valid control options are 
// *relative 			- relative coords.
// *absolute			- absolute coords.
// *asymmetric			- use asymmetric compartments
// *symmetric			- use symmetric compartments

// #	name	parent		x	y	z	d	ch	dens	ch	dens	.	.	.

*symmetric //Neuron has only symmetric compartments - however hsolve symmetrizes the Ra-s internally.
// With the present implementation, it doesn't matter whether you write symmetric or asymmetric here.
*relative

*set_global	EREST_ACT	-0.065

// reduced two compartmental granule adapted from Migliore and Shepherd

*cartesian

*start_cell
// Both Migliore and Shepherd 2008 and Egger et al 2003 say taum = 30ms, Rin=1GOhm. So set RM and CM for it. 
*set_global	RA	0.8 // Ohm-m from 80 Ohm-cm
*set_global	CM	0.04 // F/m^2
//*set_global	RM	0.75 // RM = 30e-3/CM = 0.75 Ohm-m^2
// Aditya: I reduced CM to 0.02 from 0.04. Thus time constant is now 15ms not 30ms.
*set_global	RM	0.75 // RM = 30e-3/CM = 0.75 Ohm-m^2
// soma is somagc of gc.hoc of Migliore and Shepherd 2008.
soma	none	8	0	0	8	Na_rat_ms   400 KDR_ms	60  KA_ms   40

// periphery is priden of gc.hoc of Migliore and Shepherd 2008.
// Following Migliore and Shepherd 2008, I have only KA and not KDR in dendrites! This is supposed to cause spike latency
periphery  soma	    150	0	0	0.5 Na_rat_ms   200 KA_ms   80
