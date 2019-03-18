/***********************************************************************************************

  IF_NET.. (project with Eleni on STDP + freq-dependent synapses)

DATA STRUCTURES DECLARATION
***********************************************************************************************/

typedef unsigned long int INT;  // Definition of a new LONG type, just for convenience.


FILE *fopen();                  // File pointers to be used for logging results on disk.
FILE *sample;		            // Output file, contains sample membrane voltage..
FILE *raster;  		            // Output file, contains firing times...
FILE *gfile;				    // Output file, contains weight or conductance time series...
FILE *weightmatrix;		        // Output file, contains the actual weight matrix...
FILE *frate;                    // Output file, contains the mean firing rates for each neuron...

INT N;                          // Overall number of excitatory neurons (set by the MATLAB-GUI)..


INT i, j, k;                    // Temporary useful indexes..
double tau_rec, tau_fac;		// Temporary variables.
double U, tmp_rec, tmp_fac, delta;// Temporary variables.
double udot, wdot;				// Temporary variables.
double tmp_r1, tmp_r2;			// Numerical efficiency temporary variables (i.e. calculated once for all)
double tmp_Ge, tmp_o1, tmp_o2;	// Numerical efficiency temporary variables (i.e. calculated once for all)
double tmp_u1, tmp_u2, tmp_u3;	// Numerical efficiency temporary variables (i.e. calculated once for all)
double tmp_u4;// Numerical efficiency temporary variables (i.e. calculated once for all)
double tmpnoise1, tmpnoise2, tmpnoise3; // Numerical efficiency temporary variables (i.e. calculated once for all)

double seed, fixseed;			// Imposed random seed and flag {0,1} on whether to fix the seed to it.

double t;						// Actual simulation time               [ms]
double dt, sdt;					// Integration time step and its square root [ms, ms^-1] (Forward Euler)
double Tsim;					    // Overall simulation time..			[ms]

double Io;                      // variables for the single-cell current injection(s)..
double Itot;					// (i.e. maximal amplitude, delays, actual waveform)
double off, del, period, number; //
double *mu;						// Mean of the input current - accounting for heterogeneity
double sigma;					// Std.dev. of the input current
double taui;					// Autocorrelation length of the input current fluctuations


short **CC;               // Dynamic/Static connectivity matrix: 0, iff no connection, 1, iff facilitating, -1 iff depressing
double **W;						// Synaptic weight matrix, modified by the STDP triplet rule...
double maxweight;				// Upper bound for the elements of the synaptic matrix W (ie. lower is 0)
double fU, fD, fF;				// Synaptic parameters (U, tau_D, tau_F) for facilitating synapses
double dU, dD, dF;				// Synaptic parameters (U, tau_D, tau_F) for depressing synapses
double **Gr, **Gu;              // Synaptic dynamics variables, one per synapse (since the same neuron can establish different synapses with distinct targets)
double *r1, *r2, *o1, *o2;		// Long-term, triplet-based STDP rule, pre- and postsynaptic indicators
double taur2, tauo2, taur1, tauo1; // Decay time constant for pre and postsynaptic indicators..
double A2p, A3p, A2m, A3m;		// Amplitude of weight change in the various conditions..
double eta;						// Learning rate for STDP weight changes...

INT *spiked;					// Vector containing at runtime the indexes of those neurons who fired.
INT Nspiked;					// Counter for counting how many neurons fired.
double rate;                     // Variable containing the mean firing rate of the entire network.
double *rates;                   // Vars containing, for each neurons, its mean firing rate.
double Trate;                    // Variable containing the time window 'rate' and 'rates[]' are computed over.

double *u;                      // Membrane voltage       [mV]
double *w;						// Spike-frequency adaptation variable [pA]
double *to;      // Last firing time for the i-th neuron [ms] (needed for the Tarp)
double *to_old;  // Last firing time for the i-th neuron [ms] (needed for the synaptic dynamics) // maybe useless


double *Ge;						// Incoming current to each neuron
double *Iext;                   // External current to each neuron
double Getot;					// maybe useless


double th;						// Peak value for a spike [mV]
double C;						// Membrane capacitance [pF]
double g_L;                     // Membrane conductance [nS]
double E_L;						// Resting voltage        [mV]
double V_T;						// Excitability threshold [mV]
double H;                       // Reset voltage          [mV]
double Tarp;					// Absolute refractory period  [ms]
double Delta_T;					// Excitability slope     [mV]

double tau_w, a, b;				// spike-frequency adaptation

double latency, tau_Isyn, GA;	// Synaptic parameters..

        