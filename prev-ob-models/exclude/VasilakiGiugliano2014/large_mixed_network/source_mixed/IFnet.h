/**********************  INIT  *************************************/
void init() {
 INT i,j,k;

  dt    = 0.1;						// Integration time step [ms].
  sdt   = sqrt(0.1);				// Square root of the dt [ms^0.5].

  eta   = 1.;						// Learning rate

  // (1) Model neuron parameters (i.e. exponential I&F)
  th    = 20.;						// [mV] - peak value for a spike
  C     = 281.;     					// Membrane capacitance [pF]
  g_L   = 30.;						// Membrane conductance [nS] (i.e. taum = 16.6667 ms) 
  E_L   = -70.6;						// Resting membrane potential [mV]
  V_T   = -50.4;                    // Excitability threshold [mV].
  H     = -70.6;                    // Reset potential [mV].
  Delta_T = 2.;						// Excitability slope [mV]
  Tarp  = 2.;						// Absolute refractory period [ms].

 // (2) Model neuron parameters (i.e. spike-freq. adaptation of the exp I&F)
  tau_w   = 144.;					// [ms] - decay time constant for adaptation variable
  a       = 4.;						// [nS] - voltage-dependence of adaptation variable
  b       = 0.0805;					// [nA] - spike-dependence of adaptation variable

 // (3) Synaptic coupling (and plasticity) parameters
  maxweight	= 5.;					// was 1.2
  tau_Isyn= 5.;						// Decay time constant of individual EPSC [ms]
  //GA    = 200.;						// Maximal synaptic weight (for each synapse, same value) - was 200 for facilitating and depressing synapses
  fU    = 0.1;						// Probability of release (for each synapse, same value) 
  fD    = 100.;	// was 5					// [ms] recovery from short-term depression (for each synapse, same value)
  fF    = 900.; // was 5000					// [ms] recovery from short-term facilitation (for each synapse, same value)
 
  dU    = 0.8;						// Probability of release (for each synapse, same value) 
  dD    = 900.;	// was 500					// [ms] recovery from short-term depression (for each synapse, same value)
  dF    = 100.;	// was 5					// [ms] recovery from short-term facilitation (for each synapse, same value)

  rate  = 0.;  
  Trate = 5000.;                     // [ms] temporal window to compute mean firing frequencies over
    
 // -- VISUAL CORTEX, MINIMAL MODEL (Pfister & Gerstner, 2006) ------------------------------
  taur2    = 101.;           // Decay time constant of presynaptic indicators (tau_x)
  tauo2    = 114.;           // Decay time constant of postsynaptic indicators (tau_y)
  taur1    = 16.8;           // Decay time constant of presynaptic indicators (tau^+)
  tauo1    = 33.7;           // Decay time constant of postsynaptic indicators (tau^-)
	
  A2p      = 0;				// Amplitude of weight change, in case of a pre-post
  A3p      = 6.5e-3;		// Amplitude of weight change (Triplet-term), in case of a pre-post
  A2m      = 7.1e-3;		// Amplitude of weight change, in case of a post-pre
  A3m      = 0;				// Amplitude of weight change (Triplet-term), in case of a post-pre
 //-------------------------------------------------------------------------------------------	
  tmp_Ge		=	exp(-dt/tau_Isyn);

  tmp_r1		=	exp(-dt/taur1);
  tmp_r2		=	exp(-dt/taur2);
  tmp_o1		=	exp(-dt/tauo1);
  tmp_o2		=	exp(-dt/tauo2);
  
  tmp_u1     =   (dt / C) * g_L * E_L;
  tmp_u2     =   (dt / C) * g_L;
  tmp_u3     =   (dt / C) * g_L * Delta_T;
  tmp_u4     =   (dt / C);
  
//  mu    = 0.;
//  sigma = 0.;
  taui  = 5.;

  tmpnoise1 = (1. - dt/taui);
  //tmpnoise2 = dt*mu/taui;
  tmpnoise3 = sigma * sqrt(2 * dt/taui);


for(i=0;i<N;i++) {
	to[i]    = -9999.;                  // Initialization of the last firing times (to -infty)
    to_old[i]= -99999.;                 // Initialization of the last firing times (to -infty)
    u[i]     = E_L +(V_T - E_L) * 0*drand49(); // Random initialization of the membrane voltages.
    w[i]     = 0.;						// Spike-frequency adaptation variable [pA]
    Ge[i]    = 0.;
    mu[i]    = Io * (1 + 1 * gauss());
    //mu[i]    = Io;
    Iext[i]  = mu[i];
	r1[i]    = 0.;						// Presynaptic indicator, one for each neuron
	r2[i]    = 0.;						// Presynaptic indicator, one for each neuron
	o1[i]    = 0.;						// Postsynaptic indicator, one for each neuron
	o2[i]    = 0.;						// Postsynaptic indicator, one for each neuron
    
    rates[i] = 0.;
    }
    
    for(i=0; i<N; i++)
     for(j=0; j<N; j++) {
     
         if ((i<(N/2)) & (j<(N/2)))
             W[i][j] = drand49() * maxweight;
         else if ((i>(N/2)) & (j>(N/2)))
             W[i][j] = drand49() * maxweight;
         else
             W[i][j] = 0.1* drand49() * maxweight;
       
     //    W[i][j] = drand49() * maxweight;
      
        if (CC[j][i] != 0) // The current neuron may receive an EPSC..
            if (CC[j][i] == 1)   // Facilitatory excitatory synapse
                    {
                     Gr[j][i] = 1.;
                     Gu[j][i] = fU;
                    }
            if (CC[j][i] == -1)   // Depressing excitatory synapse
                    {
                     Gr[j][i] = 1.;
                     Gu[j][i] = dU;
                    }
     }

 Nspiked    = 0;			// ???		
 return;
}// end init()



int read_connectivity() {		// UPDATED JULY 2011 !!!!!!!!!!!!!!!!!!!
	// ..reads 'connectivity.dat' (binary) 
	// short **CC; <-- global variable
	// INT N;      <-- global variable
	INT i,j,k,out,MM;					// Useful temporary variables (e.g. indexes, return values, etc.)
	FILE *fopen(), *fp;					// Needed for file I/O operations...
	char fname[20];						// This will contain the file name for connectivity...
	double Ncells, tmp;					// This will contain the number of cells defined in the file...
	
	sprintf(fname, "connectivity.dat");	// This is the file generated by the MATLAB GUI [FIXED NAME]
    //NOTE:
    // This file has been generated in MATLAB by the following statements:
    //MM    = length(find(C(:)~=0));  % Number of non-zero elements of C(i,j)
    //Ctemp = C';
    //Ctemp = Ctemp(:);   % A single vector, encoding 'C' row-wise
    //tmp   = find(Ctemp ~= 0); // Only the non-zero elements of C are considered to be written on disk
    //fname = 'connectivity.dat';
    //fp    = fopen(fname, 'wb');
    //fwrite(fp, Ncell_exc, 'double');  
    //fwrite(fp, MM, 'double');  
    //for kk=1:length(tmp)
    //    fwrite(fp, tmp(kk) * sign(Ctemp(tmp(kk))), 'double');
    //end
    //fclose(fp);
	//------------------------------------------------------------------------------------------- 
	printf("Opening connectivity file <%s> (generated by MATLAB-GUI)...\n", fname);
	fp = fopen(fname, "r");				// Try to open file 'fname' for reading...
	if (fp==NULL) {  printf("Error: unable to open the file..\n");  return -1; }
	
	out = fread(&Ncells, sizeof(double), 1, fp); 	// I expect the first element to be 'Ncells'..
	if (out!=1) {  printf("Error: corrupted or wrong format for the connectivity file - first element is NOT Ncells!\n");  return -1; }
	
	N  = (INT) Ncells;								// First element stored as the global variable 'N'.
	printf("Found connectivity information for %d cells!\n", (int) N);
	
	out = fread(&tmp, sizeof(double), 1, fp); 	// I expect the second element to be 'MM'..
	if (out!=1) {  printf("Error: corrupted or wrong format for the connectivity file! - second element is NOT MM\n");  return -1; }
	MM = (INT) tmp;
	
	//------------------------------------------------------------------------------------------- 
	printf("Allocating memory [CC and W] for this simulation..\n");
	CC  = calloc(N, sizeof(short*));			// Reserving memory for the N x N int connectivity matrix
	if (CC==NULL) {  printf("Error: unable to allocate memory [CC]..\n"); fflush(NULL); return -1; }

	W  = calloc(N, sizeof(double*));				// Reserving memory for the N x N double synaptic matrix
	if (W==NULL) {  printf("Error: unable to allocate memory {W]..\n"); fflush(NULL); return -1; }
	
	for (i=0; i<N; i++) {
		CC[i] = calloc(N, sizeof(short)); // Reserving memory for the N x N int connectivity matrix
		if (CC[i]==NULL) {  printf("Error: unable to allocate memory [CCi]..\n"); fflush(NULL); return -1; }  

		W[i] = calloc(N, sizeof(double)); // Reserving memory for the N x N double synaptic matrix
		if (W[i]==NULL) {  printf("Error: unable to allocate memory [Wi]..\n"); fflush(NULL); return -1; }  
		
		//out = fread(CC[i], sizeof(short), N, fp);
		//if (out!=N) {  printf("Error: corrupted or wrong format for the connectivity file!\n", out);  return -1; }  
	}

    for (k=0; k<MM; k++) {
     out = fread(&tmp, sizeof(double), 1, fp); 	// I expect the second element to be 'MM'..
	 if (out!=1) {  printf("Error: corrupted or wrong format for the connectivity file!\n"); fflush(NULL); return -1; }
      
        
      i = (INT) ((fabs(tmp)-1)/N);  
      j = (INT) fmod(fabs(tmp)-1,N); 
      CC[i][j] = (short) ((tmp>0) ? 1 : -1);
    }
        printf("\n");
	// PRINT THE CONNECTIVITY MATRIX FOR DEBUGGING PURPOUSES
	/*
	for (i=0;i<N;i++) {
	 for (j=0; j<N; j++)
	 	printf("%d ", CC[i][j]);
	 printf("\n");
	}
	*/
	printf("Connectivity file <%s> acquired successfully!\n\n\n", fname);
	fclose(fp);									// File is closed here.
	//------------------------------------------------------------------------------------------- 
	
	return 0;
} // end read_connectivity()


void print_output() {
    if (fmod(t, Trate)<dt) {
        rate = 1000. * rate / (N * Trate);
        printf("Mean firing rate: %.0f Hz\r", rate);
        frate = fopen("rates.x", "w");
        for (i=0;i<N;i++) {
            fprintf(frate, "%f\n", 1000. * rates[i] / Trate);
            rates[i] = 0.;
        }
        fclose(frate);
    }
    
    if (fmod(t, 10*Trate)<dt) {
        gfile  = fopen("G.x","w");
        for (i=0;i<N;i++) {
            for (j=0;j<N;j++)
                fprintf(gfile, "%f ", W[i][j]);
            fprintf(gfile, "\n");
        }
        fclose(gfile);
    }

/*    
    //  if(fmod(t,20*dt)<=dt)   // Every 20 simulation steps, the membrane voltage is wrote on disk.
    fprintf(sample,"%f %f %f %f %f %f %f\n",t, u[0], Ge[0], u[1], Ge[1], u[2], Ge[2]);
    //if(fmod(t,100*dt)<dt) fprintf(gfile,"%f %f %f %f %f %f %f\n",t, W[0][1], W[1][0], W[3][5], W[5][3], W[7][4], W[4][7]);
 */
    return;
}


void log_weights() {
	INT i,j,k,out,MM;					// Useful temporary variables (e.g. indexes, return values, etc.)
	FILE *fopen();
	char fname[20];						// This will contain the file name for connectivity...
	
	sprintf(fname, "W_%.0f.x", t);	

	weightmatrix = fopen(fname, "w");				 // Open output file for writing...
    for (i=0; i<N; i++) {
     for (j=0; j<N; j++) 
	fprintf(weightmatrix, "%f ", W[i][j]);
	fprintf(weightmatrix, "\n");
	}
return;
} // end log_weights()




/*
int read_pars() {		// UPDATED JULY 2011 !!!!!!!!!!!!!!!!!!!
	//
	// ..reads 'pars.dat' (binary)
	//
	//[Ncell_exc, Prob_ee, Prob_loop, Prob_ff, Prob_dd, Prob_fd, Prob_f, Prob_d, fixseed, seed];
	//
	// INT N;      <-- global
	
	FILE *fopen(), *fp;					// Needed for file I/O operations...
	char fname[20];						// This will contain the file name for parameters...
	INT i, out;							// Useful temporary variables (e.g. indexes, return values, etc.) 
	double par;
	
	sprintf(fname, "pars.dat"); // This is the file generated by the MATLAB GUI [FIXED NAME]
	// NOTE: 
	//This file has been created by the following MATLAB-statements
    //
	//fname = 'pars.dat';
	//fp    = fopen(fname, 'wb');
	//all_pars = [Ncell_exc, Prob_ee, Prob_loop, Prob_ff, Prob_dd, Prob_fd, Prob_f, Prob_d, fixseed, seed];
	//fwrite(fp, pars, 'double');  
	//fclose(fp);
	//
	//We aim at extracting only two elements: the 9th "fixseed" and the 10th "seed".

	
	//------------------------------------------------------------------------------------------- 
	printf("Opening parameter file <%s>..\n", fname);
	fp = fopen(fname, "r");				// Try to open file 'fname' for reading...
	if (fp==NULL) {  printf("Error: unable to open the file..\n");  return -1; }
	
	// The first 8 elementE_Ls are inessential at run-time for the simulation...
	for (i=0; i<8; i++) {							// Let's go through each of them...
		out = fread(&par, sizeof(double), 1, fp);	// and simply ignore it, while the file-access pointer moves forward...
		if (out!=1) {  printf("Error: unable to read connectivity file headers..\n");  return -1; }
	}
    
	out = fread(&par, sizeof(double), 1, fp);  fixseed = par;
	if (out!=1) {  printf("Error: unable to read connectivity file headers..\n");  return -1; }
	out = fread(&par, sizeof(double), 1, fp);  seed = par;
	if (out!=1) {  printf("Error: unable to read connectivity file headers..\n");  return -1; }
	
	printf("Parameter file <%s> acquired successfully!\n\n\n", fname);
	fclose(fp);									// File is closed here.
	//------------------------------------------------------------------------------------------- 
	
	return 0;
} // end read_pars()
*/

/*

void manage_synapses(INT i) {	// UPDATED JULY 2011 !!!!!!!!!!!!!!!!!!!
	// Neuron i-th is examined, ready to receive inputs from its presynaptic peers (j=1..N), if any...
	double delta;
	for(j=0;j<N;j++)     // Now, all the connected presynaptic neurons are checked for (recent) spiking activity.
	  if (CC[i][j] != 0) // A connection exists from neuron j-th (pre) to neuron i-th (post)!
         if ((t-(to[j]+latency))<dt && (t-(to[j]+latency))>=0) // Indeed, presynaptic neuron j-th fired - after latency of EPSPs...
            { 
            delta   = (to[j]-to_old[j]);
            if (CC[i][j] == 1)   // Facilitatory exciting synapse
                    {
                     Gr[i][j] = (Gr[i][j] * (1 - Gu[i][j]) - 1) * exp(-delta/fD) + 1;
                     Gu[i][j] = Gu[i][j] * (1 - fU) * exp(-delta/fF) + fU;
                    }
            else // (CC[i][j] == -1)   // Depressing exciting synapse
                    {
					Gr[i][j] = (Gr[i][j] * (1 - Gu[i][j]) - 1) * exp(-delta/dD) + 1;
					Gu[i][j] = Gu[i][j] * (1 - dU) * exp(-delta/dF) + dU;
                    }
			Ge[i]  += GA *Gr[j][i] * Gu[j][i];
            }
return;
} // end manage_synapses()




void manage_stdp(INT j) {		// UPDATED JULY 2011 !!!!!!!!!!!!!!!!!!!
   // Neuron j-th just fired...
   INT i;
   for(i=0;i<N;i++)    // Now, all the connected neurons, postsynaptic to the j-th must be updated..
        {
			W[i][j] -= eta * o1[i] * abs(CC[i][j]) * (A2m + A3m * r2[j]); // the j-th neuron is presynaptic to all (connected) generic neurons i-th
			W[j][i] += eta * r1[i] * abs(CC[j][i]) * (A2p + A3p * o2[j]); // the j-th neuron is postsynaptic to all (connected) generic neurons i-th
		} 
	
   r1[j]++;			// update the pre/post-synaptic detectors of the neuron j-th, which just fired
   r2[j]++;			// update the pre/post-synaptic detectors of the neuron j-th, which just fired 
   o1[j]++;			// update the pre/post-synaptic detectors of the neuron j-th, which just fired
   o2[j]++;			// update the pre/post-synaptic detectors of the neuron j-th, which just fired
	return;
} // end manage_stdp()
*/
