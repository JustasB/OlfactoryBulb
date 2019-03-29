/***********************************************************************************************

  IF_NET.. (project with Eleni on STDP + freq-dependent synapses)

SIMULATION CODE
***********************************************************************************************/

void simulate(double Tsim) {
 INT i,j,k,center=0;
 double tstop;

 printf("Simulation started - %.0f msec to go.\n", Tsim);
 tstop = t + Tsim;
 
  while (t < tstop) {      // Main simulation cycle.. (i.e t=0,dt,2dt,...)

  center = center + (fmod(t, 5.)<dt);
  center = center % N; 
  Nspiked = 0;
  for (i=0; i<N; i++) { //------------- LET's EVOLVE ALL THE STATE VARIABLES --------------------------------
      Iext[i] = Iext[i] * (1. - dt/taui) + dt*mu[i]/taui + sigma * sqrt(2 * dt/taui) * gauss();
   //Iext[i] = Iext[i] * tmpnoise1 + tmpnoise2 + tmpnoise3 * gauss();   
   // Evolution of neuronal variables (unless refractory) - one for each cell
   //udot = (t - to[i]) < Tarp ? 0 : (dt / C) * (  g_L * (E_L - u[i]) + g_L * Delta_T * exp((u[i] - V_T)/Delta_T) - w[i] + Ge[i] + Iext[i] );
   udot = (t - to[i]) < Tarp ? 0 : tmp_u1 - tmp_u2 * u[i] + tmp_u3 * exp((u[i] - V_T)/Delta_T)  + tmp_u4 * (-w[i] + Ge[i] + Iext[i]);
      
   //if ((i==center) & (fmod(t, 5.)<dt)) u[i] = th;
   if ((abs(i-center)<=0) & (fmod(t, 5.)<dt)) u[i] = th;
      
   if ((u[i] + udot) > th) {
   udot = th - u[i];
   spiked[Nspiked++] = i;
   rate++;
   rates[i]++;    
   //printf("spike in unit %d\n", (int) i);
   }
   wdot = (dt/tau_w) * (a * (u[i] - E_L) - w[i]);	
   w[i] += wdot;
   u[i] += udot;
 
   // Evolution of synaptic currents and STDP-associated variables (one for each cell)
   Ge[i] *= tmp_Ge;
   r1[i] *= tmp_r1;
   r2[i] *= tmp_r2;
   o1[i] *= tmp_o1;
   o2[i] *= tmp_o2;
  } // end for i --------------------------------------------------------------------------------------------

      print_output();
  
 if (Nspiked > 0) { //------------- WHO DID JUST FIRED ?!?? ------------------------------------------------
  for (k=0; k<Nspiked; k++) { 	// I examine only those neurons that actually just fired
  i         = spiked[k];		// For each of them: i.e. neuron i-th just fired!!!
  u[i]      = H;				// Its membrane potential is reset to an hyperpolarized value.
  w[i]     += b;				// Its adaptation mechanisms are updated accordingly.
  delta     = t - to[i];		// or = t-dt - to[i] ???	Interspike interval is computed.
  to_old[i] = to[i];			// perhaps useless
  to[i]     = t;				// or = t - dt; ????		Its firing time is recorded.
  r1[i]++;						// Its STDP presynaptic activity sensor is updated.
  r2[i]++;						// Its STDP presynaptic activity sensor is updated.
  o1[i]++;						// Its STDP postsynaptic activity sensor is updated.
  o2[i]++;						// Its STDP postsynaptic activity sensor is updated.
  fprintf(raster,"%f %ld\n", t, i);// The firing event is stored on disk.
  
  for (j=0; j<N; j++) {//------------- MANAGE SYNAPTIC INTERACTIONS AND SHORT/LONG-TERM PLASTICITY -----------
   if (CC[j][i]!=0)  { 	// Cycle over all POSTynaptic neurons that the neuron i-th targets...
    Ge[j]   += W[j][i] * GA * Gu[j][i] * Gr[j][i];	// I manage the consequence of a presynaptic spike (i.e. EPSC).
    U        = (CC[j][i] > 0) ? fU : dU;
    tau_rec  = (CC[j][i] > 0 ) ? fD : dD;
    tau_fac  = (CC[j][i] > 0) ? fF : dF;
    tmp_rec  = exp(-delta/tau_rec);
    tmp_fac  = exp(-delta/tau_fac);
         
    //printf("C[%d, %d] = %d\n", (int)j,(int)i,CC[j][i]);
    
    Gr[j][i] = Gr[j][i] * (1. - Gu[j][i]) * tmp_rec + 1 - tmp_rec;	// Short-term synaptic depression.
    Gu[j][i] = Gu[j][i] * (1. - U) * tmp_fac + U;					// Short-term synaptic facilitation.
   
    W[j][i] -= eta * o1[j] * (A2m + A3m * (r2[i]-1.));				// Long-term depression event.
    W[j][i] = (W[j][i] < 0) ? 0. : W[j][i];							// Clipping to the lower bound 0.
   } // if CC[j,i] != 0
  
   if (CC[i][j]!=0)  { 	// Cycle over all (presynaptic) neurons that target neuron i-th...
    W[i][j] += eta * r1[j] * (A2p + A3p * (o2[i]-1.));				// Long-term potentiation event.
    W[i][j] = (W[i][j] > maxweight) ? maxweight : W[i][j];			// Clipping to the upper bound maxw..
   } // if CC[j][i] != 0
   } // end for j
  } // end for k
 } // end if Nspiked > 0 ----------------------------------------------------------------------------------
 

/*
if (fmod(t,Tsim * 0.25)<=dt) {
   printf("25%% of %f - has been simulated\r", Tsim); fflush(NULL); }
  else if (fmod(t,Tsim * 0.5)<=dt) {
   printf("50%% of %f - has been simulated\r", Tsim); fflush(NULL); }
  else if (fmod(t,Tsim * 0.75)<=dt) {
   printf("75%% of %f - has been simulated\r", Tsim); fflush(NULL); }
*/
  fflush(NULL);
  
  t += dt;                  // Finally, the current simulation time is increased by dt.
  } // end while()          // End of the main simulation loop.
printf("Simulation done!\n\n");
}// end simulate()


        
