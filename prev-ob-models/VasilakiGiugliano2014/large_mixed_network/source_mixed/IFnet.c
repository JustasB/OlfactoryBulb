/***********************************************************************************************

  IF_NET.. (project with Eleni on STDP + freq-dependent synapses)
  Antwerpen, 29/7/2011 - 9/9/2011
    
  Compile with:    gcc -o IFnet IFnet.c -lm -O
  Output files:    sample.x, raster.x, G.x

************************************************************************************************/

#include <stdio.h>              // Included (needed by fprintf, fopen, etc. ).
#include <stdlib.h>             // Standard library (needed by atof(), atoi(), etc.).
#include <math.h>               // Mathematical library (needed by exp(), etc.).

#include "rando.c"                 // Random numbers generation custom library (~Numerical Recipes).
#include "IF_net_datastructures.h" // Data structures, variables, etc. 
#include "IFnet.h"				   // Utility functions.
#include "IF_net_simulate.h"	   // Actual simulation code.
#include "IF_net_memory_mngt.h"    // Allocation and deallocation routines.


/********************************  MAIN ***********************************/ 
int main (int argc,char *argv[]) { 
    int index, out;
 
  if (argc < 4)  {
      printf("USAGE: IFnet Tsim [sec]  Io [pA] GA [pA] (was 200.)\n"); 
      exit(0);
  }
    
  Tsim       = atof(argv[1]) * 1000.;           // convert it to [msec]
  Io         = atof(argv[2]);                   // [pA]  
  GA         = atof(argv[3]);
    
  printf("\nIFNet November 2011 - including STP and STDP\n\n");
    
  read_connectivity();  // read 'connectivity.dat' (including short-term dynamics information) and sets 'N' and 'CC'..
  out = allocate_mem();
    if (out) {
      fprintf(stderr, "Error in memory allocation! Aborting...\n");
        return 0;
    }
 //-------------------------------------------------------------------
  t = 0.;        // The current time is set to 0 ms.

     //mu = Io;
  sigma = 2*Io;
  init();          // All the parameters of the network are initialized. 
    
  GA         = atof(argv[3]);
 
 log_weights(); 
    
  for (index=0; index < 4; index++) {
   simulate(Tsim*0.25);      // The simulation is started.
   fflush(NULL);
//   log_weights();
  }
log_weights();
  //-------------------------------------------------------------------

  deallocate_mem();
  return 0;                                           // The software returns no error.
}//end main()
