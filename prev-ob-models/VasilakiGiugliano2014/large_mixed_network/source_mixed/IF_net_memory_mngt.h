/***********************************************************************************************

  IF_NET.. (project with Eleni on STDP + freq-dependent synapses)

MEMORY MANAGEMENT
***********************************************************************************************/
//int allocate_mem()
//int deallocate_mem()

int allocate_mem() {
INT i, j, k;
int error = 0;

    printf("Allocating memory...\n"); fflush(NULL);
  raster = fopen("raster.x","w");                // Open output file.
  sample = fopen("sample.x","w");                // Open output file.
 
  to         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..
  if (to == NULL) error = 1;
  to_old     = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (to_old == NULL) error = 1;
  u          = calloc(N, sizeof(double));         // Dynamic allocation of data structures..
  if (u == NULL) error = 1;
  w          = calloc(N, sizeof(double));         // Dynamic allocation of data structures..
  if (w == NULL) error = 1;
  Ge         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..
  if (Ge == NULL) error = 1;
  Iext       = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (Iext == NULL) error = 1;
  r1         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (r1 == NULL) error = 1;
  r2         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (r2 == NULL) error = 1;
  o1         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (o1 == NULL) error = 1;
  o2         = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
  if (o2 == NULL) error = 1;
  spiked     = calloc(N, sizeof(INT));
  if (spiked == NULL) error = 1;
  rates         = calloc(N, sizeof(double));
  if (rates == NULL) error = 1;
  mu         = calloc(N, sizeof(double));
  if (mu == NULL) error = 1;
 
  
  Gr          = calloc(N, sizeof(double*));         // Dynamic allocation of data structures..  
  if (Gr == NULL) error = 1;
  Gu         = calloc(N, sizeof(double*));         // Dynamic allocation of data structures..      
  if (Gu == NULL) error = 1;
  for (i=0; i<N; i++) {
   Gr[i]       = calloc(N, sizeof(double));         // Dynamic allocation of data structures..  
   if (Gr[i] == NULL) error = 1;
   Gu[i]      = calloc(N, sizeof(double));         // Dynamic allocation of data structures..       
   if (Gu[i] == NULL) error = 1;
   }
   
   if (error) printf("Unable to allocate memory!\n\n"); 
   else printf("Done!\n\n");
return error;
} // end allocate_mem()




void deallocate_mem() {
INT i, j, k;

   printf("Deallocating the memory...\n");
   for(i=0;i<N;i++) {     // The allocated memory is set free.
    free(CC[i]); 
    free(W[i]); 
    free(Gr[i]); 
    free(Gu[i]); 
   }           
    
	free(Gr);  free(Gu);		// The allocated memory is set free.
	free(CC);  free(W);         // The allocated memory is set free.
	free(Ge); free(Iext);		// The allocated memory is set free.
	free(to); free(to_old);		// The allocated memory is set free.
	free(w); free(u);			// The allocated memory is set free.
	free(r1); free(r2); free(o1); free(o2);
	free(spiked); free(mu); free(rates);

    // The files that were previously opened are now closed.
    fclose(raster); fclose(sample);
	printf("Done!\n\n\n\n");
return;
} // end deallocate_mem()




        