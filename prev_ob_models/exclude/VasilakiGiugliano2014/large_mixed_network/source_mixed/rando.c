//
// Bern, 26/11/2001 - First used by M. Giugliano for the Mee_Duck *generate_trial* project
// Antwerp, Aug - Sep 2011, revised and cross-checked for Eleni's joint project. 
//
//

#include <math.h>
//
// double drand49(): returns a random floating-point uniformly distributed between 0.0 and 1.0        (from Numerical Recipes - Press et al.) 
// double srand49(long): initializes the seed and returns a uniformly distributed between 0.0 and 1.0 (from Numerical Recipes - Press et al.) 
// double gauss() : returns a normally distributed random number (i.e. zero mean, unitary variance and Gaussian distribution density) (from Numerical Recipes - Press et al.) 
//
// long mysrand49(long): initializes the seed and returns the previous one  <---- NEW ADD-ON FUNCTION
//

static long rand49_idum = -77531;

#define MM 714025
#define IA 1366
#define IC 150889

//----------------------------------------------------------------------------------------------------------------
double drand49() {
        static long iy,ir[98];
        static int iff=0;
        int j;

	if (rand49_idum < 0 || iff == 0) {
            iff=1;
            if((rand49_idum=(IC-rand49_idum) % MM)<0)
                             rand49_idum=(-rand49_idum);
            for (j=1;j<=97;j++) {
                    rand49_idum=(IA*(rand49_idum)+IC) % MM;
                    ir[j]=(rand49_idum);
            }
            rand49_idum=(IA*(rand49_idum)+IC) % MM;
            iy=(rand49_idum);
        }
        j=1 + 97.0*iy/MM;
	if (j > 97 || j < 1) printf("drand49(): This cannot happen.");
        iy=ir[j];
        rand49_idum=(IA*(rand49_idum)+IC) % MM;
        ir[j]=(rand49_idum);
        return (double) iy/MM;
} // end drand49()
//----------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------
double srand49(seed)
long seed;
{
   rand49_idum=(-seed);
   return drand49();
} // end srand49()
//----------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------
long mysrand49(seed)
long seed;
{
  long temp;
  temp = rand49_idum;
  rand49_idum = (-seed);
  return temp;
} // end mysrand49()
//----------------------------------------------------------------------------------------------------------------


#undef MM
#undef IA
#undef IC

//----------------------------------------------------------------------------------------------------------------
// Returns a normally distributed random number (i.e. zero mean, unitary variance and Gaussian distribution density
double gauss() {
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;
	
	if  (iset == 0) {
		do {
			v1=2.0*drand49()-1.0;
			v2=2.0*drand49()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
} // end gauss()
//----------------------------------------------------------------------------------------------------------------


