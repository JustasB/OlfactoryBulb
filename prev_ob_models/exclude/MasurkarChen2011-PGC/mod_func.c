#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ICat1m3h_reg();
extern void _ICat2_reg();
extern void _IKt1m4h_reg();
extern void _IKt3_reg();
extern void _Ikt2m2h_reg();
extern void _ipulse2_reg();
extern void _nmitral_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ICat1m3h.mod");
fprintf(stderr," ICat2.mod");
fprintf(stderr," IKt1m4h.mod");
fprintf(stderr," IKt3.mod");
fprintf(stderr," Ikt2m2h.mod");
fprintf(stderr," ipulse2.mod");
fprintf(stderr," nmitral.mod");
fprintf(stderr, "\n");
    }
_ICat1m3h_reg();
_ICat2_reg();
_IKt1m4h_reg();
_IKt3_reg();
_Ikt2m2h_reg();
_ipulse2_reg();
_nmitral_reg();
}
