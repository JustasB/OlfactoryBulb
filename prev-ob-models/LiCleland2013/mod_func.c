#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaPN_reg();
extern void _CaT_reg();
extern void _Caint_reg();
extern void _Can_reg();
extern void _GradNMDA_reg();
extern void _GradeAMPA_reg();
extern void _GradeGABA_reg();
extern void _KCa_reg();
extern void _KDRmt_reg();
extern void _KS_reg();
extern void _LCa_reg();
extern void _NaP_reg();
extern void _Naxn_reg();
extern void _Nicotin_reg();
extern void _OdorInput_reg();
extern void _cadecay_reg();
extern void _cadecay2_reg();
extern void _hpg_reg();
extern void _kAmt_reg();
extern void _kM_reg();
extern void _kfasttab_reg();
extern void _kslowtab_reg();
extern void _nafast_reg();
extern void _nmdanet_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaPN.mod");
fprintf(stderr," CaT.mod");
fprintf(stderr," Caint.mod");
fprintf(stderr," Can.mod");
fprintf(stderr," GradNMDA.mod");
fprintf(stderr," GradeAMPA.mod");
fprintf(stderr," GradeGABA.mod");
fprintf(stderr," KCa.mod");
fprintf(stderr," KDRmt.mod");
fprintf(stderr," KS.mod");
fprintf(stderr," LCa.mod");
fprintf(stderr," NaP.mod");
fprintf(stderr," Naxn.mod");
fprintf(stderr," Nicotin.mod");
fprintf(stderr," OdorInput.mod");
fprintf(stderr," cadecay.mod");
fprintf(stderr," cadecay2.mod");
fprintf(stderr," hpg.mod");
fprintf(stderr," kAmt.mod");
fprintf(stderr," kM.mod");
fprintf(stderr," kfasttab.mod");
fprintf(stderr," kslowtab.mod");
fprintf(stderr," nafast.mod");
fprintf(stderr," nmdanet.mod");
fprintf(stderr, "\n");
    }
_CaPN_reg();
_CaT_reg();
_Caint_reg();
_Can_reg();
_GradNMDA_reg();
_GradeAMPA_reg();
_GradeGABA_reg();
_KCa_reg();
_KDRmt_reg();
_KS_reg();
_LCa_reg();
_NaP_reg();
_Naxn_reg();
_Nicotin_reg();
_OdorInput_reg();
_cadecay_reg();
_cadecay2_reg();
_hpg_reg();
_kAmt_reg();
_kM_reg();
_kfasttab_reg();
_kslowtab_reg();
_nafast_reg();
_nmdanet_reg();
}
