#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _cadecay_reg(void);
extern void _kA_reg(void);
extern void _kca3_reg(void);
extern void _kfasttab_reg(void);
extern void _kslowtab_reg(void);
extern void _lcafixed_reg(void);
extern void _nafast_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," cadecay.mod");
    fprintf(stderr," kA.mod");
    fprintf(stderr," kca3.mod");
    fprintf(stderr," kfasttab.mod");
    fprintf(stderr," kslowtab.mod");
    fprintf(stderr," lcafixed.mod");
    fprintf(stderr," nafast.mod");
    fprintf(stderr, "\n");
  }
  _cadecay_reg();
  _kA_reg();
  _kca3_reg();
  _kfasttab_reg();
  _kslowtab_reg();
  _lcafixed_reg();
  _nafast_reg();
}
