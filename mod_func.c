#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _capump_reg();
extern void _kv_reg();
extern void _na_reg();
extern void _na12_reg();
extern void _na16_reg();
extern void _spike_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," capump.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr," na12.mod");
fprintf(stderr," na16.mod");
fprintf(stderr," spike.mod");
fprintf(stderr, "\n");
    }
_capump_reg();
_kv_reg();
_na_reg();
_na12_reg();
_na16_reg();
_spike_reg();
}
