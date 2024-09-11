#ifndef TISSUE_C_CODE
#define TISSUE_C_CODE

#pragma once
#include "defines.h"
#include "mesh.h"

typedef struct tissue
{
    real_cpu Db;
    real_cpu Dapcn;
    real_cpu Dapca;
    real_cpu Dc;
    real_cpu Da;
    real_cpu rhoc;
    real_cpu chi;
    real_cpu Kts;
    real_cpu tau_b;
    real_cpu capcb;
    real_cpu cbTs;
    real_cpu cab;
    real_cpu capcc;
    real_cpu mut;
    real_cpu gammaapca;
    real_cpu gammaapcn;
    real_cpu gammaatct;
    real_cpu ommegac;
    real_cpu apcn_max;

    real_cpu C1;
    real_cpu C2;
    real_cpu C3;
    real_cpu C4;
    real_cpu C5;
    real_cpu C6;
    real_cpu C7;
    real_cpu C8;
    real_cpu C9;
    real_cpu C10;
    real_cpu C11;

    real_cpu Bn0;
    real_cpu Bm_max;

    mesh *tissue_mesh;
} ts;

void init(ts *tissue_used, const char *tissue_mesh);

#endif /* TISSUE_C_CODE */
