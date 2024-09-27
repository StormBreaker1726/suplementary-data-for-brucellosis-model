#ifndef GAUSS_C_CODE
#define GAUSS_C_CODE

#pragma once
#include "defines.h"
#include "tissue.h"

void tissue_integration(ts *tissue_here, real_cpu *Antibody_tissue, real_cpu *APC_a_tissue, real_cpu *integral);

#endif /* GAUSS_C_CODE */
