#ifndef INATE_IMMUNE_SYSTEM_SOLVER_2D_C_CODE
#define INATE_IMMUNE_SYSTEM_SOLVER_2D_C_CODE

#pragma once
#include "defines.h"
#include "mesh.h"
#include "tissue.h"

real_cpu chemotaxis(real_cpu apcn_ij, real_cpu apcn_ipj, real_cpu apcn_ijp, real_cpu apcn_imj, real_cpu apcn_ijm, real_cpu cit_ij, real_cpu cit_ipj,
                    real_cpu cit_ijp, real_cpu cit_imj, real_cpu cit_ijm, size_t i, size_t j, unsigned sx, unsigned sy, real_cpu k_t,
                    real_cpu h_s);  // function to calculate the chemotaxis term

void apply_initial_conditions_inate_ims(real_cpu **Bac, real_cpu **Cytokine, real_cpu **APC_n, real_cpu **APC_a_tissue, real_cpu **Antibody_tissue,
                                        real_cpu **Tissue, real_cpu h, unsigned sx, unsigned sy, unsigned simulation_lengths, ts *tissue_used,
                                        unsigned seed); /* function to apply the initial conditions in the inate immune system */

void solve_model(unsigned simulation_length, unsigned sx, unsigned sy, unsigned omp_threads, unsigned seed); /* solves the inate immune system equations */

#endif /* INATE_IMMUNE_SYSTEM_SOLVER_2D_C_CODE */