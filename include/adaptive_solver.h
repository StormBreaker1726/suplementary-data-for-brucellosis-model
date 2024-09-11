#ifndef INATE_IMMUNE_SYSTEM_SOLVER_C_CODE
#define INATE_IMMUNE_SYSTEM_SOLVER_C_CODE

#pragma once
#include "defines.h"
#include "gauss.h"
#include "mesh.h"
#include "tissue.h"

void apply_initial_conditions_adaptive_ims(real_cpu **B_naive, real_cpu **B_active, real_cpu **B_memory, real_cpu **Antibody_linf, real_cpu **APC_active_linf,
                                           real_cpu kt_ode, unsigned n_points); /* function to apply the initial conditions in the adaptive immune system */

// functions to help in RKC solving process
real_cpu f_APCa_reduced(real_cpu APCa, ts *tissue_used);

real_cpu f_Bnaive_reduced(real_cpu Bn, real_cpu APCa, ts *tissue_used);

real_cpu f_Bactive_reduced(real_cpu Ba, real_cpu Bn, real_cpu APCa, ts *tissue_used);

real_cpu f_Bmemory_reduced(real_cpu Bm, real_cpu Ba, ts *tissue_used);

real_cpu f_Ant_reduced(real_cpu Ant, real_cpu Bm, ts *tissue_used, real_cpu *Antibody_migration);

void solve_reduced_adaptive(real_cpu *B_naive, real_cpu *B_active, real_cpu *B_memory, real_cpu *Antibody_linf, real_cpu *APC_a, real_cpu k_t_ode,
                            unsigned solver_iterations, unsigned pos_start, ts *tissue_used, real_cpu APC_a_migrated,
                            real_cpu *Antibody_migration); /* solves the reduced adaptive immune system equations */

void solve_expanded_adaptive(); /* solves the expanded adaptive immune system equations */

void save_adaptive_immune_system(real_cpu *B_naive, real_cpu *B_active, real_cpu *B_memory, real_cpu *Antibody_linf, real_cpu *APC_a, unsigned n_points);

#endif /* INATE_IMMUNE_SYSTEM_SOLVER_C_CODE */
