#include "../include/adaptive_solver.h"

void apply_initial_conditions_adaptive_ims(real_cpu **B_naive, real_cpu **B_active, real_cpu **B_memory, real_cpu **Antibody_linf, real_cpu **APC_active_linf,
                                           real_cpu kt_ode, unsigned n_points)
{
    for (size_t i = 0; i < n_points; i++)
    {
        (*B_naive)[i]         = 0.0f;
        (*B_active)[i]        = 0.0f;
        (*B_memory)[i]        = 0.0f;
        (*Antibody_linf)[i]   = 0.0f;
        (*APC_active_linf)[i] = 0.0f;
    }
}

real_cpu f_APCa_reduced(real_cpu APCa, ts *tissue_used)
{
    return tissue_used->C1 * APCa - tissue_used->C2 * APCa;
}

real_cpu f_Bnaive_reduced(real_cpu Bn, real_cpu APCa, ts *tissue_used)
{
    return tissue_used->C3 * (tissue_used->Bn0 - Bn) - tissue_used->C4 * APCa * Bn;
}

real_cpu f_Bactive_reduced(real_cpu Ba, real_cpu Bn, real_cpu APCa, ts *tissue_used)
{
    return tissue_used->C4 * APCa * Bn - tissue_used->C5 * Ba - tissue_used->C6 * Ba;
}

real_cpu f_Bmemory_reduced(real_cpu Bm, real_cpu Ba, ts *tissue_used)
{
    return tissue_used->C5 * Ba + tissue_used->C7 * Bm * (1 - (Bm / tissue_used->Bm_max)) - tissue_used->C8 * Bm;
}

real_cpu f_Ant_reduced(real_cpu Ant, real_cpu Bm, ts *tissue_used, real_cpu *Antibody_migration)
{
    real_cpu migration_in_this_pass = tissue_used->C11 * gaussian_quadrature(0, 1);
    (*Antibody_migration) += migration_in_this_pass;
    return tissue_used->C9 * Bm - tissue_used->C10 * Ant - migration_in_this_pass;
}

void solve_reduced_adaptive(real_cpu *B_naive, real_cpu *B_active, real_cpu *B_memory, real_cpu *Antibody_linf, real_cpu *APC_a, real_cpu k_t_ode,
                            unsigned solver_iterations, unsigned pos_start, ts *tissue_used, real_cpu APC_a_migrated, real_cpu *Antibody_migration)
{
    /* 
    A maneira como estou fazendo é alocar só o numero de posições realmente necessárias para o sistema imune adaptativo (o mesmo do sim_lengh) 
    e usar o passo-adaptativo armazenando só o último resultado do passo no vetor de resultados. Dessa maneira, economizo espaço de armazenamento e ainda
    consigo usar o passo.
    */

    // TODO lembrar que tenho que somar o termo de migração em todos os passos mas migrar os anticorpos do linfonodo apenas no último passo

    real_cpu APC_a_accumulated         = 0;
    real_cpu B_naive_accumulated       = 0;
    real_cpu B_active_accumulated      = 0;
    real_cpu B_memory_accumulated      = 0;
    real_cpu Antibody_linf_accumulated = 0;

    real_cpu APC_a_now         = 0;
    real_cpu B_naive_now       = 0;
    real_cpu B_active_now      = 0;
    real_cpu B_memory_now      = 0;
    real_cpu Antibody_linf_now = 0;

    (*Antibody_migration) = 0;

    APC_a_accumulated += (APC_a[pos_start - 1] + APC_a_migrated);
    B_naive_accumulated += B_naive[pos_start - 1];
    B_active_accumulated += B_active[pos_start - 1];
    B_memory_accumulated += B_memory[pos_start - 1];
    Antibody_linf_accumulated += Antibody_linf[pos_start - 1];

    unsigned pos_end = pos_start + solver_iterations;

    real_cpu f1, f2, f3, f4;

    for (size_t i = 0; i < pos_end; i++)
    {
        // resolvendo para o APCa
        f1 = f_APCa_reduced(APC_a_accumulated, tissue_used);
        f2 = f_APCa_reduced(APC_a_accumulated + (k_t_ode / 2.0f) * f1, tissue_used);
        f3 = f_APCa_reduced(APC_a_accumulated + (k_t_ode / 2.0f) * f2, tissue_used);
        f4 = f_APCa_reduced(APC_a_accumulated + (k_t_ode)*f3, tissue_used);
        APC_a_accumulated += (k_t_ode / 6.0f) * (f1 + 2.0f * f2 + 2.0f * f3 + f4);

        // resolvendo para o B naive
        f1 = f_Bnaive_reduced(B_naive_accumulated, APC_a_accumulated, tissue_used);
        f2 = f_Bnaive_reduced(B_naive_accumulated + (k_t_ode / 2.0f) * f1, APC_a_accumulated, tissue_used);
        f3 = f_Bnaive_reduced(B_naive_accumulated + (k_t_ode / 2.0f) * f2, APC_a_accumulated, tissue_used);
        f4 = f_Bnaive_reduced(B_naive_accumulated + (k_t_ode)*f3, APC_a_accumulated, tissue_used);
        B_naive_accumulated += (k_t_ode / 6.0f) * (f1 + 2.0f * f2 + 2.0f * f3 + f4);

        // resolvendo para o B de memória
        f1 = f_Bmemory_reduced(B_memory_accumulated, B_active_accumulated, tissue_used);
        f2 = f_Bmemory_reduced(B_memory_accumulated + (k_t_ode / 2.0f) * f1, B_active_accumulated, tissue_used);
        f3 = f_Bmemory_reduced(B_memory_accumulated + (k_t_ode / 2.0f) * f2, B_active_accumulated, tissue_used);
        f4 = f_Bmemory_reduced(B_memory_accumulated + (k_t_ode)*f3, B_active_accumulated, tissue_used);
        B_memory_accumulated += (k_t_ode / 6.0f) * (f1 + 2.0f * f2 + 2.0f * f3 + f4);

        // resolvendo para o B ativa
        f1 = f_Bactive_reduced(B_active_accumulated, B_naive_accumulated, APC_a_accumulated, tissue_used);
        f2 = f_Bactive_reduced(B_active_accumulated + (k_t_ode / 2.0f) * f1, B_naive_accumulated, APC_a_accumulated, tissue_used);
        f3 = f_Bactive_reduced(B_active_accumulated + (k_t_ode / 2.0f) * f2, B_naive_accumulated, APC_a_accumulated, tissue_used);
        f4 = f_Bactive_reduced(B_active_accumulated + (k_t_ode)*f3, B_naive_accumulated, APC_a_accumulated, tissue_used);
        B_active_accumulated += (k_t_ode / 6.0f) * (f1 + 2.0f * f2 + 2.0f * f3 + f4);

        // resolvendo para o Anticorpo
        f1 = f_Ant_reduced(Antibody_linf_accumulated, B_memory_accumulated, tissue_used, Antibody_migration);
        f2 = f_Ant_reduced(Antibody_linf_accumulated + (k_t_ode / 2.0f) * f1, B_memory_accumulated, tissue_used, Antibody_migration);
        f3 = f_Ant_reduced(Antibody_linf_accumulated + (k_t_ode / 2.0f) * f2, B_memory_accumulated, tissue_used, Antibody_migration);
        f4 = f_Ant_reduced(Antibody_linf_accumulated + (k_t_ode)*f3, B_memory_accumulated, tissue_used, Antibody_migration);
        Antibody_linf_accumulated += (k_t_ode / 6.0f) * (f1 + 2.0f * f2 + 2.0f * f3 + f4);
    }

    APC_a[pos_start]         = APC_a_accumulated;
    B_naive[pos_start]       = B_naive_accumulated;
    B_active[pos_start]      = B_active_accumulated;
    B_memory[pos_start]      = B_memory_accumulated;
    Antibody_linf[pos_start] = Antibody_linf_accumulated;
}

void save_adaptive_immune_system(real_cpu *B_naive, real_cpu *B_active, real_cpu *B_memory, real_cpu *Antibody_linf, real_cpu *APC_a, unsigned n_points)
{
    FILE *b_memory_file;
    FILE *b_active_file;
    FILE *b_naive_file;
    FILE *apc_a_linf_file;
    FILE *antibody_file;

    apc_a_linf_file = fopen("results/apc_a_linf/result_in_time.dat", "w");
    b_naive_file    = fopen("results/b_naive/result_in_time.dat", "w");
    b_active_file   = fopen("results/b_active/result_in_time.dat", "w");
    b_memory_file   = fopen("results/b_memory/result_in_time.dat", "w");
    antibody_file   = fopen("results/antibody_linf/result_in_time.dat", "w");

    for (size_t i = 0; i < n_points; i++)
    {
        fprintf(apc_a_linf_file, "%f", APC_a[i]);
        fprintf(b_naive_file, "%f", B_naive[i]);
        fprintf(b_active_file, "%f", B_active[i]);
        fprintf(b_memory_file, "%f", B_memory[i]);
        fprintf(antibody_file, "%f", Antibody_linf[i]);
    }
}