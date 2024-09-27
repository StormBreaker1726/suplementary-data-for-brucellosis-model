#include "../include/inate_solver.h"

real_cpu chemotaxis(real_cpu apcn_ij, real_cpu apcn_ipj, real_cpu apcn_ijp, real_cpu apcn_imj, real_cpu apcn_ijm, real_cpu cit_ij, real_cpu cit_ipj,
                    real_cpu cit_ijp, real_cpu cit_imj, real_cpu cit_ijm, size_t i, size_t j, unsigned sx, unsigned sy, real_cpu k_t, real_cpu h_s)
{
    real_cpu ax, chemotaxis_x, chemotaxis_y, term;
    ax = chemotaxis_x = chemotaxis_y = term = 0;

    // tratando a quimiotaxia no eixo x

    ax   = (cit_ipj - cit_imj) / (2 * h_s);
    term = (ax * k_t / h_s);
    if ((i == 0 || i == sx - 1) || (j == 0 || j == sy - 1) || (ax == 0))
    {
        chemotaxis_x += 0.0f;
    }
    else if (ax > 0)
    {
        chemotaxis_x += term * (apcn_ij - apcn_imj) / h_s;
    }
    else if (ax < 0)
    {
        chemotaxis_x += term * (apcn_ipj - apcn_ij) / h_s;
    }

    ax   = (cit_ijp - cit_ijm) / (2 * h_s);
    term = (ax * k_t / h_s);
    if ((i == 0 || i == sx - 1) || (j == 0 || j == sy - 1) || (ax == 0))
    {
        chemotaxis_y += 0.0f;
    }
    else if (ax > 0)
    {
        chemotaxis_y += term * (apcn_ij - apcn_ijm) / h_s;
    }
    else if (ax < 0)
    {
        chemotaxis_y += term * (apcn_ijp - apcn_ij) / h_s;
    }

    return (chemotaxis_x + chemotaxis_y);
}

void apply_initial_conditions_inate_ims(real_cpu **Bac, real_cpu **Cytokine, real_cpu **APC_n, real_cpu **APC_a_tissue, real_cpu **Antibody_tissue,
                                        real_cpu **Tissue, real_cpu h, unsigned sx, unsigned sy, unsigned simulation_length, ts *tissue_used, unsigned seed)
{
    unsigned jump  = floor(0.1 * sx);
    unsigned jump2 = floor(0.1 * sx);

    tinymt64_t prng;
    tinymt64_init(&prng, seed);

    real_cpu chance_having_apcn = 0.25;

    for (size_t i = 0; i < sx; i++)
    {
        for (size_t j = 0; j < sy; j++)
        {
            (*APC_a_tissue)[i * sy + j]    = 0.0f;
            (*Antibody_tissue)[i * sy + j] = 0.0f;
            (*Tissue)[i * sy + j]          = 1.0f;
            (*Cytokine)[i * sy + j]        = 0.0f;

            if (i >= floor(sx / 2) && i <= floor(sx / 2) + jump && j >= floor(sy / 2) && j <= floor(sy / 2) + jump)
            {
                (*Bac)[i * sx + j] = 0.0100f;
                // (*Bac)[i * sx + j] = 100.0f;
            }
            else
                (*Bac)[i * sx + j] = 0.0f;

            real_cpu apcn_conc_point = 1.0f;
            real_cpu num             = tinymt64_generate_double(&prng);

            (*APC_n)[i * sx + j] = (num <= chance_having_apcn) ? apcn_conc_point : 0.0f;
        }
    }
}

void solve_model(unsigned simulation_length, unsigned sx, unsigned sy, unsigned omp_threads, unsigned seed)
{
    real_cpu start_simulation = omp_get_wtime();

    real_cpu h_s;
    real_cpu k_t;
    real_cpu k_t_ode;

    ts *tissue_used = (ts *)malloc(sizeof(ts));
    init(tissue_used, "meshes/1cmBy1cmTissue/100_100_0.msh");

    h_s     = (real_cpu)((real_cpu)(tissue_used->tissue_mesh->length_x) / (real_cpu)(sx));
    k_t     = pow(10, -6);
    k_t_ode = pow(10, -3);

    printf("Simulation for mesh with %d x %d points and %d days of duration\n", sx, sy, simulation_length);

    real_cpu malloc_time_init = omp_get_wtime();
    /* Arrays to store the inate immune system concentrations */
    real_cpu *Bac               = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *Cytokine          = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *APC_naive         = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *APC_active_tissue = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *Antibody_tissue   = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *Tissue            = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));

    real_cpu *old_Bac               = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *old_Cytokine          = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *old_APC_naive         = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *old_APC_active_tissue = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *old_Antibody_tissue   = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu *old_Tissue            = (real_cpu *)malloc((sx * sy) * sizeof(real_cpu));
    real_cpu  malloc_time_end       = omp_get_wtime();

    /* Arrays to store the adaptive immune system concentrations */
    real_cpu *B_naive;
    real_cpu *B_active;
    real_cpu *B_memory;
    real_cpu *Antibody_linf;
    real_cpu *APC_active_linf;

    real_cpu u_a, u_b, u_c, u_d;

    u_a = u_b = u_c = u_d = 0;

    unsigned sim_points = ceil((simulation_length) / k_t);

    unsigned sim_jump = ceil(sim_jump / 60);

    printf("Simulation time points = %d\n", sim_points);

    real_cpu hs_sq = (h_s * h_s);

    real_cpu init_conditions_time_init = omp_get_wtime();
    apply_initial_conditions_inate_ims(&old_Bac,
                                       &old_Cytokine,
                                       &old_APC_naive,
                                       &old_APC_active_tissue,
                                       &old_Antibody_tissue,
                                       &old_Tissue,
                                       h_s,
                                       sx,
                                       sy,
                                       simulation_length,
                                       tissue_used,
                                       seed);
    real_cpu init_conditions_time_end = omp_get_wtime();

    real_cpu malloc_time_adaptive_init = omp_get_wtime();
    real_cpu malloc_time_adaptive_end  = omp_get_wtime();

    real_cpu save_time_total = 0;
    real_cpu save_time_init, save_time_end;

    real_cpu init_conditions_adaptive_time_init = omp_get_wtime();
    real_cpu init_conditions_adaptive_time_end  = omp_get_wtime();

    save_time_init = omp_get_wtime();
    save_result_bacteria(old_Bac, 0, sx, sy, h_s);
    save_result_tissue(old_Tissue, 0, sx, sy, h_s);
    save_result_apc_naive(old_APC_naive, 0, sx, sy, h_s);
    save_result_cytokine(old_Cytokine, 0, sx, sy, h_s);
    save_result_antibody_tct(old_Antibody_tissue, 0, sx, sy, h_s);
    save_result_apc_active_tct(old_APC_active_tissue, 0, sx, sy, h_s);
    save_time_end = omp_get_wtime();
    save_time_total += (save_time_end - save_time_init);

    unsigned saving_rate = (unsigned)ceil(sim_points / 500.0f);

    unsigned iterations_ratio = (unsigned)ceil(k_t_ode / k_t);

    printf("Saving rate = %d\n", saving_rate);

    real_cpu solving_adaptive, adaptive_solver_init, adaptive_solver_end;
    solving_adaptive = 0.0f;

    unsigned  save = 1;
    size_t    i, j, t;
    real_cpu *integral_tissue = (real_cpu *)malloc(2 * sizeof(real_cpu));

    real_cpu APC_active_bar = 0;
    real_cpu Antibody_bar   = 0;
    printf("%d\n", __LINE__);
#pragma omp parallel num_threads(omp_threads) default(none) private(i, j, t) shared(old_Bac,                   \
                                                                                        old_Cytokine,          \
                                                                                        old_APC_naive,         \
                                                                                        old_APC_active_tissue, \
                                                                                        old_Antibody_tissue,   \
                                                                                        old_Tissue,            \
                                                                                        Bac,                   \
                                                                                        Cytokine,              \
                                                                                        APC_naive,             \
                                                                                        APC_active_tissue,     \
                                                                                        Antibody_tissue,       \
                                                                                        Antibody_bar,          \
                                                                                        APC_active_bar,        \
                                                                                        Tissue,                \
                                                                                        h_s,                   \
                                                                                        save,                  \
                                                                                        hs_sq,                 \
                                                                                        sx,                    \
                                                                                        sy,                    \
                                                                                        simulation_length,     \
                                                                                        tissue_used,           \
                                                                                        k_t_ode,               \
                                                                                        k_t,                   \
                                                                                        saving_rate,           \
                                                                                        iterations_ratio,      \
                                                                                        solving_adaptive,      \
                                                                                        adaptive_solver_init,  \
                                                                                        adaptive_solver_end,   \
                                                                                        integral_tissue,       \
                                                                                        B_naive,               \
                                                                                        B_active,              \
                                                                                        B_memory,              \
                                                                                        Antibody_linf,         \
                                                                                        APC_active_linf,       \
                                                                                        sim_points,            \
                                                                                        u_a,                   \
                                                                                        u_b,                   \
                                                                                        u_c,                   \
                                                                                        u_d,                   \
                                                                                        save_time_init,        \
                                                                                        save_time_end,         \
                                                                                        save_time_total)
    {
        for (t = 0; t < sim_points; t++)
        {
            APC_active_bar = 0.0f;
#pragma omp for
            for (i = 0; i < sx; i++)
            {
                for (j = 0; j < sy; j++)
                {
                    // treating the border conditions and storing the necessaire data
                    real_cpu bac_ijp = (j == sy - 1) ? 2 * h_s * u_c + old_Bac[i * sy + (sy - 1)] : old_Bac[i * sy + (j + 1)];
                    real_cpu bac_ijm = (j == 0) ? 2 * h_s * u_a + old_Bac[i * sy + 1] : old_Bac[i * sy + (j - 1)];
                    real_cpu bac_ipj = (i == sx - 1) ? 2 * h_s * u_d + old_Bac[(sx - 1) * sy + j] : old_Bac[(i + 1) * sy + j];
                    real_cpu bac_imj = (i == 0) ? 2 * h_s * u_b + old_Bac[(1) * sy + j] : old_Bac[(i - 1) * sy + j];
                    real_cpu bac_ij  = old_Bac[i * sy + j];

                    real_cpu cytokine_ijp = (j == sy - 1) ? 2 * h_s * u_c + old_Cytokine[i * sy + (sy - 1)] : old_Cytokine[i * sy + (j + 1)];
                    real_cpu cytokine_ijm = (j == 0) ? 2 * h_s * u_a + old_Cytokine[i * sy + 1] : old_Cytokine[i * sy + (j - 1)];
                    real_cpu cytokine_ipj = (i == sx - 1) ? 2 * h_s * u_d + old_Cytokine[(sx - 1) * sy + j] : old_Cytokine[(i + 1) * sy + j];
                    real_cpu cytokine_imj = (i == 0) ? 2 * h_s * u_b + old_Cytokine[(1) * sy + j] : old_Cytokine[(i - 1) * sy + j];
                    real_cpu cytokine_ij  = old_Cytokine[i * sy + j];

                    real_cpu apc_naive_ijp = (j == sy - 1) ? 2 * h_s * u_c + old_APC_naive[i * sy + (sy - 1)] : old_APC_naive[i * sy + (j + 1)];
                    real_cpu apc_naive_ijm = (j == 0) ? 2 * h_s * u_a + old_APC_naive[i * sy + 1] : old_APC_naive[i * sy + (j - 1)];
                    real_cpu apc_naive_ipj = (i == sx - 1) ? 2 * h_s * u_d + old_APC_naive[(sx - 1) * sy + j] : old_APC_naive[(i + 1) * sy + j];
                    real_cpu apc_naive_imj = (i == 0) ? 2 * h_s * u_b + old_APC_naive[(1) * sy + j] : old_APC_naive[(i - 1) * sy + j];
                    real_cpu apc_naive_ij  = old_APC_naive[i * sy + j];

                    real_cpu apc_active_ijp =
                        (j == sy - 1) ? 2 * h_s * u_c + old_APC_active_tissue[i * sy + (sy - 1)] : old_APC_active_tissue[i * sy + (j + 1)];
                    real_cpu apc_active_ijm = (j == 0) ? 2 * h_s * u_a + old_APC_active_tissue[i * sy + 1] : old_APC_active_tissue[i * sy + (j - 1)];
                    real_cpu apc_active_ipj =
                        (i == sx - 1) ? 2 * h_s * u_d + old_APC_active_tissue[(sx - 1) * sy + j] : old_APC_active_tissue[(i + 1) * sy + j];
                    real_cpu apc_active_imj = (i == 0) ? 2 * h_s * u_b + old_APC_active_tissue[(1) * sy + j] : old_APC_active_tissue[(i - 1) * sy + j];
                    real_cpu apc_active_ij  = old_APC_active_tissue[i * sy + j];

                    real_cpu antibody_ijp = (j == sy - 1) ? 2 * h_s * u_c + old_Antibody_tissue[i * sy + (sy - 1)] : old_Antibody_tissue[i * sy + (j + 1)];
                    real_cpu antibody_ijm = (j == 0) ? 2 * h_s * u_a + old_Antibody_tissue[i * sy + 1] : old_Antibody_tissue[i * sy + (j - 1)];
                    real_cpu antibody_ipj = (i == sx - 1) ? 2 * h_s * u_d + old_Antibody_tissue[(sx - 1) * sy + j] : old_Antibody_tissue[(i + 1) * sy + j];
                    real_cpu antibody_imj = (i == 0) ? 2 * h_s * u_b + old_Antibody_tissue[(1) * sy + j] : old_Antibody_tissue[(i - 1) * sy + j];
                    real_cpu antibody_ij  = old_Antibody_tissue[i * sy + j];

                    real_cpu tissue_ijp = (j == sy - 1) ? 2 * h_s * u_c + old_Tissue[i * sy + (sy - 1)] : old_Tissue[i * sy + (j + 1)];
                    real_cpu tissue_ijm = (j == 0) ? 2 * h_s * u_a + old_Tissue[i * sy + 1] : old_Tissue[i * sy + (j - 1)];
                    real_cpu tissue_ipj = (i == sx - 1) ? 2 * h_s * u_d + old_Tissue[(sx - 1) * sy + j] : old_Tissue[(i + 1) * sy + j];
                    real_cpu tissue_imj = (i == 0) ? 2 * h_s * u_b + old_Tissue[(1) * sy + j] : old_Tissue[(i - 1) * sy + j];
                    real_cpu tissue_ij  = old_Tissue[i * sy + j];

                    real_cpu bacteria_now, cytokine_now, APCn_now, APCa_now, Antibody_now, Tissue_now;
                    bacteria_now = cytokine_now = APCn_now = APCa_now = Antibody_now = Tissue_now = 0;

                    bacteria_now = k_t * (tissue_used->Db * ((bac_ipj + bac_ijp - 4.0f * bac_ij + bac_imj + bac_ijm) / (hs_sq)) + tissue_used->tau_b * bac_ij -
                                          tissue_used->capcb * bac_ij * apc_naive_ij - tissue_used->cab * bac_ij * antibody_ij) +
                                   bac_ij;
                    bacteria_now = bacteria_now < pow(10, -5) ? 0.0f : bacteria_now;
                    if (bacteria_now < 0.0f && bacteria_now != -0.0f)
                    {
                        printf("Deu ruim nas bactérias!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    cytokine_now = k_t * (tissue_used->Dc * ((cytokine_ipj + cytokine_ijp - 4.0f * cytokine_ij + cytokine_imj + cytokine_ijm) / (hs_sq)) +
                                          tissue_used->rhoc * bac_ij * (tissue_ij + apc_naive_ij) - tissue_used->ommegac * cytokine_ij) +
                                   cytokine_ij;
                    cytokine_now = cytokine_now < pow(10, -7) ? 0.0 : cytokine_now;
                    if (cytokine_now < 0.0f && cytokine_now != -0.0f)
                    {
                        printf("Deu ruim nas citocinas!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    real_cpu chemotaxis_term = 0;
                    chemotaxis_term          = chemotaxis(apc_naive_ij,
                                                 apc_naive_ipj,
                                                 apc_naive_ijp,
                                                 apc_naive_imj,
                                                 apc_naive_ijm,
                                                 cytokine_ij,
                                                 cytokine_ipj,
                                                 cytokine_ijp,
                                                 cytokine_imj,
                                                 cytokine_ijm,
                                                 i,
                                                 j,
                                                 sx,
                                                 sy,
                                                 k_t,
                                                 h_s);

                    real_cpu font_term = (tissue_used->tissue_mesh->cells[i * sy + j].type == 0) ?
                                             (tissue_used->capcc * cytokine_ij * (tissue_used->apcn_max - apc_naive_ij)) :
                                             0.0f;
                    // real_cpu font_term = (tissue_used->capcc * cytokine_ij * (tissue_used->apcn_max - apc_naive_ij));
                    APCn_now = k_t * (tissue_used->Dapcn * ((apc_naive_ipj + apc_naive_ijp - 4.0f * apc_naive_ij + apc_naive_imj + apc_naive_ijm) / (hs_sq)) -
                                      tissue_used->chi * chemotaxis_term + font_term - tissue_used->capcb * apc_naive_ij * bac_ij -
                                      tissue_used->gammaapcn * apc_naive_ij) +
                               apc_naive_ij;
                    APCn_now = APCn_now < pow(10, -7) ? 0.0f : APCn_now;
                    if (APCn_now < 0.0f && APCn_now != -0.0f)
                    {
                        printf("Deu ruim nas APC naive!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    real_cpu m_term_apc = (tissue_used->tissue_mesh->cells[i * sy + j].type == 1) ? APC_active_bar : 0.0f;
                    APCa_now = k_t * (tissue_used->Dapca * ((apc_naive_ipj + apc_naive_ijp - 4.0f * apc_naive_ij + apc_naive_imj + apc_naive_ijm) / (hs_sq)) +
                                      tissue_used->capcb * apc_naive_ij * bac_ij - tissue_used->gammaapca * apc_active_ij - m_term_apc) +
                               apc_active_ij;
                    APCa_now = APCa_now < pow(10, -7) ? 0.0f : APCa_now;
                    if (APCa_now < 0.0f && APCa_now != -0.0f)
                    {
                        printf("Deu ruim nas APC active!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    real_cpu m_term_anti = (tissue_used->tissue_mesh->cells[i * sy + j].type == 0) ? Antibody_bar : 0.0f;

                    Antibody_now = k_t * (tissue_used->Da * ((antibody_ipj + antibody_ijp - 4.0f * antibody_ij + antibody_imj + antibody_ijm) / (hs_sq)) -
                                          tissue_used->gammaatct * antibody_ij + m_term_anti) +
                                   antibody_ij;
                    Antibody_now = Antibody_now < pow(10, -7) ? 0.0f : Antibody_now;
                    if (Antibody_now < 0.0f && Antibody_now != -0.0f)
                    {
                        printf("Deu ruim nos anticorpos!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    Tissue_now =
                        k_t * (tissue_used->mut * (1.0f - (tissue_ij / tissue_used->Kts)) * tissue_ij - tissue_used->cbTs * tissue_ij * bac_ij) + tissue_ij;
                    Tissue_now = Tissue_now < pow(10, -7) ? 0.0f : Tissue_now;
                    if (Tissue_now < 0.0f && Tissue_now != -0.0f)
                    {
                        printf("Deu ruim no tecido!\n");
                        unsigned k = 1;
                        while (k)
                            ;
                    }

                    // atualizando nos vetores
                    Bac[i * sx + j]               = bacteria_now;
                    Cytokine[i * sx + j]          = cytokine_now;
                    APC_naive[i * sx + j]         = APCn_now;
                    APC_active_tissue[i * sx + j] = APCa_now;
                    Antibody_tissue[i * sx + j]   = Antibody_now;
                    Tissue[i * sx + j]            = Tissue_now;
                }
            }

            if (omp_get_thread_num() == 0)
            {
                // printf("Time = %.2e \n\n", t * k_t);
                real_cpu *tmp_Bac;
                real_cpu *tmp_Cytokine;
                real_cpu *tmp_APC_naive;
                real_cpu *tmp_APC_active_tissue;
                real_cpu *tmp_Antibody_tissue;
                real_cpu *tmp_Tissue;

                tmp_Bac               = Bac;
                tmp_Cytokine          = Cytokine;
                tmp_APC_naive         = APC_naive;
                tmp_APC_active_tissue = APC_active_tissue;
                tmp_Antibody_tissue   = Antibody_tissue;
                tmp_Tissue            = Tissue;

                Bac               = old_Bac;
                Cytokine          = old_Cytokine;
                APC_naive         = old_APC_naive;
                APC_active_tissue = old_APC_active_tissue;
                Antibody_tissue   = old_Antibody_tissue;
                Tissue            = old_Tissue;

                old_Bac               = tmp_Bac;
                old_Cytokine          = tmp_Cytokine;
                old_APC_naive         = tmp_APC_naive;
                old_APC_active_tissue = tmp_APC_active_tissue;
                old_Antibody_tissue   = tmp_Antibody_tissue;
                old_Tissue            = tmp_Tissue;

                adaptive_solver_init = omp_get_wtime();
                adaptive_solver_end  = omp_get_wtime();

                solving_adaptive += (adaptive_solver_end - adaptive_solver_init);
                if (t % saving_rate == 0)
                {
                    printf("\n\t\tSaving time = %f\n", t * k_t);
                    save_time_init = omp_get_wtime();
                    save_result_bacteria(Bac, save, sx, sy, h_s);
                    save_result_tissue(Tissue, save, sx, sy, h_s);
                    save_result_apc_naive(APC_naive, save, sx, sy, h_s);
                    save_result_cytokine(Cytokine, save, sx, sy, h_s);
                    save_result_antibody_tct(Antibody_tissue, save, sx, sy, h_s);
                    save_result_apc_active_tct(APC_active_tissue, save, sx, sy, h_s);
                    save++;
                    save_time_end = omp_get_wtime();
                    save_time_total += (save_time_end - save_time_init);
                }
            }
        }
    }

    save_time_init = omp_get_wtime();
    save_result_bacteria(Bac, save, sx, sy, h_s);
    save_result_tissue(Tissue, save, sx, sy, h_s);
    save_result_apc_naive(APC_naive, save, sx, sy, h_s);
    save_result_cytokine(Cytokine, save, sx, sy, h_s);
    save_result_antibody_tct(Antibody_tissue, save, sx, sy, h_s);
    save_result_apc_active_tct(APC_active_tissue, save, sx, sy, h_s);

    save_time_end = omp_get_wtime();
    save_time_total += (save_time_end - save_time_init);

    real_cpu free_time_init = omp_get_wtime();

    free(old_Bac);
    free(old_Cytokine);
    free(old_APC_naive);
    free(old_APC_active_tissue);
    free(old_Antibody_tissue);
    free(old_Tissue);

    free(B_naive);
    free(B_active);
    free(B_memory);
    free(Antibody_linf);
    free(APC_active_linf);

    free(Bac);
    free(Cytokine);
    free(APC_naive);
    free(APC_active_tissue);
    free(Antibody_tissue);
    free(Tissue);

    real_cpu free_time_end = omp_get_wtime();

    real_cpu end_simulation = omp_get_wtime();
    printf("Total time to simulate the model: %.8e seconds\n", (end_simulation - start_simulation));
    printf("\tTime to solve the adaptive immune system: %.8e seconds\n", (solving_adaptive));
    printf("\tTime to solve the inate immune system: %.8e seconds\n", (end_simulation - start_simulation) - (save_time_total) - (solving_adaptive));
    printf("\tTime to save results: %.8e seconds\n", (save_time_total));
    printf("\tTime to allocate the arrays: %.8e seconds\n", (malloc_time_adaptive_end - malloc_time_adaptive_init) + (malloc_time_end - malloc_time_init));
    printf("\tTime to free the arrays: %.8e seconds\n", (free_time_end - free_time_init));
    printf("\tTime to apply initial condition: %.8e seconds\n",
           (init_conditions_adaptive_time_end - init_conditions_adaptive_time_init) + (init_conditions_time_end - init_conditions_time_init));
    printf("\tTime to simulate without saving: %8e seconds\n\n", (end_simulation - start_simulation) - (save_time_total));
}
