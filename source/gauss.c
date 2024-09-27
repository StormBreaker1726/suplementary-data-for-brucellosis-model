#include "../include/gauss.h"

void tissue_integration(ts *tissue_here, real_cpu *Antibody_tissue, real_cpu *APC_a_tissue, real_cpu *integral)
{
    integral[0] = 0.0f;
    integral[1] = 0.0f;

    real_cpu integral_blood, integral_linf;

    integral_blood = 0.0;
    integral_linf  = 0.0f;

    for (size_t i = 0; i < tissue_here->tissue_mesh->sx; i++)
    {
        for (size_t j = 0; j < tissue_here->tissue_mesh->sy; j++)
        {
            integral_blood += (tissue_here->tissue_mesh->cells[i * tissue_here->tissue_mesh->sy + j].type == 1) ?
                                  Antibody_tissue[i * tissue_here->tissue_mesh->sy + j] :
                                  0.0f;

            integral_linf +=
                (tissue_here->tissue_mesh->cells[i * tissue_here->tissue_mesh->sy + j].type == 0) ? APC_a_tissue[i * tissue_here->tissue_mesh->sy + j] : 0.0f;
        }
    }

    integral_blood = integral_blood * (1.0f / (tissue_here->tissue_mesh->qtd_blood));
    integral_linf  = integral_linf * (1.0f / (tissue_here->tissue_mesh->qtd_linf));

    integral[0] = integral_blood;
    integral[1] = integral_linf;
}