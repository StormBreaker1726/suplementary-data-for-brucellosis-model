#ifndef MESH_C_CODE
#define MESH_C_CODE

#pragma once
#include "defines.h"
#include "tinymt64.h"

typedef struct
{
    unsigned type;  // define se é vaso linfático (1) ou vaso sanguíneo (0) ou nenhum dos dois (2) // TODO ter chance de nao ser nenhum nem outro tb
    unsigned ix;    // coordenada x da célula
    unsigned jy;    // coordenada y da célula
} points;

typedef struct
{
    points  *cells;
    unsigned sx;
    unsigned sy;
} mesh;

void save_result_bacteria(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);
void save_result_cytokine(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);
void save_result_apc_naive(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);
void save_result_apc_active_tct(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);
void save_result_antibody_tct(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);
void save_result_tissue(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h);

void generate_mesh(unsigned sx, unsigned sy, unsigned dif);
void read_mesh_from_file(mesh *to_save, const char *mesh_file_path);

#endif /* MESH_C_CODE */
