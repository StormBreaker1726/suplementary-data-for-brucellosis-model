#include "../include/mesh.h"

void generate_mesh(unsigned sx, unsigned sy, unsigned dif)
{
    char mesh_file_path[80];

    FILE *msh_file;

    tinymt64_t prng;
    long       seed = (unsigned)(time(NULL) * time(NULL) * time(NULL) * 2.5 * 2.67564 * 7.8734546);
    tinymt64_init(&prng, seed);

    sprintf(mesh_file_path, "meshes/%d_%d_%ld.msh", sx, sy, dif);
    msh_file = fopen(mesh_file_path, "w");
    if (msh_file == NULL)
    {
        printf("%s\n", mesh_file_path);
        printf("nao conseguiu abrir...\n");
        return;
    }

    fprintf(msh_file, "mesh for brucella abortus simulation\nsx = %d\nsy = %d\n", sx, sy);

    mesh *new_mesh  = (mesh *)malloc(sizeof(mesh));
    new_mesh->cells = (points *)malloc((sx * sy) * sizeof(points));

    for (size_t i = 0; i < sx; i++)
    {
        for (size_t j = 0; j < sy; j++)
        {
            real_cpu num                   = tinymt64_generate_double(&prng);
            new_mesh->cells[i * sy + j].ix = i;
            new_mesh->cells[i * sy + j].jy = j;
            if (num <= 0.15)
            {
                new_mesh->cells[i * sy + j].type = 0;
            }
            else if (num > 0.15 && num <= 0.17)
            {
                new_mesh->cells[i * sy + j].type = 1;
            }
            else
            {
                new_mesh->cells[i * sy + j].type = 2;
            }
        }
    }

    for (size_t i = 0; i < sx; i++)
    {
        for (size_t j = 0; j < sy; j++)
        {
            fprintf(msh_file,
                    "ix = %d, jy = %d, type = %d\n",
                    new_mesh->cells[i * sy + j].ix,
                    new_mesh->cells[i * sy + j].jy,
                    new_mesh->cells[i * sy + j].type);
        }
    }

    fclose(msh_file);

    free(new_mesh->cells);
    free(new_mesh);
}

void read_mesh_from_file(mesh *to_save, const char *mesh_file_path)
{
    FILE *file = fopen(mesh_file_path, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        exit(EXIT_FAILURE);
    }
    // printf("Line = %d\n", __LINE__);
    // Leitura das dimensões sx e sy
    fscanf(file, "mesh for brucella abortus simulation\n");
    // printf("Line = %d\n", __LINE__);
    fscanf(file, "sx = %u\n", &(to_save->sx));
    // printf("Line = %d\n", __LINE__);
    fscanf(file, "sy = %u\n", &(to_save->sy));
    // printf("Line = %d\n", __LINE__);
    // printf("Line = %d\n", __LINE__);
    // Alocação dinâmica de memória para o vetor cells
    to_save->cells = (points *)malloc(to_save->sx * to_save->sy * sizeof(points));
    if (to_save->cells == NULL)
    {
        perror("Erro ao alocar memória para cells");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Leitura das coordenadas e do tipo para cada célula
    for (unsigned i = 0; i < to_save->sx * to_save->sy; ++i)
    {
        fscanf(file, "ix = %u, jy = %u, type = %u\n", &(to_save->cells[i].ix), &(to_save->cells[i].jy), &(to_save->cells[i].type));
    }

    fclose(file);
}

void save_result_bacteria(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/bacteria/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}

void save_result_cytokine(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/cytokine/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}

void save_result_apc_naive(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/apc_n/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}

void save_result_apc_active_tct(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/apc_a/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}

void save_result_antibody_tct(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/antibody/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}

void save_result_tissue(real_cpu *concentration_matrix, unsigned step, unsigned sx, unsigned sy, real_cpu h)
{
    FILE *vtk_file;
    char  filename[80];
    sprintf(filename, "results/tissue/concentration%d.vtk", step);

    vtk_file = fopen(filename, "w");

    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "results.vtk\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET RECTILINEAR_GRID\n");
    fprintf(vtk_file, "DIMENSIONS %d %d 1\n", sx, sy);

    fprintf(vtk_file, "X_COORDINATES %d double\n", sx);
    for (int i = 0; i < sx; i++)
    {
        fprintf(vtk_file, "%f ", i * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Y_COORDINATES %d double\n", sy);
    for (int j = 0; j < sy; j++)
    {
        fprintf(vtk_file, "%f ", j * (h / 1));
    }
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "Z_COORDINATES 1 double\n");
    fprintf(vtk_file, "0");
    fprintf(vtk_file, "\n");

    fprintf(vtk_file, "POINT_DATA %d \n", sy * sx * 1);
    fprintf(vtk_file, "FIELD FieldData 1 \n");
    fprintf(vtk_file, "Concentração 1 %d double \n", sy * sx * 1);
    for (int j = 0; j < sx; j++)
    {
        for (int i = 0; i < sy; i++)
        {
            fprintf(vtk_file, "%f \n", concentration_matrix[j * sy + i]);
        }
    }
    fprintf(vtk_file, "\n");
    fclose(vtk_file);
}
