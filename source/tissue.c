#include "../include/tissue.h"

void init(ts *tissue_used, const char *tissue_mesh)
{
    tissue_used->Db          = 2.00f * pow(10, -5);
    tissue_used->Dapcn       = 1.21f * pow(10, -4);
    tissue_used->Dapca       = 1.21f * pow(10, -4);
    tissue_used->Dc          = 9.22f * pow(10, -5);
    tissue_used->Da          = 1.52f * pow(10, -4);
    tissue_used->rhoc        = 1.00f * pow(10, 0);
    tissue_used->chi         = 1.44f * pow(10, -4);
    tissue_used->Kts         = 1.00f * pow(10, 0);
    tissue_used->tau_b       = 0.85f * pow(10, 0);
    tissue_used->capcb       = 0.55f * pow(10, 0);
    tissue_used->cbTs        = 0.27f * pow(10, 0);
    tissue_used->cab         = 0.27f * pow(10, 0);
    tissue_used->capcc       = 0.57f * pow(10, 0);
    tissue_used->mut         = 0.50f * pow(10, 0);
    tissue_used->gammaapca   = 3.43f * pow(10, 0);
    tissue_used->gammaapcn   = 3.43f * pow(10, 0);
    tissue_used->gammaatct   = 1.00f * pow(10, 0);
    tissue_used->ommegac     = 7.50f * pow(10, 0);
    tissue_used->apcn_max    = 80.00f * pow(10, 0);
    tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));

    read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
}
