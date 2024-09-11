#include "../include/tissue.h"

// funcionando mas extinguindo as APCn muito rápido (25x25 e 40s) (50x50 e 40s as bactérias morrem muito rápido)
// void init(ts *tissue_used, const char *tissue_mesh)
// {
//     tissue_used->Db          = 0.001f;
//     tissue_used->Dapcn       = 0.001f;
//     tissue_used->Dapca       = 0.001f;
//     tissue_used->Dc          = 0.0001f;
//     tissue_used->Da          = 0.0001f;
//     tissue_used->rhoc        = 0.0020f;
//     tissue_used->chi         = 0.0010f;
//     tissue_used->Kts         = 0.0800f;
//     tissue_used->tau_b       = 0.8000f;
//     tissue_used->capcb       = 0.005f;
//     tissue_used->cbTs        = 0.0070f;
//     tissue_used->cab         = 0.0080f;
//     tissue_used->capcc       = 0.78f;
//     tissue_used->mut         = 0.0800f;
//     tissue_used->gammaapca   = 0.0060f;
//     tissue_used->gammaapcn   = 0.00060f;
//     tissue_used->gammaatct   = 0.0060f;
//     tissue_used->ommegac     = 0.0010f;
//     tissue_used->apcn_max    = 0.5500f;
//     tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));
//     read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
// }

// Problema: quanto maior a taxa de difusão, mais rápido a bactéria se extingue

// funcionando... mas bactéria cresce bem rápido
// void init(ts *tissue_used, const char *tissue_mesh)
// {
//     tissue_used->Db          = 0.0007f;
//     tissue_used->Dapcn       = 0.001f;
//     tissue_used->Dapca       = 0.001f;
//     tissue_used->Dc          = 0.0001f;
//     tissue_used->Da          = 0.0001f;
//     tissue_used->rhoc        = 0.0020f;
//     tissue_used->chi         = 0.0010f;
//     tissue_used->Kts         = 0.0800f;
//     tissue_used->tau_b       = 0.780f;
//     tissue_used->capcb       = 0.005f;
//     tissue_used->cbTs        = 0.0070f;
//     tissue_used->cab         = 0.0080f;
//     tissue_used->capcc       = 0.78f;
//     tissue_used->mut         = 0.0800f;
//     tissue_used->gammaapca   = 0.0060f;
//     tissue_used->gammaapcn   = 0.00060f;
//     tissue_used->gammaatct   = 0.0060f;
//     tissue_used->ommegac     = 0.0010f;
//     tissue_used->apcn_max    = 0.5500f;
//     tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));
//     read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
// }

// funciona para ver a ação da citocina, bactéria e APCn
void init(ts *tissue_used, const char *tissue_mesh)
{
    tissue_used->Db          = 0.0007f;
    tissue_used->Dapcn       = 0.0002f;
    tissue_used->Dapca       = 0.001f;
    tissue_used->Dc          = 0.0001f;
    tissue_used->Da          = 0.0001f;
    tissue_used->rhoc        = 0.020f;
    tissue_used->chi         = 0.00010f;
    tissue_used->Kts         = 0.0800f;
    tissue_used->tau_b       = 0.780f;
    tissue_used->capcb       = 0.05f;
    tissue_used->cbTs        = 0.0070f;
    tissue_used->cab         = 0.0080f;
    tissue_used->capcc       = 0.99f;
    tissue_used->mut         = 0.0800f;
    tissue_used->gammaapca   = 0.0060f;
    tissue_used->gammaapcn   = 0.00060f;
    tissue_used->gammaatct   = 0.0060f;
    tissue_used->ommegac     = 0.00010f;
    tissue_used->apcn_max    = 1.0000f;
    tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));
    read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
}

// void init(ts *tissue_used, const char *tissue_mesh)
// {
//     tissue_used->Db          = 0.0007f;
//     tissue_used->Dapcn       = 0.0002f;
//     tissue_used->Dapca       = 0.001f;
//     tissue_used->Dc          = 0.0001f;
//     tissue_used->Da          = 0.0001f;
//     tissue_used->rhoc        = 0.90f;
//     tissue_used->chi         = 0.00010f;
//     tissue_used->Kts         = 0.7000f;
//     tissue_used->tau_b       = 0.700f;
//     tissue_used->capcb       = 0.65f;
//     tissue_used->cbTs        = 0.0070f;
//     tissue_used->cab         = 0.0080f;
//     tissue_used->capcc       = 0.99f;
//     tissue_used->mut         = 0.0800f;
//     tissue_used->gammaapca   = 0.0060f;
//     tissue_used->gammaapcn   = 0.00060f;
//     tissue_used->gammaatct   = 0.0060f;
//     tissue_used->ommegac     = 0.00010f;
//     tissue_used->apcn_max    = 1.0000f;
//     tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));
//     read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
// }

// void init(ts *tissue_used, const char *tissue_mesh)
// {
//     tissue_used->Db          = 0.01f;
//     tissue_used->Dapcn       = 0.01f;
//     tissue_used->Dapca       = 0.01f;
//     tissue_used->Dc          = 0.01f;
//     tissue_used->Da          = 0.01f;
//     tissue_used->rhoc        = 0.1f;
//     tissue_used->chi         = 0.1f;
//     tissue_used->Kts         = 1.0f;
//     tissue_used->tau_b       = 0.5f;
//     tissue_used->capcb       = 0.01f;
//     tissue_used->cbTs        = 0.01f;
//     tissue_used->cab         = 0.01f;
//     tissue_used->capcc       = 0.01f;
//     tissue_used->mut         = 0.01f;
//     tissue_used->gammaapca   = 0.01f;
//     tissue_used->gammaapcn   = 0.01f;
//     tissue_used->gammaatct   = 0.01f;
//     tissue_used->ommegac     = 0.01f;
//     tissue_used->apcn_max    = 1.0f;
//     tissue_used->tissue_mesh = (mesh *)malloc(sizeof(mesh));
//     read_mesh_from_file(tissue_used->tissue_mesh, tissue_mesh);
// }
