#include "TF1.h"
#include "TGraph.h"
#include "TColor.h"
#include "../Algorithms/landauGaussian.h"

void Theory_dEdx_resrange();
Double_t landauGaussian();

TGraph *g_ThdEdxRR_Proton;
TGraph *g_ThdEdxRR_Kaon;
TGraph *g_ThdEdxRR_Pion;
TGraph *g_ThdEdxRR_Muon;


