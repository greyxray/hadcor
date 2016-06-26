#include <TMath.h>
#include "selector.h"

Double_t selector::q2_reweighting(Float_t Mc_q2, TString mc_type)
{
  Double_t res = 1.;
  if(mc_type == "mc_bg_norad") {
    res = 1.02464 - 0.00062131*Mc_q2;
  }
  if(mc_type == "mc_prph") {
    res = 1.29184 - 0.00553934*Mc_q2 + 9.1547e-06*Mc_q2*Mc_q2;
  }
  return res;
}
