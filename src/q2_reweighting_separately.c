#include <TMath.h>
#include "selector.h"

Double_t selector::q2_reweighting_separately(Float_t Mc_q2, TString mc_type, TString period)
{
  Double_t res = 1.;
  if(mc_type == "mc_bg_norad") {
    if(period == "0405e") {
      res =  1.09097e+00 - 2.01635e-03*Mc_q2;
    }
    if(period == "06e") {
      res = 1.05911 - 0.00217231*Mc_q2;
    }
    if(period == "0607p") {
      res = 1.23541 - 0.00404778*Mc_q2;
    }
  }

  if(mc_type == "mc_prph") {
    if(period == "0405e") {
      res = 1.23541-0.00404778*Mc_q2;
    }
    if(period == "06e") {
      res = 1.1776-0.00315364*Mc_q2;
    }
  }
  return res;
}
