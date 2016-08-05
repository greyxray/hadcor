#include <TMath.h>
#include "selector.h"

Double_t selector::etaph_reweighting(Double_t eta, TString mc_type)
{
  Double_t res = 1.;

  //  if(mc_type == "mc_bg_norad") {
  //    res = 0.996392 + eta * 0.223788;
  //  }
  
  if(mc_type == "mc_bg_norad") {
    res = 0.967547 + eta*0.453864;
  }
  if(mc_type == "mc_prph") {
    res = 1.04285 - eta*0.486939;
  }

  if(res < 0) res = 0;
  
  return res;
}
