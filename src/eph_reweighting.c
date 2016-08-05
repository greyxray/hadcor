#include <TMath.h>
#include "selector.h"

Double_t selector::eph_reweighting(Double_t e, TString mc_type)
{
  Double_t res = 1.;
  Double_t e_bins[12] = {4, 5.575, 7.15, 8.725, 10.3, 11.875, 13.45, 15.025, 16.6 ,18.175, 19.75, 21.325};
  Double_t prph_rew[11] = {1.19762, 1.00293, 0.883233, 0.833037, 0.798075, 0.869911, 1.024, 1.26381, 0.624993, 0.303846, 0.867265};
  Double_t norad_rew[11] = {1.07293, 1.07388, 0.932607, 0.74603, 0.759989, 0.716628, 0.871822, 0.426204, 1.88063, 2.50714, 1.};
  if(mc_type == "mc_bg_norad") {
    for(Int_t i=0; i<11; i++) 
      if(e > e_bins[i] && e < e_bins[i+1])
	res = norad_rew[i];
  }

  if(mc_type == "mc_prph") {
    for(Int_t i=0; i<11; i++) 
      if(e > e_bins[i] && e < e_bins[i+1])
	res = prph_rew[i];
  }

/*   if(mc_type == "mc_bg_norad") { */
/*     res = 1.32075 - 0.0487086*e; */
/*   } */
/*   if(mc_type == "mc_prph") { */
/*     res = 1.2998 - 0.0423466*e; */
/*   } */


  return res;
}
