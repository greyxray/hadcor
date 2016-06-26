#include "TH1D.h"
#include "selector.h"
#include "constants.h"
#include <iostream>
using namespace std;


Double_t selector::jet_en_corr(Double_t eta, Double_t et, TString period, TString mc_type)
{
  const Int_t           num_eta_bins = 16;
  Double_t        eta_min = -1.5;
  Double_t        eta_max = 1.8;
  Double_t        eta_bins[num_eta_bins+1] = {-1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5};
  Double_t a0_corr040506e[num_eta_bins] = {
    0.0876955,
    0.0185162,
    0.147459,
    0.217467,
    0.182831,
    0.164446,
    0.149922,
    0.202732,
    0.221769,
    0.292769,
    0.371231,
    0.313208,
    0.233267,
    -0.0273743,
    -0.114927,
    -0.216469
  };
  Double_t a1_corr040506e[num_eta_bins] = {
    0.84324, 
    0.887504, 
    0.856274,
    0.85811,
    0.878717,
    0.893775,
    0.905287,
    0.895028,
    0.891519,
    0.869106,
    0.856903,
    0.869514,
    0.876799,
    0.909305,
    0.914599,
    0.964785
  };
  Double_t a0_corr0607p[num_eta_bins] = {
    0.378383,
    0.1652,
    0.131539,
    0.2056,
    0.209384,
    0.156203,
    0.166381,
    0.227577,
    0.22005,
    0.30041,
    0.410341,
    0.337955,
    0.137334,
    -0.00163451,
    -0.099333,
    -0.08347
  };
  Double_t a1_corr0607p[num_eta_bins] = {
    0.731974,
    0.831666,
    0.86174,
    0.860342,
    0.871313,
    0.89491,
    0.902218,
    0.889853,
    0.892855,
    0.868315,
    0.849522,
    0.86249,
    0.894076,
    0.909222,
    0.914492,
    0.935125
  };

  Double_t a0_corr_prph[num_eta_bins] = {
    0.170948,
    0.0178937,
    -0.00947929,
    -0.0118412,
    0.222901,
    -0.0626655,
    -0.0653755,
    0.386311,
    0.0354895,
    -0.0463198,
    0.181643,
    -0.401992,
    0.0429563,
    -0.192041,
    -0.662895,
    0.620101
  };
  Double_t a1_corr_prph[num_eta_bins] = {
    0.846053,
    0.898699,
    0.88945,
    0.933213,
    0.882172,
    0.986765,
    1.00725,
    0.912114,
    0.952197,
    0.948717,
    0.844007,
    1.0041,
    0.97342,
    1.01482,
    1.0375,
    0.931508
  };

  Double_t a1_corr, a0_corr;
  Double_t et_corr = 0.;
  Double_t corfac = 1.;

  for(Int_t i=0; i<num_eta_bins; i++)
    {
      if(eta>eta_bins[i] && eta<eta_bins[i+1])
	{
	  if(period == "0405e" || period == "06e"){
	    a1_corr = a1_corr040506e[i];
	    a0_corr = a0_corr040506e[i];
	  }
	  if(period == "0607p" || period == "04p") {
	    a1_corr = a1_corr0607p[i];
	    a0_corr = a0_corr0607p[i];
	  }
	  if(!Data && mc_type == "mc_prph") {
	    a1_corr = a1_corr_prph[i];
	    a0_corr = a0_corr_prph[i];
	  }
	  et_corr = (et-a0_corr)/a1_corr;
	  //	  et_corr = (-a1_corr[i]+sqrt(a1_corr[i]*a1_corr[i]-4*a2_corr[i]*(a0_corr[i]-et)))*(0.5/a2_corr[i]);
	    ///	  corfac = (1./a1_corr[i]) - (a0_corr[i]/(a1_corr[i] * et));
	}
    }
  //  cout << "return " << et_corr << endl;
  if(check_cuts)
    cout << "return value: " << eta << ", " << et  << ", " << period << ", " << mc_type << ", " << et_corr << ", " << a0_corr << ", " << a1_corr << endl;
  return et_corr;
}

Double_t selector::jet_en_corr(Double_t eta, Double_t et, Double_t* C_scale, Double_t* a0_corr, Double_t* a1_corr, Double_t* a2_corr)
{
  Double_t res = 0.;
  return res;
}
