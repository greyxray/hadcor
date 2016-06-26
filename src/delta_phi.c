#include "TH1D.h"
#include "selector.h"
#include <iostream>
using namespace std;

Double_t selector::delta_phi(Double_t phi1, Double_t phi2)
{
  Double_t res = phi1-phi2;
  while(res > TMath::Pi())
    res -= 2.*TMath::Pi();
  
  while(res < -TMath::Pi())
    res += 2.*TMath::Pi();

  res = TMath::Abs(res);
  return res;
}
