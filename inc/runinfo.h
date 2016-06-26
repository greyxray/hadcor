#ifndef RUNINFO_H
#define RUNINFO_H
#include <TH1D.h>
class RunInfo
{
 public:
 	RunInfo(Int_t, Double_t, Double_t, Double_t, Double_t);
 	RunInfo(); 
 	~RunInfo();
	Int_t GetRunnr();
	Double_t GetLpol_lumi();
	Double_t GetTpol_lumi();
	Double_t GetLpol_pol();
	Double_t GetTpol_pol();
	void SetRunnr(Int_t);
	void SetLpol_lumi(Double_t);
	void SetTpol_lumi(Double_t);
	void SetLpol_pol(Double_t);
	void SetTpol_pol(Double_t);
 private:
 	Int_t itsRunnr;
 	Double_t itsLpol_lumi;
 	Double_t itsTpol_lumi;
 	Double_t itsLpol_pol;
 	Double_t itsTpol_pol;
 };

#endif

