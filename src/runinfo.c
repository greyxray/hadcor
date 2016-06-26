#include "runinfo.h"
RunInfo::RunInfo(Int_t runnr, 
	Double_t lpol_lumi,
	Double_t tpol_lumi,
	Double_t lpol_pol,
	Double_t tpol_pol)
{
	itsLpol_lumi = lpol_lumi;
	itsTpol_lumi = tpol_lumi;
	itsLpol_pol = lpol_pol;
	itsTpol_pol = tpol_pol;
	itsRunnr = runnr;
}

RunInfo::RunInfo()
{
	itsLpol_lumi = 0;
	itsTpol_lumi = 0;
	itsLpol_pol = 0;
	itsTpol_pol = 0;
	itsRunnr = 0;
}

RunInfo::~RunInfo()
{
}

Int_t RunInfo::GetRunnr()
{
	return itsRunnr;
}

Double_t RunInfo::GetLpol_lumi()
{
	return itsLpol_lumi;
}

Double_t RunInfo::GetTpol_lumi()
{
	return itsTpol_lumi;
}

Double_t RunInfo::GetLpol_pol()
{
	return itsLpol_pol;
}

Double_t RunInfo::GetTpol_pol()
{
	return itsTpol_pol;
}

void RunInfo::SetRunnr(Int_t runnr)
{
	itsRunnr = runnr;
}

void RunInfo::SetLpol_lumi(Double_t lpol_lumi)
{
	itsLpol_lumi = lpol_lumi;
}

void RunInfo::SetTpol_lumi(Double_t tpol_lumi)
{
	itsTpol_lumi = tpol_lumi;
}

void RunInfo::SetLpol_pol(Double_t lpol_pol)
{
	itsLpol_pol = lpol_pol;
}

void RunInfo::SetTpol_pol(Double_t lpol_pol)
{
	itsTpol_pol = lpol_pol;
}
