#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <iterator>
using namespace std;

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TF1.h>
#include <TStyle.h>
#include <TString.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TText.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>

#include "plot_style_utils.h"
#include "../inc/constants.h"
#include "AFG.h"
//#include "plot_style_utils.h"

void dout() 
{
    std::cout << std::endl; 
}
template <typename Head, typename... Tail>
void dout(Head H, Tail... T) 
{
  std::cout << H << ' ';
  dout(T...);
}

int main(int argc, char *argv[])
{
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleH(0.08);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleY(0.99);
	gStyle->SetTitleX(0.15);
	TString q2_sufix = "";

	////////////////////////////////////////////////////////////////////////////////////
	if (argc > 1) 
		{
			q2_sufix = (TString)argv[1];
			q2_sufix = "_" + q2_sufix;
		}
	cout << "q2_sufix:" << q2_sufix << endl;

	TString filename = TString("../mc_prph0405e_parton") + q2_sufix + TString(".root");
	cout << "filename:" << filename << endl;
	TFile *file2 = new TFile(filename, "read"); 
	//TFile *file2 = new TFile("../temp_root2/temp2.root", "read"); 
	
	//TFile *file2 = new TFile("../mc_prph0405e_parton_really_Ian_without_EM_conserv.root", "read"); // before the constrain on EM cons hist_part: 62370 62370 hist_had: 52412 52412
	//TFile *file2 = new TFile("../mc_prph0405e_parton_really_Ian_with_EM_conserv.root", "read"); //Ian? hist_part: 46993 46993 hist_had: 39647 39647
	
	//TFile *file2 = new TFile("../mc_prph0405e_parton_really_ne_ne_without_EM_conserv.root", "read");//TFile *file2 = new TFile("../new_new_MC/mc_prph0405e_parton_new_mc.root", "read");  //new-new .root
	//TFile *file2 = new TFile("../mc_prph0405e_parton_really_ne_ne_with_EM_conserv.root", "read"); //TFile *file2 = new TFile("../mc_prph0405e_parton_new_new_with_EM_conserv.root", "read");  //new-new .root with EM conserv
	////////////////////////////////////////////////////////////////////////////////////
	
	cout << "files attached" << endl;

	const Int_t n_hist = 12;

	TString s_var[n_hist] = {"et", "eta", "Q2", "x", "et_jet", "eta_jet", "xgamma", "xp", "dphi", "deta", "dphi_e_ph", "deta_e_ph"};
  TString s_var_cxx[n_hist] = {"et", "eta", "Q2", "x", "et_jet", "eta_jet", "xgamma", "xp", "dphi", "deta", "dphi_e_gamma", "deta_e_gamma"};
	TString s_dim[n_hist] = {"E_{T}^{#gamma} (GeV)", "#eta^{#gamma}", "Q^{2} (GeV^{2})", "x", "E_{T}^{jet} (GeV)", "#eta^{jet}", "x_{#gamma}", "x_{p}", "#Delta#phi", "#Delta#eta", "#Delta#phi_{e,#gamma}", "#Delta#eta_{e,#gamma}"};
	TString s_nbins[n_hist] = {"number_etbins", "number_etabins", "number_Q2bins", "number_xbins", "number_et_jetbins", "number_eta_jetbins", "number_xgamma_bins", "number_xp_bins", "number_dphi_bins", "number_deta_bins", "number_dphi_e_ph_bins", "number_deta_e_ph_bins"};
  TString s_hist_had[n_hist];
	TString s_hist_part[n_hist];
	TString s_hist_had_nopart[n_hist];
	TString s_hist_part_nohad[n_hist];

	for(Int_t i = 0; i < n_hist; i++) 
	{
		s_hist_had[i] = "Cross_Sections_Histograms/h_had_cross_" + s_var[i];
		s_hist_part[i] = "Cross_Sections_Histograms/h_part_cross_" + s_var[i];
		s_hist_had_nopart[i] = "Cross_Sections_Histograms/h_had_nopart_cross_" + s_var[i];
		s_hist_part_nohad[i] = "Cross_Sections_Histograms/h_part_nohad_cross_" + s_var[i];
	}
	//  TString s_hist2[n_hist] = {"h_ana_DetkLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_JetSel_EtaBreitJets"};

	TString s_legend[n_hist] = 
	{
		"electron method",
		"Double-Angle (zufos)",
		"Double-Angle (cells)",
		"Jacquet-Blondel (zufos)"//,
		//			      "Jacquet-Blondel (cells)",
		//			      "Q^{2} at generated level"
	};
	TH1D *hist_had[n_hist];
	TH1D *hist_had_nopart[n_hist];
	TH1D *hist_part[n_hist];
	TH1D *hist_part_nohad[n_hist];
	TH1D *hist_had_to_part[n_hist];

	TH1D *hist_had_prph[n_hist];
	TH1D *hist_had_nopart_prph[n_hist];
	TH1D *hist_part_prph[n_hist];
	TH1D *hist_part_nohad_prph[n_hist];
	TH1D *hist_had_to_part_prph[n_hist];
	for(Int_t i = 0; i < n_hist; i++)
	{
		file2->cd();
		cout << "reading "  << i << " ..." << endl;
		hist_had_prph[i] = (TH1D*)file2->Get(s_hist_had[i])->Clone();
		hist_had_prph[i]->SetName( (TString)hist_had_prph[i]->GetName() + "_prph");
		hist_had_nopart_prph[i] = (TH1D*) file2->Get(s_hist_had_nopart[i])->Clone();
		hist_had_nopart_prph[i]->SetName( (TString)hist_had_nopart_prph[i]->GetName() + "_prph");
		cout << "\tOK\t" ;

		hist_part_prph[i] = (TH1D*)file2->Get(s_hist_part[i])->Clone();
		hist_part_prph[i]->SetName( (TString)hist_part_prph[i]->GetName() + "_prph");
		hist_part_nohad_prph[i] = (TH1D*) file2->Get(s_hist_part_nohad[i])->Clone();
		hist_part_nohad_prph[i]->SetName( (TString)hist_part_nohad_prph[i]->GetName() + "_prph");
		cout << "\tOK\t";

		hist_had_to_part_prph[i] = (TH1D*)hist_had_prph[i]->Clone();
		hist_had_to_part_prph[i]->SetName((TString)hist_had_prph[i]->GetTitle() + "_ratio");
		hist_had_to_part_prph[i]->Divide(hist_part_prph[i]);
		hist_had_to_part_prph[i]->SetLineWidth(1);
		hist_had_to_part_prph[i]->SetLineColor(kBlue);
		cout << "\tOK" << endl;

		hist_had_prph[i]->Sumw2();
		hist_had_nopart_prph[i]->Sumw2();
		hist_part_prph[i]->Sumw2();
		hist_part_nohad_prph[i]->Sumw2();
		hist_had_to_part_prph[i]->Sumw2();
	}
	cout << "histos succesfully read" << endl;

 	cout << "Integrals" << endl;
 	cout << hist_part_prph[0]->Integral() << endl;
	cout << hist_had_prph[0]->Integral() << endl;
 	cout << hist_part_prph[0]->GetEntries() << endl;
	cout << hist_had_prph[0]->GetEntries() << endl;

	//
	// calculate uncertainties of the ratio according to the Gavin's McCance's thesis (appendix A) - ratio of the correlated histograms
	//
	// C = hadron && parton 
	// D = hadron && !parton (hist_had_nopart)
	// E = parton && !hadron (hist_part_nohad)
	// C + D = hist_had, C + E = hist_part

	// for QQ
	for(Int_t j = 0; j < n_hist; j++)
	{
		Double_t CE_sum = 0, CD_sum = 0, E_sum = 0, D_sum = 0, C_sum = 0, C_copy_sum = 0;
		for(Int_t i = 0; i<hist_had_to_part_prph[j]->GetNbinsX(); i++) {
			Double_t CE = hist_part_prph[j]->GetBinContent(i+1)* hist_part_prph[j]->GetBinWidth(i+1);
			Double_t CD = hist_had_prph[j]->GetBinContent(i+1)* hist_had_prph[j]->GetBinWidth(i+1);
			Double_t E = hist_part_nohad_prph[j]->GetBinContent(i+1)* hist_part_nohad_prph[j]->GetBinWidth(i+1);
			Double_t D = hist_had_nopart_prph[j]->GetBinContent(i+1)* hist_had_nopart_prph[j]->GetBinWidth(i+1);
			Double_t C = CE - E;
			Double_t C_copy = CD - D;
			CE_sum += CE;
			CD_sum += CD;
			E_sum += E;
			D_sum += D;
			C_sum += C;
			C_copy_sum += C_copy;
			cout << "hist " << j << " bin " << i << ": " << C << " " << C_copy << endl;
			if(C != C_copy) {
				cout << C << " != " << C_copy << endl;
				//	exit(-1);
			}
			cout << s_var[j] <<  "\t\tCE = " << CE << ", CD = " << CD << ", E = " << E << ", D = " << D << endl;

			Double_t err = 0;
			Double_t err1 = 1, err2 = 1;
			err1 = C*C*(D+E) + D*D*CE + E*E*CD + 2*C*D*E;
			err2 = TMath::Power(CE, 4);
			if(err2 != 0)
				err = TMath::Sqrt(err1 / err2);
			else {
				//	cout << "C + E = 0! exiting..." << endl; exit(-1);
			}
			hist_had_to_part_prph[j]->SetBinError(i+1, err);
		}
		cout << "integral of hist " << j << endl;
		cout << "hist_part: ";
		cout <<  hist_part_prph[j]->Integral(1, hist_part_prph[j]->GetNbinsX()) << " ";
		cout <<  hist_part_prph[j]->Integral(0, hist_part_prph[j]->GetNbinsX()+1) << endl; //With overunderflow
		cout << "hist_had: ";
		cout <<  hist_had_prph[j]->Integral(1, hist_had_prph[j]->GetNbinsX()) << " ";
		cout <<  hist_had_prph[j]->Integral(0, hist_had_prph[j]->GetNbinsX()+1) << endl;

		cout << "hist_part_nohad: ";
		cout <<  hist_part_nohad_prph[j]->Integral(1, hist_part_nohad_prph[j]->GetNbinsX()) << " ";
		cout <<  hist_part_nohad_prph[j]->Integral(0, hist_part_nohad_prph[j]->GetNbinsX()+1) << endl;
		cout << "hist_had_nopart: ";
		cout <<  hist_had_nopart_prph[j]->Integral(1, hist_had_nopart_prph[j]->GetNbinsX()) << " ";
		cout <<  hist_had_nopart_prph[j]->Integral(0, hist_had_nopart_prph[j]->GetNbinsX()+1) << endl;

		cout << " ==== sums ==== " << endl;
		cout << "CE_sum = " << CE_sum << ", CD_sum = " << CD_sum << ", E_sum = " << E_sum << ", D_sum = " << D_sum 
			<< ", C_sum = " << C_sum << ", C_copy_sum = " << C_copy_sum << endl;
	}

	cout << "before calling hist_maximum: " << endl;
	cout << "============ Errors ================" << endl;
	cout << "hist_had_to_part_prph[i]->GetNbinsX() ";
	cout << hist_had_to_part_prph[0]->GetNbinsX() << endl;
	Double_t max_hadcor = - 1;
	for(Int_t j = 0; j < n_hist; j++)
	{
		cout << "in bins of " << s_var[j] << ": " << endl;
		for(Int_t i = 0; i < hist_had_to_part_prph[j]->GetNbinsX(); i++) 
		{
			cout << "\t" << i << ")";
			cout << "\t\t" << hist_had_prph[j]->GetBinContent(i+1) << " +- " << hist_had_prph[j]->GetBinError(i+1) << " / " <<
			 hist_part_prph[j]->GetBinContent(i+1) << " +- " << hist_part_prph[j]->GetBinError(i+1) 
			 << " = " << hist_had_to_part_prph[j]->GetBinContent(i+1) << " +- " << hist_had_to_part_prph[j]->GetBinError(i+1) 
			 << " || "<< hist_had_to_part_prph[j]->GetBinContent(i+1) \
			 * sqrt(pow(hist_had_prph[j]->GetBinError(i+1) / hist_had_prph[j]->GetBinContent(i+1), 2) \
			 	+ pow(hist_part_prph[j]->GetBinError(i+1) / hist_part_prph[j]->GetBinContent(i+1),2))  << endl;
			 /*
			  << " || "<< hist_had_to_part_prph[j]->GetBinContent(i+1) \
			 * sqrt(pow(hist_had_prph[j]->GetBinError(i+1)/hist_had_prph[j]->GetBinContent(i+1), 2) \
			 	+ pow(hist_part_prph[j]->GetBinError(i+1) /hist_part_prph[j]->GetBinContent(i+1),2))
			 */
			 if (hist_had_to_part_prph[j]->GetBinContent(i+1) > max_hadcor) max_hadcor = hist_had_to_part_prph[j]->GetBinContent(i+1);

			 if (hist_part_prph[j]->GetBinContent(i+1) == 0)
			 {
			 	cout <<"!!! WARNING: division by zero. Will force had.corr. to be zero." << endl;
			 	hist_part_prph[j]->SetBinContent(i+1, 0);
			 }
		}
	}

  dout("====For .cxx");
  for(Int_t j = 0; j < n_hist; j++)
  {
    //cout << "in bins of " << s_var_cxx[j] << ": " << endl;
    cout << "\thadcor_" << s_var_cxx[j] << q2_sufix <<  "[" << s_nbins[j] << "] = {";
    for(Int_t i = 0; i < hist_had_to_part_prph[j]->GetNbinsX(); i++) 
    {
      cout << hist_had_to_part_prph[j]->GetBinContent(i+1) ; //<< hist_had_to_part_prph[j]->GetBinError(i+1) 
      if (i+1 < hist_had_to_part_prph[j]->GetNbinsX()) cout << ", ";
    }
    cout << " }; \n\thadcor_" << s_var_cxx[j] << "_err" << q2_sufix <<  "[" << s_nbins[j] << "] = {";
    for(Int_t i = 0; i < hist_had_to_part_prph[j]->GetNbinsX(); i++) 
    {
      cout << hist_had_to_part_prph[j]->GetBinError(i+1);
      if (i+1 < hist_had_to_part_prph[j]->GetNbinsX()) cout << ", ";
    }
    cout << " }; \n";
  }


	cout << "TH2D *h_window[n_hist]" << endl;
	cout <<"max_hadcor: " << max_hadcor << endl;
	TH2D *h_window[n_hist];
	//max_hadcor = 1.25;
	h_window[0] = new TH2D("h_window_et", "", number_etbins, et_bin, 10, 0.1, max_hadcor);
	h_window[1] = new TH2D("h_window_eta", "", number_etabins, eta_bin_crosssec, 10, 0.1, max_hadcor);
	h_window[2] = new TH2D("h_window_Q2", "", number_Q2bins, Q2_bin, 10, 0.1, max_hadcor);
	h_window[3] = new TH2D("h_window_x", "", number_xbins, x_bin, 10, 0.1, max_hadcor);
	h_window[4] = new TH2D("h_window_et_jet", "", number_et_jetbins, et_jet_bin, 10, 0.1, max_hadcor);
	// h_window[4] = new TH2D("h_window_et_jet", "", number_et_jetbins, 2.5, 28., 10, 0.1, 2.);
	h_window[5] = new TH2D("h_window_eta_jet", "", number_eta_jetbins, eta_jet_bin, 10, 0.1, max_hadcor);
	h_window[6] = new TH2D("h_window_xgamma", "", number_xgamma_bins, xgamma_bin, 10, 0.1, 1.5 * max_hadcor);
	h_window[7] = new TH2D("h_window_xp", "", number_xp_bins, xp_bin, 10, 0.1, max_hadcor);
	h_window[8] = new TH2D("h_window_dphi", "", number_dphi_bins, dphi_bin, 10, 0.1, max_hadcor);
	h_window[9] = new TH2D("h_window_deta", "", number_deta_bins, deta_bin, 10, 0.1, max_hadcor);
	h_window[10] = new TH2D("h_window_dphi_e_ph", "", number_dphi_e_ph_bins, dphi_e_ph_bin, 10, 0.1, max_hadcor);
	h_window[11] = new TH2D("h_window_deta_e_ph", "", number_deta_e_ph_bins, deta_e_ph_bin, 10, 0.1, max_hadcor);
	cout << "TH2D *h_window[n_hist] finished" << endl;

	TCanvas *c1, *c2;
	c1 = new TCanvas("c1", "c1", 800, 600);
	c1->Divide(3, 2);
	make_clean_pads(c1, 6, 0, 0);
	c2 = new TCanvas("c2", "c2", 800, 600);
	c2->Divide(3, 2);
	make_clean_pads(c2, 6, 0, 0);

	for(Int_t i = 0; i < n_hist; i++)
		if (i < 6) sign_window(c1->GetPad(i+1), h_window[i], s_dim[i], "hadron / parton", "", "middle");
		else sign_window(c2->GetPad(i+1 - 6), h_window[i], s_dim[i], "hadron / parton", "", "middle");


	//c1->GetPad(3)->SetLogx();
	//c1->GetPad(4)->SetLogx();
	//c1->GetPad(5)->SetLogx();
	//c2->GetPad(1)->SetLogx();
	//c2->GetPad(2)->SetLogx();

	TF1 *f_unity = new TF1("unity", "1", -10000., 10000.);
	f_unity->SetLineColor(34);
	TLegend *leg = new TLegend(0.3, 0.68, 0.9, 0.93);//0.2, 0.2, 0.8, 0.6
	TLegend *leg2 = new TLegend(0.3, 0.68, 0.9, 0.93);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->AddEntry(hist_had_to_part_prph[0], "QQ", "lf");
	leg2->SetBorderSize(0);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(0);
	leg2->AddEntry(hist_had_to_part_prph[6], "QQ", "lf");
	for(Int_t i = 0; i < n_hist; i++)
	{
		TCanvas* c_temp;
		if (i < 6) c_temp = c1;
		else c_temp = c2; 
		cout << "c_temp = " << c_temp <<" c1, c2: " << c1 << " " << c2 << endl;

		c_temp->GetPad(i+1 - 6 * (i >5))->cd();
		cout << "here\n";
		h_window[i]->DrawClone();    
		f_unity->DrawClone("SAME");
		//hist_had_to_part[i]->DrawClone("][ HIST SAME");
		//hist_had_to_part[i]->SetFillColor(kRed);
		//hist_had_to_part[i]->SetFillStyle(3335);

		cout << "C_QQ in bins of " << s_var[i] << ": " << endl;
		for(Int_t j = 0; j < hist_had_to_part_prph[i]->GetNbinsX(); j++)
			cout << hist_had_to_part_prph[i]->GetBinContent(j+1)  << ", " ;
		cout << endl;

		//hist_had_to_part[i]->DrawClone("][ E2 SAME");
		hist_had_to_part_prph[i]->DrawClone("][ HIST SAME");
		hist_had_to_part_prph[i]->SetFillColor(kBlue);
		hist_had_to_part_prph[i]->SetFillStyle(3353);
		hist_had_to_part_prph[i]->DrawClone("][ E2 SAME");

		//    hist2[i]->DrawClone("E P X0 SAME");
		if(i==0)
			leg->Draw();
		else if (i == 6)
			{leg2->Draw();c1->Print("hadronisation_corrections.png");}
		else if (i == n_hist - 1 )
			c_temp->Print("hadronisation_corrections2.png");
	}
	cout << "I want my hadronisation_corrections" << endl;
	
	cout << "got it" << endl;


	for(Int_t i = 0; i < n_hist; i++) 
	{
		hist_had_prph[i]->Scale(1., "width");
		// hist_had[i]->Scale(1., "width");
		hist_part_prph[i]->Scale(1., "width");
		// hist_part[i]->Scale(1., "width");
	}
	TH2D *h_window3[12];
	h_window3[0] = new TH2D("h_window3_et", "", number_etbins, et_bin, 10, 0.1, 1.25 * hist_part_prph[0]->GetMaximum());
	h_window3[1] = new TH2D("h_window3_eta", "", number_etabins, eta_bin_crosssec, 10, 0.1, 1.25 * hist_part_prph[1]->GetMaximum());
	h_window3[2] = new TH2D("h_window3_Q2", "", number_Q2bins, Q2_bin, 10, 0.1, 1.25 * hist_part_prph[2]->GetMaximum());
	h_window3[3] = new TH2D("h_window3_x", "", number_xbins, x_bin, 10, 0.1, 1.5 * hist_part_prph[3]->GetMaximum());
	h_window3[4] = new TH2D("h_window3_et_jet", "", number_et_jetbins, et_jet_bin, 10, 0.1, 1.8 * hist_part_prph[4]->GetMaximum());
	h_window3[5] = new TH2D("h_window3_eta_jet", "", number_eta_jetbins, eta_jet_bin, 10, 0.1, 1.25 * hist_part_prph[5]->GetMaximum());
	h_window3[6] = new TH2D("h_window3_xgamma", "", number_xgamma_bins, xgamma_bin, 10, 0.1, 1.25 * hist_part_prph[6]->GetMaximum());
	h_window3[7] = new TH2D("h_window3_xp", "", number_xp_bins, xp_bin, 10, 0.1, 1.25 * hist_part_prph[7]->GetMaximum());
	h_window3[8] = new TH2D("h_window3_dphi", "", number_dphi_bins, dphi_bin, 10, 0.1, 1.25 * hist_part_prph[8]->GetMaximum());
	h_window3[9] = new TH2D("h_window3_deta", "", number_deta_bins, deta_bin, 10, 0.1, 1.25 * hist_part_prph[9]->GetMaximum());
	h_window3[10] = new TH2D("h_window3_dphi_e_ph", "", number_dphi_e_ph_bins, dphi_e_ph_bin, 10, 0.1, 1.25 * hist_part_prph[10]->GetMaximum());
	h_window3[11] = new TH2D("h_window3_deta_e_ph", "", number_deta_e_ph_bins, deta_e_ph_bin, 10, 0.1, 1.25 * hist_part_prph[11]->GetMaximum());
	
	TCanvas* c3, *c3_2;
	c3 = new TCanvas("c3", "c3", 800, 600);
	c3->Divide(3, 2);
	make_clean_pads(c3, 6, 0, 0);
	// for(Int_t i = 0; i < 6; i++)
	// 	sign_window(c3->GetPad(i+1), h_window3[i], s_dim[i], "Events / bin width", "", "middle");
	c3_2 = new TCanvas("c3_2", "c3_2", 800, 600);
	c3_2->Divide(3, 2);
	make_clean_pads(c3_2, 6, 0, 0);

	for(Int_t i = 0; i < n_hist; i++)
		if (i < 6) sign_window(c3->GetPad(i+1), h_window3[i], s_dim[i], "Events / bin width", "", "middle");
		else sign_window(c3_2->GetPad(i+1 - 6), h_window3[i], s_dim[i], "Events / bin width", "", "middle");

	//c3->GetPad(3)->SetLogx();
	//c3->GetPad(4)->SetLogx();
	//c3->GetPad(5)->SetLogx();
	//c3_2->GetPad(1)->SetLogx();
	//c3_2->GetPad(2)->SetLogx();

	TLegend *leg3 = new TLegend(0.3, 0.68, 0.9, 0.93, "QQ Pythia");
	leg3->SetBorderSize(0);
	leg3->SetFillColor(0);

	leg3->SetFillStyle(0);
	leg3->AddEntry(hist_had_prph[4], "hadron level", "lf");
	leg3->AddEntry(hist_part_prph[4], "parton level", "lf");

	cout << "\t\tdelta_eta:\n";
	Double_t temp_sum(0), temp_err_sum(0), temp_sum_mult(0), temp_err_sum_mult(0); 
	cout << "\t\tdelta_phi:\n";
	for(Int_t i = 1; i <= hist_part_prph[8]->GetNbinsX(); i++)
	{
		temp_sum += hist_part_prph[8]->GetBinContent(i)* hist_part_prph[8]->GetBinWidth(i);
		temp_err_sum += hist_part_prph[8]->GetBinError(i) * hist_part_prph[8]->GetBinWidth(i);
		cout << i << ")" << hist_part_prph[8]->GetBinContent(i) << " +/- " << hist_part_prph[8]->GetBinError(i) 
					<< " width: " << hist_part_prph[8]->GetBinWidth(i) <<" : "
					<<  "\n";
	}
	cout << "sum:"<< temp_sum<<" +- "<<temp_err_sum << "\n";


	 temp_sum=0;
	  temp_err_sum=0;
	  temp_sum_mult=0;
	  temp_err_sum_mult=0;
	for(Int_t i = 1; i <= hist_part_prph[9]->GetNbinsX(); i++)
	{
		temp_sum += hist_part_prph[9]->GetBinContent(i)* hist_part_prph[9]->GetBinWidth(i);
		temp_err_sum += hist_part_prph[9]->GetBinError(i) * hist_part_prph[9]->GetBinWidth(i);
		cout << i << ")" << hist_part_prph[9]->GetBinContent(i) << " +/- " << hist_part_prph[9]->GetBinError(i) 
					<< " width: " << hist_part_prph[9]->GetBinWidth(i) <<" : "
					<<  "\n";
	}
	cout << "sum:"<< temp_sum<<" +- "<<temp_err_sum << "\n";


	cout << "\t\tdelta_phi_e_g:\n";
	 temp_sum=0;
	  temp_err_sum=0;
	  temp_sum_mult=0;
	  temp_err_sum_mult=0;
	for(Int_t i = 1; i <= hist_part_prph[10]->GetNbinsX(); i++)
	{
		temp_sum += hist_part_prph[10]->GetBinContent(i) * hist_part_prph[10]->GetBinWidth(i);
		temp_err_sum += hist_part_prph[10]->GetBinError(i) * hist_part_prph[10]->GetBinWidth(i);
		cout << i << ")" << hist_part_prph[10]->GetBinContent(i) << " +/- " << hist_part_prph[10]->GetBinError(i) 
					<< " width: " << hist_part_prph[10]->GetBinWidth(i) <<" : "
					<<  "\n";
	}
	cout << "sum:"<< temp_sum<<" +- "<<temp_err_sum << "\n";


	for(Int_t i = 0; i < n_hist; i++)
	{
		cout << i << endl;
		TCanvas* c_temp;
		if (i < 6) c_temp = c3;
		else c_temp = c3_2; 
		cout << "c_temp = " << c_temp <<" c3, c3_2: " << c3 << " " << c3_2 << endl;

		

		c_temp->GetPad(i+1 - 6 * (i >5))->cd();
		h_window3[i]->DrawClone();
		hist_part_prph[i]->SetLineColor(kRed);
		hist_part_prph[i]->DrawClone("][ HIST SAME");
		hist_part_prph[i]->SetFillColor(kRed);
		hist_part_prph[i]->SetFillStyle(3335);
		hist_part_prph[i]->DrawClone("][ E2 SAME");
		hist_had_prph[i]->SetLineColor(kBlue);
		hist_had_prph[i]->DrawClone("][ HIST SAME");
		hist_had_prph[i]->SetFillColor(kBlue);
		hist_had_prph[i]->SetFillStyle(3353);
		hist_had_prph[i]->DrawClone("][ E2 SAME");


		//    hist2[i]->DrawClone("E P X0 SAME");
		if(i==0) leg3->Draw();
		else if (i == 6)
		{
			leg3->Draw();
			c3->Print("had_part_QQ" + q2_sufix + ".png");
		}
		else if (i == n_hist - 1 ) c_temp->Print("had_part_QQ2" + q2_sufix + ".png");
	}



	dout("Producing reweighting weights", s_var[6]);
  	if (q2_sufix.Contains("lt"))
    {
        std::copy(std::begin(all_theory_cs_font_pt25_Q2lt30), std::end(all_theory_cs_font_pt25_Q2lt30), std::begin(all_theory_cs_font_pt25));
        std::copy(std::begin(all_theory_cs_font_pt25_pos_Q2lt30), std::end(all_theory_cs_font_pt25_pos_Q2lt30), std::begin(all_theory_cs_font_pt25_pos));
        std::copy(std::begin(all_theory_cs_font_pt25_neg_Q2lt30), std::end(all_theory_cs_font_pt25_neg_Q2lt30), std::begin(all_theory_cs_font_pt25_neg));
    }
    else if (q2_sufix.Contains("gt"))
    {
        std::copy(std::begin(all_theory_cs_font_pt25_Q2gt30), std::end(all_theory_cs_font_pt25_Q2gt30), std::begin(all_theory_cs_font_pt25));
        std::copy(std::begin(all_theory_cs_font_pt25_pos_Q2gt30), std::end(all_theory_cs_font_pt25_pos_Q2gt30), std::begin(all_theory_cs_font_pt25_pos));
        std::copy(std::begin(all_theory_cs_font_pt25_neg_Q2gt30), std::end(all_theory_cs_font_pt25_neg_Q2gt30), std::begin(all_theory_cs_font_pt25_neg));
    }

    double total_cs = 0;
    for(Int_t i = 1; i <= hist_part_prph[6]->GetNbinsX(); i++)
  		total_cs += all_theory_cs_font_pt25[6][i-1] * hist_part_prph[6]->GetBinWidth(i);

  	temp_sum = 0;
    dout ("hist_part_prph[6]->GetEntries():", hist_part_prph[6]->GetEntries());
  	for(Int_t i = 1; i <= hist_part_prph[6]->GetNbinsX(); i++)
  	{
      float binWidth = hist_part_prph[6]->GetBinWidth(i);
  		double p_i = hist_part_prph[6]->GetBinContent(i) * binWidth / hist_part_prph[6]->GetEntries();
  		//p_i - is a fraction of events following in bin i from parton level
  		//a_i - is a fraction of events following in bin i from AFG 
  		//a_i = sigma_i * bin_width_i / [sum over sigma_j * bin_width_j for j running over all bins] [%]
  		double a_i = all_theory_cs_font_pt25[6][i-1] * binWidth / total_cs ;
  		dout("bin", i, a_i,"/", p_i,"=", a_i/p_i);
  		temp_sum += hist_part_prph[6]->GetBinContent(i) * binWidth ;
  		//dout("bin", i, hist_part_prph[6]->GetBinContent(i)*binWidth , hist_part_prph[6]->Integral(i, i));
  	}
  	dout(hist_part_prph[6]->GetEntries());
  	cout << "sum2:"<< temp_sum<<" +- " << "\n";


  dout("Make a test plot");
  {
    bool separate_plots = false;
    TH1D * hist_part_rew[n_hist];
    TH2D * h_win[n_hist];
    TH1D * AFG[n_hist];
    for(Int_t i = 0; i < n_hist; i++)
    {
      dout("i",i);
      if (all_theory_cs_font_pt25[i] == 0)
      {
        AFG[i] = 0;
        h_win[i] = 0;
        continue;
      }
      AFG[i] = new TH1D("h_AFG_" + TString(i), "", hist_part_prph[i]->GetNbinsX(), all_bins_arrays[i]);
      for(Int_t bin = 0; bin < hist_part_prph[i]->GetNbinsX(); bin++)
      {
          AFG[i]->SetBinContent(bin + 1, all_theory_cs_font_pt25[i][bin]);
          AFG[i]->SetBinError(bin + 1, all_theory_cs_font_pt25_pos[i][bin]);
          dout("\tsetting bin", bin, AFG[i]->GetBinContent(bin + 1), "VS", all_theory_cs_font_pt25[i][bin]);
      }
      dout("1.5 * AFG[i]->GetMaximum()", 1.5 * AFG[i]->GetMaximum());
      if (i != 8 && i != 10) h_win[i] = new TH2D("h_win_" + TString(i), "", hist_part_prph[i]->GetNbinsX(), all_bins_arrays[i], 10, 0.1, 1.5 * AFG[i]->GetMaximum());
      else h_win[i] = new TH2D("h_win_" + TString(i), "", hist_part_prph[i]->GetNbinsX(), all_bins_arrays[i], 10, 0.0, 1.5 * AFG[i]->GetMaximum());
    }

    TCanvas *c1;
    c1 = new TCanvas("c1", "c1", 800, 600);
    c1->Divide(3, 2);
    make_clean_pads(c1, 6, 0, 0);
    for(Int_t j = 0; j < 6; j++)
    {
      dout("j",j);
      sign_window(c1->GetPad(j+1), h_win[j+6], s_dim[j+6], "cross section, pb", "", "middle");
    }

    dout("Proceed");
    for(Int_t i = 0; i < n_hist; i++)
    {
      hist_part_rew[i] = new TH1D("h_cs_" + TString(i), "", hist_part_prph[i]->GetNbinsX(), all_bins_arrays[i]);//(TH1D*) hist_part_prph[i]->Clone();
      hist_part_rew[i]->SetName( (TString)hist_had_prph[i]->GetName() + "_prph_rew");
      hist_part_rew[i]->SetLineWidth(1);
      hist_part_rew[i]->SetLineColor(kBlue);
      hist_part_rew[i]->Sumw2();
      int nbins  = hist_part_prph[i]->GetNbinsX();
      //if (i != 6) continue;
      if (all_theory_cs_font_pt25[i] == 0) continue;
      //TH1D * AFG = new TH1D("h_AFG_" + TString(i), "", hist_part_prph[i]->GetNbinsX(), all_bins_arrays[i]);
      hist_part_rew[i]->SetLineWidth(1);
      hist_part_rew[i]->SetLineColor(kRed);
      dout("scale");
      float entriess =  hist_part_prph[i]->GetEntries();
      float total_cs = 0;
      for(Int_t j = 1; j <= hist_part_prph[i]->GetNbinsX(); j++)
        total_cs += all_theory_cs_font_pt25[i][j - 1] * hist_part_prph[i]->GetBinWidth(j);
      float w = total_cs / entriess;

      dout("canvas");
      TCanvas *c2;
      c2 = new TCanvas("c2", "c2", 800, 600);
      //make_clean_pads(c1->GetPad(0), 1, 0, 0);
      //make_clean_canv(c1);
      if (separate_plots) c2->cd();
      else c1->GetPad(i+1 - 6)->cd();

      dout("go", i);
      for(Int_t bin = 0; bin < nbins; bin++)
      {
       
        Double_t
         part_rew = hist_part_prph[i]->GetBinContent(bin + 1),
         bin_width = hist_part_prph[i]->GetBinWidth(bin + 1),
         acceptance = 1;
         float Lumi = 3552.40;
        /*  //s_var
            lumi_data[0] = 134.003 * 1.01;//  lumi_data[0] = 134.15997;
            lumi_data[1] = 52.4195 * 1.01;//  lumi_data[1] = 54.79574;
            lumi_data[2] = 136.219 * 1.01;//  lumi_data[2] = 142.93778;
            lumi_mc_bg[0] = 271.36;///703;//271.36;
            lumi_mc_bg[1] = 165.26;
            lumi_mc_bg[2] = 364.6;
            lumi_mc_prph = 3552.40;
            total_luminosity_data = 0.;
        */
        float part_QQ_rew_cross_sec = part_rew * w; // / (acceptance * bin_width * Lumi);
        hist_part_rew[i]->SetBinContent(bin + 1, part_QQ_rew_cross_sec);
        hist_part_rew[i]->SetBinError(bin + 1, hist_part_prph[i]->GetBinError(bin + 1) * w);
        if (i == 8)
        {
          //dout("bin", bin, AFG[i]->GetBinContent(bin + 1), "vs", part_QQ_rew_cross_sec);
          dout("bin", bin, AFG[i]->GetBinContent(bin + 1), "vs", hist_part_rew[i]->GetBinContent(bin + 1));
        }
          
      }
      TLegend *leg;
      if (i != 8 && i != 10) leg = new TLegend(0.6, 0.68, 0.9, 0.93);
      else leg = new TLegend(0.15, 0.68, 0.45, 0.93);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->AddEntry(hist_part_rew[i], "Part., scaled", "l");
      leg->AddEntry(AFG[i], "AFG", "l");

      if (separate_plots)
      {
        //sign_window(c1->GetPad(0), hist_part_rew[i], s_dim[i], "cs", s_var[i], "middle");
        /*
          {
            h->GetXaxis()->SetTitleOffset(1.2);
            h->GetYaxis()->SetLabelSize(fontsizey);
            
            if(p->GetLogy())
                  h->GetYaxis()->SetTitleOffset(2.5);
            else
                  h->GetYaxis()->SetTitleOffset(1.5);
            h->GetXaxis()->SetNdivisions(507, kTRUE);
            h->GetYaxis()->SetNdivisions(507, kTRUE);
          }
        */

        Double_t height = 1 - c1->GetPad(0)->GetTopMargin() - c1->GetPad(0)->GetBottomMargin();
        Double_t width = 1 - c1->GetPad(0)->GetLeftMargin() - c1->GetPad(0)->GetRightMargin();
        Double_t fontsizex = 0.08/width; //0.06
        Double_t fontsizey = 0.08/height;
        cout << "pad width is: " << width << endl;
        cout << "pad height is: " << height << endl;
        //h_win[i] = new TH2D("h_win", "", nbins, all_bins_arrays[i], 10, 0.1, 1.5 * AFG[i]->GetMaximum());
        if (i != 8 && i != 10) h_win[i] = new TH2D("h_win_" + TString(i), "", nbins, all_bins_arrays[i], 10, 0.1, 1.5 * AFG[i]->GetMaximum());
        else h_win[i] = new TH2D("h_win_" + TString(i), "", nbins, all_bins_arrays[i], 10, 0.0, 1.25 * AFG[i]->GetMaximum());
    
        h_win[i]->SetStats(kFALSE);
        h_win[i]->SetTitle("Reweighting to AFG control plot," + q2_sufix);
        h_win[i]->GetXaxis()->SetTitle(s_dim[i]);
        h_win[i]->GetYaxis()->SetTitle("cross section, pb");
        h_win[i]->GetXaxis()->CenterTitle();
        h_win[i]->GetYaxis()->CenterTitle();
        //h_win[i]->GetXaxis()->SetTitleSize(fontsizex);
        h_win[i]->GetXaxis()->SetTitleFont(42);
        h_win[i]->GetYaxis()->SetTitleFont(42);
        //h_win[i]->GetXaxis()->SetLabelSize(fontsizex);
        h_win[i]->GetXaxis()->SetLabelFont(42);
        h_win[i]->GetYaxis()->SetLabelFont(42);
        dout(h_win[i]->GetNbinsX(), "and", h_win[i]->GetNbinsY());        
      }
      
      h_win[i]->DrawClone();
      hist_part_rew[i]->Draw("HIST SAME");
      AFG[i]->Draw("SAME");
      leg->Draw();

      if (!separate_plots && i == 11) 
      {
        // TPaveText t(0.1, 0.9, 0.9, 0.99, "NDC"); // left-up
        // //t.SetTextSize(t_size);
        // t.SetFillStyle(0);//0 - transparent 3001 - dashed
        // t.SetTextAlign(22);//center bottom
        //     //t.SetX1(0.4);
        //     //t.SetX2(0.6);
        //     //t.SetY1(0.90);
        //     //t.SetY2(1.03);
        // t.AddText("Reweighting to AFG control plot," + q2_sufix);
        // t.SetFillColor(0);
        // t.SetBorderSize(0);
        // t.Draw("SAME");
        c1->Print("test_plot_rew_cs_" + s_var[i] + q2_sufix + ".png");
      }
      else if (separate_plots) c2->Print("test_plot_rew_cs_" + s_var[i] + q2_sufix + ".png");
    }


    
  }
	return 0;
}


