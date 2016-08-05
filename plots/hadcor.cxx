#include <iostream>
#include <fstream>
#include <vector>
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
//#include "plot_style_utils.h"


int main(int argc, char *argv[])
{
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleH(0.08);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleY(0.99);
	gStyle->SetTitleX(0.15);

	TFile *file2 = new TFile("../mc_prph0405e_parton.root", "read"); 
	cout << "files attached" << endl;

	const Int_t n_hist = 6;
	TString s_var[n_hist] = {"et", "eta", "Q2", "x", "et_jet", "eta_jet"};
	TString s_dim[n_hist] = {"E_{T}^{#gamma} (GeV)", "#eta^{#gamma}", "Q^{2} (GeV^{2})", "x", "E_{T}^{jet} (GeV)", "#eta^{jet}"};
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
		cout << "3 - OK, " ;
		cout << "reading "  << i << " ..." << endl;
		hist_part_prph[i] = (TH1D*)file2->Get(s_hist_part[i])->Clone();
		hist_part_prph[i]->SetName( (TString)hist_part_prph[i]->GetName() + "_prph");
		hist_part_nohad_prph[i] = (TH1D*) file2->Get(s_hist_part_nohad[i])->Clone();
		hist_part_nohad_prph[i]->SetName( (TString)hist_part_nohad_prph[i]->GetName() + "_prph");
		cout << "4 - OK, " ;
		hist_had_to_part_prph[i] = (TH1D*)hist_had_prph[i]->Clone();
		hist_had_to_part_prph[i]->SetName((TString)hist_had_prph[i]->GetTitle() + "_ratio");
		hist_had_to_part_prph[i]->Divide(hist_part_prph[i]);
		hist_had_to_part_prph[i]->SetLineWidth(1);
		hist_had_to_part_prph[i]->SetLineColor(kBlue);
	}
	cout << "histos succesfully read" << endl;

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
			Double_t CE = hist_part_prph[j]->GetBinContent(i+1);
			Double_t CD = hist_had_prph[j]->GetBinContent(i+1);
			Double_t E = hist_part_nohad_prph[j]->GetBinContent(i+1);
			Double_t D = hist_had_nopart_prph[j]->GetBinContent(i+1);
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
			cout << "CE = " << CE << ", CD = " << CD << ", E = " << E << ", D = " << D << endl;

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
		cout <<  hist_part_prph[j]->Integral(0, hist_part_prph[j]->GetNbinsX()+1) << endl;
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
	for(Int_t i = 0; i < n_hist; i++)
	{
		cout << "in bins of " << s_var[i] << ": " << endl;
		for(Int_t i = 0; i < hist_had_to_part_prph[i]->GetNbinsX(); i++) 
		{
			cout << "\t\tbin numerator " << i << ": " << hist_had_prph[0]->GetBinContent(i+1) << " +- " << hist_had_prph[0]->GetBinError(i+1) << endl;
			cout << "\t\t%" << endl;
			cout << "\t\tbin denomenator " << i << ": " << hist_part_prph[0]->GetBinContent(i+1) << " +- " << hist_part_prph[0]->GetBinError(i+1) << endl;
			cout << "\t\t = bin ratio " << i << ": " << hist_had_to_part_prph[0]->GetBinContent(i+1) << " +- " << hist_had_to_part_prph[0]->GetBinError(i+1) << endl;

		}
	}
	TH2D *h_window[n_hist];
	h_window[0] = new TH2D("h_window_et", "", number_etbins, et_bin, 10, 0.1, 1.25);
	h_window[1] = new TH2D("h_window_eta", "", number_etabins, eta_bin_crosssec, 10, 0.1, 1.25);
	h_window[2] = new TH2D("h_window_Q2", "", number_Q2bins, Q2_bin, 10, 0.1, 1.25);
	h_window[3] = new TH2D("h_window_x", "", number_xbins, x_bin, 10, 0.1, 1.25);
	h_window[4] = new TH2D("h_window_et_jet", "", number_et_jetbins, et_jet_bin, 10, 0.1, 1.25);
	// h_window[4] = new TH2D("h_window_et_jet", "", number_et_jetbins, 2.5, 28., 10, 0.1, 2.);
	h_window[5] = new TH2D("h_window_eta_jet", "", number_eta_jetbins, eta_jet_bin, 10, 0.1, 1.25);
	TCanvas* c1;
	c1 = new TCanvas("c1", "c1", 800, 600);
	c1->Divide(3, 2);
	make_clean_pads(c1, 6, 0, 0);
	for(Int_t i = 0; i < n_hist; i++)
		sign_window(c1->GetPad(i+1), h_window[i], s_dim[i], "hadron / parton", "", "middle");


	c1->GetPad(3)->SetLogx();
	c1->GetPad(4)->SetLogx();
	c1->GetPad(5)->SetLogx();
	TF1 *f_unity = new TF1("unity", "1", -10000., 10000.);
	f_unity->SetLineColor(34);
	TLegend *leg = new TLegend(0.2, 0.2, 0.8, 0.6);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(hist_had_to_part_prph[0], "QQ", "lf");
	for(Int_t i = 0; i < n_hist; i++)
	{
		c1->GetPad(i+1)->cd();
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
	}
	cout << "I want my hadronisation_corrections" << endl;
	c1->Print("hadronisation_corrections.png");
	cout << "got it" << endl;


	for(Int_t i = 0; i < n_hist; i++) 
	{
		hist_had_prph[i]->Scale(1., "width");
		// hist_had[i]->Scale(1., "width");
		hist_part_prph[i]->Scale(1., "width");
		// hist_part[i]->Scale(1., "width");
	}

	TH2D *h_window3[6];
	h_window3[0] = new TH2D("h_window3_et", "", number_etbins, et_bin, 10, 0.1, 1.25 * hist_part_prph[0]->GetMaximum());
	h_window3[1] = new TH2D("h_window3_eta", "", number_etabins, eta_bin_crosssec, 10, 0.1, 1.25 * hist_part_prph[1]->GetMaximum());
	h_window3[2] = new TH2D("h_window3_Q2", "", number_Q2bins, Q2_bin, 10, 0.1, 1.25 * hist_part_prph[2]->GetMaximum());
	h_window3[3] = new TH2D("h_window3_x", "", number_xbins, x_bin, 10, 0.1, 1.25 * hist_part_prph[3]->GetMaximum());
	h_window3[4] = new TH2D("h_window3_et_jet", "", number_et_jetbins, et_jet_bin, 10, 0.1, 1.25 * hist_part_prph[4]->GetMaximum());
	h_window3[5] = new TH2D("h_window3_eta_jet", "", number_eta_jetbins, eta_jet_bin, 10, 0.1, 1.25 * hist_part_prph[5]->GetMaximum());
	TCanvas* c3;
	c3 = new TCanvas("c3", "c3", 800, 600);
	c3->Divide(3, 2);
	make_clean_pads(c3, 6, 0, 0);
	for(Int_t i = 0; i<6; i++)
		sign_window(c3->GetPad(i+1), h_window3[i], s_dim[i], "Events / bin width", "", "middle");

	c3->GetPad(3)->SetLogx();
	c3->GetPad(4)->SetLogx();
	c3->GetPad(5)->SetLogx();

	TLegend *leg3 = new TLegend(0.3, 0.68, 0.9, 0.93, "QQ Pythia");
	leg3->SetBorderSize(0);
	leg3->SetFillColor(0);
	leg3->AddEntry(hist_had_prph[4], "hadron level", "lf");
	leg3->AddEntry(hist_part_prph[4], "parton level", "lf");

	for(Int_t i = 0; i<n_hist; i++){
		c3->GetPad(i+1)->cd();
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
		if(i==0) leg3->Draw();
	}

	c3->Print("had_part_QQ.png");

	// Double_t corr_et[number_etbins] = {0.83640628,
	// 	0.88917156,
	// 	0.94435876,
	// 	1.00351536
	// };
	// Double_t corr_eta[number_etabins] = {0.94166076,
	// 	0.90925161,
	// 	0.85739744,
	// 	0.81262133
	// };
	// Double_t corr_Q2[number_Q2bins] = {0.92690714,
	// 	0.89150959,
	// 	0.84141361,
	// 	0.84980063,
	// 	0.87854045
	// };
	// Double_t corr_x[number_xbins] = {0.96148003,
	// 	0.86479910,
	// 	0.81532121,
	// 	0.84837804
	// };
	// Double_t corr_eta_jet[number_eta_jetbins] = {0.65408845,
	// 	0.80294960,
	// 	0.97016302,
	// 	1.09009874
	// };
	// Double_t corr_et_jet[number_et_jetbins] = {0.82223966,
	// 	0.82569344,
	// 	0.87374437,
	// 	0.99634799,
	// 	1.04819327,
	// 	1.08043254
	// };
	// TH1D* hist_had_to_part_sum[n_hist];
	// hist_had_to_part_sum[0] = new TH1D("corr_sum" + s_var[0], "", number_etbins, et_bin);
	// hist_had_to_part_sum[1] = new TH1D("corr_sum" + s_var[1], "", number_etabins, eta_bin_crosssec);
	// hist_had_to_part_sum[2] = new TH1D("corr_sum" + s_var[2], "", number_Q2bins, Q2_bin);
	// hist_had_to_part_sum[3] = new TH1D("corr_sum" + s_var[3], "", number_xbins, x_bin);
	// hist_had_to_part_sum[4] = new TH1D("corr_sum" + s_var[4], "", number_et_jetbins, et_jet_bin);
	// hist_had_to_part_sum[5] = new TH1D("corr_sum" + s_var[5], "", number_eta_jetbins, eta_jet_bin);

	// for(Int_t i = 0; i<number_etbins; i++) {
	// 	hist_had_to_part_sum[0]->SetBinContent(i+1, corr_et[i]);
	// 	hist_had_to_part_sum[0]->SetBinError(i+1, 0);
	// }
	// for(Int_t i = 0; i<number_etabins; i++) {
	// 	hist_had_to_part_sum[1]->SetBinContent(i+1, corr_eta[i]);
	// 	hist_had_to_part_sum[1]->SetBinError(i+1, 0);
	// }
	// for(Int_t i = 0; i<number_Q2bins; i++) {
	// 	hist_had_to_part_sum[2]->SetBinContent(i+1, corr_Q2[i]);
	// 	hist_had_to_part_sum[2]->SetBinError(i+1, 0);
	// }
	// for(Int_t i = 0; i<number_xbins; i++) {
	// 	hist_had_to_part_sum[3]->SetBinContent(i+1, corr_x[i]);
	// 	hist_had_to_part_sum[3]->SetBinError(i+1, 0);
	// }
	// for(Int_t i = 0; i<number_et_jetbins; i++) {
	// 	hist_had_to_part_sum[4]->SetBinContent(i+1, corr_et_jet[i]);
	// 	hist_had_to_part_sum[4]->SetBinError(i+1, 0);
	// }
	// for(Int_t i = 0; i<number_eta_jetbins; i++) {
	// 	hist_had_to_part_sum[5]->SetBinContent(i+1, corr_eta_jet[i]);
	// 	hist_had_to_part_sum[5]->SetBinError(i+1, 0);
	// }
	// TCanvas* c4;
	// c4 = new TCanvas("c4", "c4", 800, 600);
	// c4->Divide(3, 2);
	// make_clean_pads(c4, 6, 0, 0);
	// for(Int_t i = 0; i<6; i++)
	// 	sign_window(c4->GetPad(i+1), h_window[i], s_dim[i], "hadron / parton", "", "middle");

	// c4->GetPad(3)->SetLogx();
	// c4->GetPad(4)->SetLogx();
	// c4->GetPad(5)->SetLogx();

	// for(Int_t i = 0; i<n_hist; i++){
	// 	c4->GetPad(i+1)->cd();
	// 	h_window[i]->DrawClone();
	// 	f_unity->DrawClone("SAME");
	// 	hist_had_to_part_sum[i]->SetLineColor(kBlue);
	// 	hist_had_to_part_sum[i]->DrawClone("][ HIST SAME");
	// 	//    hist_had_to_part_sum[i]->SetFillColor(kBlue);
	// 	//    hist_had_to_part_sum[i]->SetFillStyle(3335);
	// 	hist_had_to_part_sum[i]->DrawClone("][ E2 SAME");
	// }

	// c4->Print("hadronisation.eps");


	return 0;
}
