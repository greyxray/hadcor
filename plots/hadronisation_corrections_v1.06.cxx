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
    
  //    TFile *file1 = new TFile("../data0405e_nocorr.root");
  //  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/data0607p_hanno.root");
  //  TFile *file1 = new TFile("/data/zenith234b/yegor/cross/nt_out/60211_commonNT.root", "read"); 
  //  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/data0607p_v04b.root");
  //  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/data0607p.root");

  TFile *file1 = new TFile("../mc_bg_rad0405e_parton.root", "read"); 
  TFile *file2 = new TFile("../mc_prph0405e_parton.root", "read"); 
  cout << "files attached" << endl;
  
  const Int_t n_hist = 6;
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el", "h_Q2_el", "h_Zvtx"};
  //  TString s_hist2[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist2[n_hist] = {"h_ana_DetLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_EvtSel_e_el", "h_ana_DetLvl_EvtSel_q2_el", "h_ana_DetLvl_EvtSel_zvtx"};

  TString s_var[n_hist] = {"et", "eta", "Q2", "x", "et_jet", "eta_jet"};
  TString s_dim[n_hist] = {"E_{T}^{#gamma} (GeV)", "#eta^{#gamma}", "Q^{2} (GeV^{2})", "x", "E_{T}^{jet} (GeV)", "#eta^{jet}"};
  TString s_hist_had[n_hist];
  TString s_hist_part[n_hist];
  TString s_hist_had_nopart[n_hist];
  TString s_hist_part_nohad[n_hist];
  for(Int_t i=0; i<n_hist; i++) {
    s_hist_had[i] = "Cross_Sections_Histograms/h_had_cross_" + s_var[i];
    s_hist_part[i] = "Cross_Sections_Histograms/h_part_cross_" + s_var[i];
    s_hist_had_nopart[i] = "Cross_Sections_Histograms/h_had_nopart_cross_" + s_var[i];
    s_hist_part_nohad[i] = "Cross_Sections_Histograms/h_part_nohad_cross_" + s_var[i];
  }
  //  TString s_hist2[n_hist] = {"h_ana_DetkLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_JetSel_EtaBreitJets"};

  TString s_legend[n_hist] = {"electron method",
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

  for(Int_t i=0; i<n_hist; i++)
    {
      file1->cd();
      cout << "reading "  << i << " ..." << endl;
      hist_had[i] = (TH1D*)file1->Get(s_hist_had[i])->Clone();
      hist_had_nopart[i] = (TH1D*) file1->Get(s_hist_had_nopart[i])->Clone();
      cout << "1 - OK, " ;
      cout << "reading "  << i << " ..." << endl;
      hist_part[i] = (TH1D*)file1->Get(s_hist_part[i])->Clone();
      hist_part_nohad[i] = (TH1D*) file1->Get(s_hist_part_nohad[i])->Clone();
      cout << "2 - OK, " ;
      hist_had_to_part[i] = (TH1D*)hist_had[i]->Clone();
      hist_had_to_part[i]->SetName((TString)hist_had[i]->GetTitle() + "_ratio");
      hist_had_to_part[i]->Divide(hist_part[i]);
      hist_had_to_part[i]->SetLineWidth(1);
      hist_had_to_part[i]->SetLineColor(kRed);

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

  // for LL
  for(Int_t j=0; j<n_hist; j++){
    Double_t CE_sum = 0, CD_sum = 0, E_sum = 0, D_sum = 0, C_sum = 0, C_copy_sum = 0;
    for(Int_t i=0; i<hist_had_to_part[j]->GetNbinsX(); i++) {
      Double_t CE = hist_part[j]->GetBinContent(i+1);
      Double_t CD = hist_had[j]->GetBinContent(i+1);
      Double_t E = hist_part_nohad[j]->GetBinContent(i+1);
      Double_t D = hist_had_nopart[j]->GetBinContent(i+1);
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
      hist_had_to_part[j]->SetBinError(i+1, err);
    }
    cout << "integral of hist " << j << endl;
    cout << "hist_part: ";
    cout <<  hist_part[j]->Integral(1, hist_part[j]->GetNbinsX()) << " ";
    cout <<  hist_part[j]->Integral(0, hist_part[j]->GetNbinsX()+1) << endl;
    cout << "hist_had: ";
    cout <<  hist_had[j]->Integral(1, hist_had[j]->GetNbinsX()) << " ";
    cout <<  hist_had[j]->Integral(0, hist_had[j]->GetNbinsX()+1) << endl;

    cout << "hist_part_nohad: ";
    cout <<  hist_part_nohad[j]->Integral(1, hist_part_nohad[j]->GetNbinsX()) << " ";
    cout <<  hist_part_nohad[j]->Integral(0, hist_part_nohad[j]->GetNbinsX()+1) << endl;
    cout << "hist_had_nopart: ";
    cout <<  hist_had_nopart[j]->Integral(1, hist_had_nopart[j]->GetNbinsX()) << " ";
    cout <<  hist_had_nopart[j]->Integral(0, hist_had_nopart[j]->GetNbinsX()+1) << endl;

    cout << " ==== sums ==== " << endl;
    cout << "CE_sum = " << CE_sum << ", CD_sum = " << CD_sum << ", E_sum = " << E_sum << ", D_sum = " << D_sum 
	 << ", C_sum = " << C_sum << ", C_copy_sum = " << C_copy_sum << endl;
  }

  // for QQ
  for(Int_t j=0; j<n_hist; j++){
    Double_t CE_sum = 0, CD_sum = 0, E_sum = 0, D_sum = 0, C_sum = 0, C_copy_sum = 0;
    for(Int_t i=0; i<hist_had_to_part_prph[j]->GetNbinsX(); i++) {
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
  for(Int_t i=0; i<n_hist; i++)
    sign_window(c1->GetPad(i+1), h_window[i], s_dim[i], "hadron / parton", "", "middle");
  

  c1->GetPad(3)->SetLogx();
  c1->GetPad(4)->SetLogx();
  c1->GetPad(5)->SetLogx();
  TF1 *f_unity = new TF1("unity", "1", -10000., 10000.);
  f_unity->SetLineColor(34);
  TLegend *leg = new TLegend(0.2, 0.2, 0.8, 0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist_had_to_part[0], "LL", "lf");
  leg->AddEntry(hist_had_to_part_prph[0], "QQ", "lf");
  for(Int_t i=0; i<n_hist; i++){
    c1->GetPad(i+1)->cd();
    h_window[i]->DrawClone();    
    f_unity->DrawClone("SAME");
    hist_had_to_part[i]->DrawClone("][ HIST SAME");
    hist_had_to_part[i]->SetFillColor(kRed);
    hist_had_to_part[i]->SetFillStyle(3335);
    cout << "C_LL in bins of " << s_var[i] << ": " << endl;
    for(Int_t j=0; j<hist_had_to_part[i]->GetNbinsX(); j++)
      cout << hist_had_to_part[i]->GetBinContent(j+1)  << ", " ;
    cout << endl;
    cout << "C_QQ in bins of " << s_var[i] << ": " << endl;
    for(Int_t j=0; j<hist_had_to_part_prph[i]->GetNbinsX(); j++)
      cout << hist_had_to_part_prph[i]->GetBinContent(j+1)  << ", " ;
    cout << endl;
    hist_had_to_part[i]->DrawClone("][ E2 SAME");
    hist_had_to_part_prph[i]->DrawClone("][ HIST SAME");
    hist_had_to_part_prph[i]->SetFillColor(kBlue);
    hist_had_to_part_prph[i]->SetFillStyle(3353);
    hist_had_to_part_prph[i]->DrawClone("][ E2 SAME");

    //    hist2[i]->DrawClone("E P X0 SAME");
    if(i==0)
      leg->Draw();
  }

  c1->Print("hadronisation_corrections.eps");

  TH2D *h_window2[2];
  h_window2[0] = new TH2D("h_window20", "", number_et_jetbins, et_jet_bin, 10, 0.1, 200.);
  h_window2[1] = new TH2D("h_window21", "", number_et_jetbins, et_jet_bin, 10, 0.1, 7000.);
  TCanvas* c2;
  c2 = new TCanvas("c2", "c2", 800, 500);
  c2->Divide(2, 1);
  make_clean_pads(c2, 2, 0, 0);
  for(Int_t i=0; i<2; i++)
    sign_window(c2->GetPad(i+1), h_window2[i], "E_{T}^{jet} [GeV]", "Events", "", "middle");

  c2->GetPad(1)->cd();
  h_window2[0]->DrawClone();

  hist_part[4]->SetLineColor(kRed);
  hist_part[4]->DrawClone("][ HIST SAME");
  hist_part[4]->SetFillColor(kRed);
  hist_part[4]->SetFillStyle(3335);
  hist_part[4]->DrawClone("][ E2 SAME");
  hist_had[4]->SetLineColor(kBlue);
  hist_had[4]->DrawClone("][ HIST SAME");
  hist_had[4]->SetFillColor(kBlue);
  hist_had[4]->SetFillStyle(3353);
  hist_had[4]->DrawClone("][ E2 SAME");

  TLegend *leg2 = new TLegend(0.2, 0.68, 0.8, 0.93, "LL Ariadne");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(hist_had[4], "hadron level", "lf");
  leg2->AddEntry(hist_part[4], "parton level", "lf");
  leg2->Draw();
  c2->GetPad(2)->cd();
  h_window2[1]->DrawClone();

  hist_part_prph[4]->SetLineColor(kRed);
  hist_part_prph[4]->DrawClone("][ HIST SAME");
  hist_part_prph[4]->SetFillColor(kRed);
  hist_part_prph[4]->SetFillStyle(3335);
  hist_part_prph[4]->DrawClone("][ E2 SAME");
  hist_had_prph[4]->SetLineColor(kBlue);
  hist_had_prph[4]->DrawClone("][ HIST SAME");
  hist_had_prph[4]->SetFillColor(kBlue);
  hist_had_prph[4]->SetFillStyle(3353);
  hist_had_prph[4]->DrawClone("][ E2 SAME");

  TLegend *leg3 = new TLegend(0.2, 0.68, 0.8, 0.93, "QQ Pythia");
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->AddEntry(hist_had_prph[4], "hadron level", "lf");
  leg3->AddEntry(hist_part_prph[4], "parton level", "lf");
  leg3->Draw();

  c2->Print("jet_et_had_part.eps");


  TCanvas* c3;
  c3 = new TCanvas("c3", "c3", 800, 500);
  c3->Divide(1, 1);
  make_clean_pads(c3, 1, 0, 0);
  for(Int_t i=0; i<1; i++)
    sign_window(c3->GetPad(i+1), h_window[4], "E_{T}^{jet} [GeV]", "hadron / parton", "", "middle");

  c3->GetPad(1)->cd();
  h_window[4]->DrawClone();
  //  hist_had_to_part[4]->DrawClone("][ HIST SAME");   
  hist_had_to_part[4]->SetMarkerStyle(20);
  hist_had_to_part[4]->SetMarkerColor(kRed);
  hist_had_to_part[4]->SetMarkerSize(0.8);
  hist_had_to_part[4]->DrawClone("E1 P SAME");
  //  hist_had_to_part_prph[4]->DrawClone("][ HIST SAME");
  hist_had_to_part_prph[4]->DrawClone("E1 SAME");

  c3->Print("jet_et_had_to_part.eps");
  return 0;
}
