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

  TFile *file1 = new TFile("../mc_bg_rad0405e_parton_all_copy.root", "read"); 
  //  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_v04b_RD.root");
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
  for(Int_t i=0; i<n_hist; i++) {
    s_hist_had[i] = "Cross_Sections_Histograms/h_had_cross_" + s_var[i];
    s_hist_part[i] = "Cross_Sections_Histograms/h_part_cross_" + s_var[i];
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
  TH1D *hist_had_to_part[n_hist];
  for(Int_t i=0; i<n_hist; i++)
    {
      file1->cd();
      cout << "reading "  << i << " ..." << endl;
      hist_had[i] = (TH1D*)file1->Get(s_hist_had[i])->Clone();
      cout << "1 - OK, " ;
      cout << "reading "  << i << " ..." << endl;
      hist_part[i] = (TH1D*)file1->Get(s_hist_part[i])->Clone();
      cout << "2 - OK, " ;
      hist_had_to_part[i] = (TH1D*)hist_had[i]->Clone();
      hist_had_to_part[i]->SetTitle((TString)hist_had[i]->GetTitle() + "_ratio");
      hist_had_to_part[i]->Divide(hist_part[i]);
      hist_had_to_part[i]->SetLineWidth(2);
      hist_had_to_part[i]->SetLineColor(kBlue);
    }
  cout << "histos succesfully read" << endl;

  //
  // calculate uncertainties of the ratio according to the Gavin's McCance's thesis (appendix A) - ratio of the correlated histograms
  //
  // C = hadron && parton
  // D = hadron && !parton
  // E = parton && !hadron

  
  cout << "before calling hist_maximum: " << endl;
  TH2D *h_window[n_hist];
  h_window[0] = new TH2D("h_window_et", "", number_etbins, et_bin, 10, 0.1, 1.25);
  h_window[1] = new TH2D("h_window_eta", "", number_etabins, eta_bin_crosssec, 10, 0.1, 1.25);
  h_window[2] = new TH2D("h_window_Q2", "", number_Q2bins, Q2_bin, 10, 0.1, 1.25);
  h_window[3] = new TH2D("h_window_x", "", number_xbins, x_bin, 10, 0.1, 1.25);
  h_window[4] = new TH2D("h_window_et_jet", "", number_et_jetbins, et_jet_bin, 10, 0.1, 1.25);
  h_window[5] = new TH2D("h_window_eta_jet", "", number_eta_jetbins, eta_jet_bin, 10, 0.1, 1.25);
  TCanvas* c1;
  c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(3, 2);
  make_clean_pads(c1, 6, 0, 0);
  for(Int_t i=0; i<n_hist; i++)
    sign_window(c1->GetPad(i+1), h_window[i], s_dim[i], "C_{hadronisation} = N_{had.} / N_{part.}", "", "middle");
  
  TLegend *leg = new TLegend(0.2, 0.7, 0.85, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  //  leg->AddEntry(hist1[0], "Orange Ntuples", "l");
  //  leg->AddEntry(hist2[0], "Common Ntuples", "p");
  //  for(Int_t i=0; i<n_hist; i++)
  //    {
  //      leg->AddEntry(hist[i], s_legend[i], "l");
  //    }
  //  c1->GetPad(1)->cd();
  c1->GetPad(3)->SetLogx();
  c1->GetPad(4)->SetLogx();
  c1->GetPad(5)->SetLogx();
  TF1 *f_unity = new TF1("unity", "1", -10000., 10000.);
  f_unity->SetLineColor(34);
  for(Int_t i=0; i<n_hist; i++){
    c1->GetPad(i+1)->cd();
    h_window[i]->DrawClone();    
    f_unity->DrawClone("SAME");
    hist_had_to_part[i]->DrawClone("E X0 HIST SAME");
    //    hist2[i]->DrawClone("E P X0 SAME");
    //    if(i==0)
    //      leg->Draw();
  }

  c1->Print("hadronisation_corrections.eps");
  //  c1->Print("Q2.png");
  return 0;
}
