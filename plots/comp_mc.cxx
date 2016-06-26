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

  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_HN_mc_joerg.root", "read"); 
  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_v04b_ari_joerg.root");
  cout << "files attached" << endl;
  
  const Int_t n_hist = 24;
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el", "h_Q2_el", "h_Zvtx"};
  //  TString s_hist2[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist2[n_hist] = {"h_ana_DetLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_EvtSel_e_el", "h_ana_DetLvl_EvtSel_q2_el", "h_ana_DetLvl_EvtSel_zvtx"};

  TString s_hist1[n_hist] = {"h_det_et_LAB", "h_det_eta_LAB", "h_det_et_Breit", "h_det_eta_Breit",
			     "h_E_el", "h_y_e_calc", "h_Q2_el", "h_Q_el_log",
			     "h_had_px_sum_short", "h_had_py_sum_short", "h_had_pz_sum", "h_had_e_sum",
			     "h_had_px", "h_had_py", "h_had_pz", "h_had_e",
			     "h_part_px_sum_short", "h_part_py_sum_short", "h_part_pz_sum", "h_part_e_sum",
			     "h_part_px", "h_part_py", "h_part_pz", "h_part_e"  };
  TString s_hist2[n_hist] = {"h_det_et_LAB", "h_det_eta_LAB", "h_det_et_Breit", "h_det_eta_Breit",
			     "h_E_el", "h_y_e_calc", "h_Q2_el", "h_Q_el_log",
			     "h_had_px_sum_short", "h_had_py_sum_short", "h_had_pz_sum", "h_had_e_sum", 
			     "h_had_px", "h_had_py", "h_had_pz", "h_had_e",
			     "h_part_px_sum_short", "h_part_py_sum_short", "h_part_pz_sum", "h_part_e_sum",
			     "h_part_px", "h_part_py", "h_part_pz", "h_part_e"  };
  //  TString s_hist2[n_hist] = {"h_ana_DetkLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_JetSel_EtaBreitJets"};

  TString s_legend[n_hist] = {"electron method",
			      "Double-Angle (zufos)",
			      "Double-Angle (cells)",
			      "Jacquet-Blondel (zufos)"//,
			      //			      "Jacquet-Blondel (cells)",
			      //			      "Q^{2} at generated level"
  };
  TH1D *hist1[n_hist];
  TH1D *hist2[n_hist];

  for(Int_t i=0; i<n_hist; i++)
    {
      file1->cd();
      cout << "reading "  << i << " ..." << endl;
      hist1[i] = (TH1D*)file1->Get(s_hist1[i])->Clone();
      hist1[i]->SetName(s_hist1[i] + "1");
      cout << "1 - OK, " ;
      hist1[i]->SetLineColor(kBlue);
      hist1[i]->SetLineWidth(1);

      file2->cd();
      cout << "reading "  << i << " ..." << endl;
      hist2[i] = (TH1D*)file2->Get(s_hist2[i])->Clone();
      hist2[i]->SetName(s_hist2[i] + "2");
      cout << "2 - OK, " ;
      hist2[i]->SetMarkerColor(kBlack);
      hist2[i]->SetMarkerStyle(20);
      hist2[i]->SetMarkerSize(0.5);
      //      hist2[i]->SetLineWidth(1);

      //      hist[i]->SetLineStyle(i);
    }
  cout << "histos succesfully read" << endl;

  cout << "before calling hist_maximum: " << endl;
  for(Int_t i=0; i<n_hist; i++) {
    cout << i << " " << hist1[i]->GetMaximum() << endl;
    Double_t integral1 = hist1[i]->Integral();
    Double_t integral2 = hist2[i]->Integral();
    for(Int_t j=0; j<hist2[i]->GetNbinsX(); j++){
      Double_t bin_old = hist2[i]->GetBinContent(j+1);
      Double_t binError_old = hist2[i]->GetBinError(j+1);
      hist2[i]->SetBinContent(j+1, bin_old * integral1/integral2);
      hist2[i]->SetBinError(j+1, binError_old * integral1/integral2);
    }
      
    //    hist2[i]->Scale();
  }
  //  Double_t y_max = hist_maximum(hist, n_hist);
  TH2D *h_window[n_hist];
  h_window[0] = new TH2D("h_window0", "", 10, 0., 80., 10, 0.1, 1.e7/*2 * hist2[0]->GetMaximum()*/);
  h_window[1] = new TH2D("h_window1", "", 10, -1., 3., 10, 0.1, 1.2 * hist2[1]->GetMaximum());
  h_window[2] = new TH2D("h_window2", "", 10, 0., 80., 10, 0.1, 1.2 * hist2[2]->GetMaximum());
  h_window[3] = new TH2D("h_window3", "", 10, -1., 3., 10, 0.1, 1.2 * hist2[3]->GetMaximum());

  h_window[4] = new TH2D("h_window4", "", 10, 8., 30., 10, 0.1, 1.5 * hist2[4]->GetMaximum());
  h_window[5] = new TH2D("h_window5", "", 10, 0.1, 0.7, 10, 0.1, 1.2 * hist2[5]->GetMaximum());
  h_window[6] = new TH2D("h_window6", "", 10, 0., 110., 10, 0.1, 1.2 * hist2[6]->GetMaximum());
  h_window[7] = new TH2D("h_window7", "", 10, 0.8, 2.2, 10, 0.1, 1.2 * hist2[7]->GetMaximum());

  h_window[8] = new TH2D("h_window8", "", 10, -10., 10., 10, 0.1, 1.2e4 * hist2[8]->GetMaximum());
  h_window[9] = new TH2D("h_window9", "", 10, -10., 10., 10, 0.1, 12. * hist2[9]->GetMaximum());
  h_window[10] = new TH2D("h_window10", "", 10, -1000., 1000., 10, 0.1, 12. * hist2[10]->GetMaximum());
  h_window[11] = new TH2D("h_window11", "", 10, -1000., 1000., 10, 0.1, 12. * hist2[11]->GetMaximum());

  h_window[12] = new TH2D("h_window12", "", 10, -100., 100., 10, 0.1, 1.2e4 * hist2[12]->GetMaximum());
  h_window[13] = new TH2D("h_window13", "", 10, -100., 100., 10, 0.1, 12. * hist2[13]->GetMaximum());
  h_window[14] = new TH2D("h_window14", "", 10, -50., 1000., 10, 0.1, 12. * hist2[14]->GetMaximum());
  h_window[15] = new TH2D("h_window15", "", 10, -50., 1000., 10, 0.1, 12. * hist2[15]->GetMaximum());

  h_window[16] = new TH2D("h_window16", "", 10, -100., 100., 10, 0.1, 1.2e4 * hist2[16]->GetMaximum());
  h_window[17] = new TH2D("h_window17", "", 10, -100., 100., 10, 0.1, 12. * hist2[17]->GetMaximum());
  h_window[18] = new TH2D("h_window18", "", 10, -1000., 300., 10, 0.1, 12. * hist2[18]->GetMaximum());
  h_window[19] = new TH2D("h_window19", "", 10, -1000., 300., 10, 0.1, 12. * hist2[19]->GetMaximum());

  h_window[20] = new TH2D("h_window20", "", 10, -100., 100., 10, 0.1, 1.2e4 * hist2[20]->GetMaximum());
  h_window[21] = new TH2D("h_window21", "", 10, -100., 100., 10, 0.1, 12. * hist2[21]->GetMaximum());
  h_window[22] = new TH2D("h_window22", "", 10, -50., 1000., 10, 0.1, 12. * hist2[22]->GetMaximum());
  h_window[23] = new TH2D("h_window23", "", 10, -50., 1000., 10, 0.1, 12. * hist2[23]->GetMaximum());
  //  h_window[4] = new TH2D("h_window4", "", 10, ., 150., 10, 3.e3, 1.2 * hist2[4]->GetMaximum());
  //  h_window[5] = new TH2D("h_window5", "", 10, -50., 50., 10, 0.1, 1.2 * hist2[5]->GetMaximum());
  TCanvas* c1[n_hist/4+1];
  for(Int_t i=0; i<n_hist/4+1; i++)
    {
      TString s; s.Form("c_%i", i);
      c1[i] = new TCanvas(s, s, 800, 600);
      c1[i]->Divide(2, 2);
      make_clean_pads(c1[i], 4, 0, 0);
    }

  sign_window(c1[0]->GetPad(1), h_window[0], "E_{T, Lab} (GeV)", "Events", "", "middle");
  sign_window(c1[0]->GetPad(2), h_window[1], "#eta_{Lab}", "Events", "", "middle");
  sign_window(c1[0]->GetPad(3), h_window[2], "E_{T, Breit} (GeV)", "Events", "", "middle");
  sign_window(c1[0]->GetPad(4), h_window[3], "#eta_{Breit}", "Events", "", "middle");

  sign_window(c1[1]->GetPad(1), h_window[4], "E_{e} (GeV)", "Events", "", "middle");
  sign_window(c1[1]->GetPad(2), h_window[5], "y_{el}", "Events", "", "middle");
  sign_window(c1[1]->GetPad(3), h_window[6], "Q^{2}_{el} (GeV^{2})", "Events", "", "middle");
  sign_window(c1[1]->GetPad(4), h_window[7], "log_{10} Q^{2}_{el}", "Events", "", "middle");

  sign_window(c1[2]->GetPad(1), h_window[8], "P_{x, sum}^{had} (GeV)", "Events", "", "middle");
  sign_window(c1[2]->GetPad(2), h_window[9], "P_{y, sum}^{had} (GeV)", "Events", "", "middle");
  sign_window(c1[2]->GetPad(3), h_window[10], "P_{z, sum}^{had} (GeV)", "Events", "", "middle");
  sign_window(c1[2]->GetPad(4), h_window[11], "E_{sum}^{had} (GeV)", "Events", "", "middle");

  sign_window(c1[3]->GetPad(1), h_window[12], "P_{x}^{had} (GeV)", "Particles", "", "middle");
  sign_window(c1[3]->GetPad(2), h_window[13], "P_{y}^{had} (GeV)", "Particles", "", "middle");
  sign_window(c1[3]->GetPad(3), h_window[14], "P_{z}^{had} (GeV)", "Particles", "", "middle");
  sign_window(c1[3]->GetPad(4), h_window[15], "E^{had} (GeV)", "Particles", "", "middle");

  sign_window(c1[4]->GetPad(1), h_window[16], "P_{x, sum}^{part} (GeV)", "Events", "", "middle");
  sign_window(c1[4]->GetPad(2), h_window[17], "P_{y, sum}^{part} (GeV)", "Events", "", "middle");
  sign_window(c1[4]->GetPad(3), h_window[18], "P_{z, sum}^{part} (GeV)", "Events", "", "middle");
  sign_window(c1[4]->GetPad(4), h_window[19], "E_{sum}^{part} (GeV)", "Events", "", "middle");

  sign_window(c1[5]->GetPad(1), h_window[20], "P_{x}^{part} (GeV)", "Particles", "", "middle");
  sign_window(c1[5]->GetPad(2), h_window[21], "P_{y}^{part} (GeV)", "Particles", "", "middle");
  sign_window(c1[5]->GetPad(3), h_window[22], "P_{z}^{part} (GeV)", "Particles", "", "middle");
  sign_window(c1[5]->GetPad(4), h_window[23], "E^{part} (GeV)", "Particles", "", "middle");
  //  sign_window(c1->GetPad(5), h_window[4], "Q^{2} (GeV^{2})", "Events", "", "middle");
  //  sign_window(c1->GetPad(6), h_window[5], "Z_{vtx} (cm)", "Events", "", "middle");
  
  TLegend *leg = new TLegend(0.2, 0.7, 0.85, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist1[0], "Orange Ntuples", "l");
  leg->AddEntry(hist2[0], "Common Ntuples", "p");
  //  for(Int_t i=0; i<n_hist; i++)
  //    {
  //      leg->AddEntry(hist[i], s_legend[i], "l");
  //    }
  //  c1[0]->GetPad(1)->cd();
  c1[0]->GetPad(1)->SetLogy();
  c1[0]->GetPad(3)->SetLogy();

  c1[2]->GetPad(1)->SetLogy();
  c1[2]->GetPad(2)->SetLogy();
  c1[2]->GetPad(3)->SetLogy();
  c1[2]->GetPad(4)->SetLogy();

  c1[3]->GetPad(1)->SetLogy();
  c1[3]->GetPad(2)->SetLogy();
  c1[3]->GetPad(3)->SetLogy();
  c1[3]->GetPad(4)->SetLogy();

  c1[4]->GetPad(1)->SetLogy();
  c1[4]->GetPad(2)->SetLogy();
  c1[4]->GetPad(3)->SetLogy();
  c1[4]->GetPad(4)->SetLogy();

  c1[5]->GetPad(1)->SetLogy();
  c1[5]->GetPad(2)->SetLogy();
  c1[5]->GetPad(3)->SetLogy();
  c1[5]->GetPad(4)->SetLogy();
  //  c1->GetPad(5)->SetLogy();
  //  c1->GetPad(1)->SetLogx();
  //  h_window[0]->GetXaxis()->SetRangeUser(0., 35.);

  for(Int_t i=0; i<n_hist; i++){
    c1[i/4]->GetPad(i%4+1)->cd();
    h_window[i]->DrawClone();    
    hist1[i]->DrawClone("E X0 HIST SAME");
    hist2[i]->DrawClone("E P X0 SAME");
    if(i%4==0)
      leg->DrawClone();
  }
  for(Int_t i=0; i<n_hist/4+1; i++) {
    TString s; s.Form("mc_hanno_common%i.eps", i);
    c1[i]->Print(s);
  }
  //  c1->Print("Q2.png");
  return 0;
}
