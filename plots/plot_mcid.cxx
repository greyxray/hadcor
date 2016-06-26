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
    

  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_v04b_ari_h1first30partons.root", "read"); 
  //  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_HN_mc_joerg.lt.60.root", "read"); 
  //  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_v04b_ari_joerg.lt.60.root");
  cout << "files attached" << endl;
  
  const Int_t n_hist = 2;
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist1[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el", "h_Q2_el", "h_Zvtx"};
  //  TString s_hist2[n_hist] = {"h_et_LAB_det", "h_eta_LAB_det", "h_et_Breit_det", "h_E_el_bef", "h_Q2_el_bef", "h_Zvtx_bef"};
  //  TString s_hist2[n_hist] = {"h_ana_DetLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_EvtSel_e_el", "h_ana_DetLvl_EvtSel_q2_el", "h_ana_DetLvl_EvtSel_zvtx"};

  TString s_hist1[n_hist] = {"h_had_id", "h_part_id"};
  TString s_hist2[n_hist] = {"h_had_id", "h_part_id"};
  //  TString s_hist2[n_hist] = {"h_ana_DetkLvl_JetSel_EtLabJets", "h_ana_DetLvl_JetSel_EtaLabJets", "h_ana_DetLvl_JetSel_EtBreitJets", "h_ana_DetLvl_JetSel_EtaBreitJets"};

  TString s_legend[n_hist] = {"electron method",
			      "Double-Angle (zufos)"		
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

      //      file2->cd();
      //      cout << "reading "  << i << " ..." << endl;
      //      hist2[i] = (TH1D*)file2->Get(s_hist2[i])->Clone();
      //      hist2[i]->SetName(s_hist2[i] + "2");
      //      cout << "2 - OK, " ;
      //      hist2[i]->SetMarkerColor(kBlack);
      //      hist2[i]->SetMarkerStyle(20);
      //      hist2[i]->SetMarkerSize(0.5);
      //      hist2[i]->SetLineWidth(1);

      //      hist[i]->SetLineStyle(i);
    }
  cout << "histos succesfully read" << endl;

//   cout << "before calling hist_maximum: " << endl;
//   for(Int_t i=0; i<n_hist; i++) {
//     cout << i << " " << hist1[i]->GetMaximum() << endl;
//     Double_t integral1 = hist1[i]->Integral();
//     Double_t integral2 = hist2[i]->Integral();
//     for(Int_t j=0; j<hist2[i]->GetNbinsX(); j++){
//       Double_t bin_old = hist2[i]->GetBinContent(j+1);
//       Double_t binError_old = hist2[i]->GetBinError(j+1);
//       hist2[i]->SetBinContent(j+1, bin_old * integral1/integral2);
//       hist2[i]->SetBinError(j+1, binError_old * integral1/integral2);
//     }
//   }

//  for(Int_t i=0; i<hist2[1]->GetNbinsX()+1; i++) {
//    cout << i << " " << hist2[1]->GetBinContent(i) << " " << hist1[1]->GetBinContent(i) << endl;    
//  }

  //  Double_t y_max = hist_maximum(hist, n_hist);
  TH2D *h_window[n_hist];
  h_window[0] = new TH2D("h_window0", "", 10, 15., 600., 10, 0.1, 10. * hist1[0]->GetMaximum());
  h_window[1] = new TH2D("h_window1", "", 10, 0., 50., 10, 0.1, 10. * hist1[1]->GetMaximum());

  //  h_window[4] = new TH2D("h_window4", "", 10, 0., 150., 10, 3.e3, 1.2 * hist2[4]->GetMaximum());
  //  h_window[5] = new TH2D("h_window5", "", 10, -50., 50., 10, 0.1, 1.2 * hist2[5]->GetMaximum());
  TCanvas* c1;
  c1 = new TCanvas("c1", "c1", 1400, 600);
  c1->Divide(1, 2);
  make_clean_pads(c1, 2, 0, 0);

  sign_window(c1->GetPad(1), h_window[0], "Idfmckin", "Particles", "", "middle");
  sign_window(c1->GetPad(2), h_window[1], "idpart", "Particles", "", "middle");

  //  sign_window(c1->GetPad(5), h_window[4], "Q^{2} (GeV^{2})", "Events", "", "middle");
  //  sign_window(c1->GetPad(6), h_window[5], "Z_{vtx} (cm)", "Events", "", "middle");
  
  TLegend *leg = new TLegend(0.7, 0.7, 0.85, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist1[0], "Orange Ntuples", "l");
  //  leg->AddEntry(hist2[0], "Common Ntuples", "p");
  //  for(Int_t i=0; i<n_hist; i++)
  //    {
  //      leg->AddEntry(hist[i], s_legend[i], "l");
  //    }
  c1->GetPad(1)->cd();
  c1->GetPad(1)->SetLogy();
  c1->GetPad(1)->SetLogx();
  c1->GetPad(2)->SetLogy();
  //  c1->GetPad(5)->SetLogy();
  //  c1->GetPad(1)->SetLogx();
  //  h_window[0]->GetXaxis()->SetRangeUser(0., 35.);

  hist1[0]->GetXaxis()->SetMoreLogLabels();
  //  hist2[0]->GetXaxis()->SetMoreLogLabels();
  for(Int_t i=0; i<n_hist; i++){
    c1->GetPad(i+1)->cd();
    h_window[i]->DrawClone();    
    hist1[i]->DrawClone("E X0 HIST SAME");
    //    hist2[i]->DrawClone("E P X0 SAME");
    if(i==0) {
      leg->Draw();
      hist1[i]->GetXaxis()->SetMoreLogLabels();
    }
  }

  c1->Print("mc_id_hanno_common_first_30_partons.eps");
  //  c1->Print("mc_id_hanno_common.lt.60.C");
  //  c1->Print("Q2.png");
  return 0;
}
