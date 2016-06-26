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
#include <TPaveStats.h>
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

void set_hist_atributes(Int_t number, TString name, TString xaxis, TString yaxis, Double_t xmin, Double_t xmax, Double_t ymin, Bool_t logy, Int_t rebin);
void set_hist_atributes(Int_t number, TString name, TString title, TString xaxis, TString yaxis, Double_t xmin, Double_t xmax, Double_t ymin, Bool_t logy, Int_t rebin);

const Int_t n_hist = 32;
TString s_hist1[n_hist];
TString s_hist2[n_hist];
TString s_xaxis[n_hist];
TString s_yaxis[n_hist];
TString s_title[n_hist];
Double_t hist_xmin[n_hist];
Double_t hist_xmax[n_hist];
Double_t hist_ymin[n_hist];
Bool_t hist_logy[n_hist];
Int_t hist_rebin[n_hist];

int main(int argc, char *argv[])
{
  set_hist_atributes(0, "h_Zvtx", "inclusive DIS", "Z_{vtx} (cm)", "Events", -40, 40, 0.1, kTRUE, 1);
  set_hist_atributes(1, "h_empz", "inclusive DIS", "E-P_{z} (GeV)", "Events", 40, 70, 0.1, kTRUE, 1);
  set_hist_atributes(2, "h_dis_e_en", "Inclusive DIS", "E_{el} (GeV)", "Events", 8., 27., 0.2, kTRUE, 5);
  set_hist_atributes(3, "h_dis_e_en_corr", "Inclusive DIS", "E_{el, corr} (GeV)", "Events", 8., 27., 0.2, kTRUE, 5);
  //  set_hist_atributes(2, "h_E_el", "inclusive DIS", "E_{e, EL} (GeV)", "Events", 8, 26, 0.1, kTRUE, 1);
  //  set_hist_atributes(3, "h_E_da", "inclusive DIS", "E_{e, DA} (GeV)", "Events", 0, 40, 0.1, kTRUE, 1);

  set_hist_atributes(4, "h_y_e_calc", "inclusive DIS", "y_{el}", "Events", 0.1, 0.7, 0.1, kTRUE, 1);
  set_hist_atributes(5, "h_Q2_el", "inclusive DIS", "Q^{2} (GeV^{2})", "Events", 0, 100, 0.1, kTRUE, 1);
  //  set_hist_atributes(6, "h_Q_el_log", "inclusive DIS", "log_{10}Q^{2}", "Events", 0.9, 2.1, 0.1, kTRUE, 1);
  set_hist_atributes(6, "h_cosgamma_had", "inclusive DIS", "cos #gamma_{had}", "Events", -1., 1., 0.1, kTRUE, 1);
  set_hist_atributes(7, "h_Chivtx", "inclusive DIS", "#chi_{vtx}", "Events", -1, 20, 0.1, kTRUE, 1);

  set_hist_atributes(8, "h_cosgamma_had", "inclusive DIS", "cos #gamma_{had}", "Events", -1., 1., 0.1, kTRUE, 1);
  set_hist_atributes(9, "h_cosgamma_had", "inclusive DIS", "cos #gamma_{had}", "Events", -1., 1., 0.1, kTRUE, 1);
  set_hist_atributes(10, "h_cosgamma_had", "inclusive DIS", "cos #gamma_{had}", "Events", -1., 1., 0.1, kTRUE, 1);
  set_hist_atributes(11, "h_cosgamma_had", "inclusive DIS", "cos #gamma_{had}", "Events", -1., 1., 0.1, kTRUE, 1);

  set_hist_atributes(12, "h_nocuts_e_en", "No cuts", "E_{el} (GeV)", "Events", 0., 70., 0.2, kTRUE, 5);
  set_hist_atributes(13, "h_nocuts_e_en_corr", "No cuts", "E_{el, corr} (GeV)", "Events", 0., 70., 0.2, kTRUE, 5);
  set_hist_atributes(14, "h_nocuts_zvtx", "No cuts", "Z_{vtx} (cm)", "Events", -200., 200., 0.2, kTRUE, 5);
  set_hist_atributes(15, "h_nocuts_empz", "No cuts", "E-P_{z}(GeV)", "Events", 0., 150., 0.2, kTRUE, 5);

  set_hist_atributes(16, "h_nocuts_sincand", "No cuts", "Sincand", "Events", 0., 4., 0.2, kTRUE, 1);
  set_hist_atributes(17, "h_nocuts_trk_ntracks", "No cuts", "Trk_ntracks", "Events", 0., 120., 0.2, kTRUE, 1);
  set_hist_atributes(18, "h_nocuts_trk_pt", "No cuts", "P_{T, track} (GeV)", "Events", 0., 400., 0.2, kTRUE, 10);
  set_hist_atributes(19, "h_nocuts_cal_e", "No cuts", "E_{cal} (GeV)", "Events", 0., 500., 0.2, kTRUE, 10);

  set_hist_atributes(20, "h_nocuts_q2", "No cuts", "Q^{2} (GeV^{2})", "Events", 0., 100., 0.2, kTRUE, 1);
  set_hist_atributes(21, "h_nocuts_x", "No cuts", "x_{el} (GeV)", "Events", 1.e-5, 5.e-2, 0.2, kTRUE, 10);
  set_hist_atributes(22, "h_dis_q2", "Inclusive DIS", "Q^{2} (GeV^{2})", "Events", 0., 100., 0.2, kTRUE, 1);
  set_hist_atributes(23, "h_dis_x", "Inclusive DIS", "x_{el} (GeV)", "Events", 1.e-5, 5.e-2, 0.2, kTRUE, 10);

  set_hist_atributes(24, "h_dis_e_en", "Inclusive DIS", "E_{el} (GeV)", "Events", 0., 70., 0.2, kTRUE, 5);
  set_hist_atributes(25, "h_dis_e_en_corr", "Inclusive DIS", "E_{el, corr} (GeV)", "Events", 0., 70., 0.2, kTRUE, 5);
  set_hist_atributes(26, "h_dis_zvtx", "Inclusive DIS", "Z_{vtx} (cm)", "Events", -200., 200., 0.2, kTRUE, 5);
  set_hist_atributes(27, "h_dis_empz", "Inclusive DIS", "E-P_{z}(GeV)", "Events", 0., 150., 0.2, kTRUE, 5);

  set_hist_atributes(28, "h_dis_sincand", "Inclusive DIS", "Sincand", "Events", 0., 4., 0.2, kTRUE, 1);
  set_hist_atributes(29, "h_dis_trk_ntracks", "Inclusive DIS", "Trk_ntracks", "Events", 0., 120., 0.2, kTRUE, 1);
  set_hist_atributes(30, "h_dis_trk_pt", "Inclusive DIS", "P_{T, track} (GeV)", "Events", 0., 400., 0.2, kTRUE, 10);
  set_hist_atributes(31, "h_dis_cal_e", "Inclusive DIS", "E_{cal} (GeV)", "Events", 0., 500., 0.2, kTRUE, 10);



//   set_hist_atributes(0, "h_det_et_LAB", "inclusive jets, det. level", "E_{T,Lab}^{jet}", "Jets", 0, 80, 0.1, kTRUE, 1);
//   set_hist_atributes(1, "h_det_eta_LAB", "inclusive jets, det. level", "#eta_{Lab}^{jet}", "Jets", -1, 3, 0.1, kFALSE, 1);
//   set_hist_atributes(2, "h_det_et_Breit", "inclusive jets, det. level", "E_{T, Breit} (GeV)", "Jets", 0., 80., 0.1, kTRUE, 1);
//   set_hist_atributes(3, "h_det_eta_Breit", "inclusive jets, det. level", "#eta_{Breit}", "Jets", -1., 3., 0.1, kFALSE, 1);
//   set_hist_atributes(4, "h_E_el", "inclusive DIS, det. level", "E_{e} (GeV)", "Events",  8., 30., 0.1, kFALSE, 1);
//   set_hist_atributes(5, "h_y_e_calc", "inclusive DIS, det. level", "y_{el}", "Events", 0.1, 0.7, 0.1, kFALSE, 1);
//   set_hist_atributes(6, "h_Q2_el", "inclusive DIS, det. level", "Q^{2}_{el} (GeV^{2})", "Events", 0., 110., 0.1, kTRUE, 1);
//   set_hist_atributes(7, "h_Q_el_log", "inclusive DIS, det.level",  "log_{10} Q^{2}_{el}", "Events", 0.8, 2.2, 0.1, kFALSE, 1);
//   set_hist_atributes(8, "h_had_px_sum_short", "without cuts, had. level", "P_{x, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(9, "h_had_py_sum_short", "without cuts, had. level", "P_{y, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(10, "h_had_pz_sum", "without cuts, had. level", "P_{z, sum}^{had} (GeV)", "Events", -1000., 1000., 0.1, kTRUE, 1);
//   set_hist_atributes(11, "h_had_e_sum", "without cuts, had. level", "E_{sum}^{had} (GeV)", "Events", -1000., 1000., 0.1, kTRUE, 1);
//   set_hist_atributes(12, "h_had_px", "without cuts, had. level", "P_{x}^{had} (GeV)", "Particles", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(13, "h_had_py", "without cuts, had. leel",  "P_{y}^{had} (GeV)", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(14, "h_had_pz", "without cuts, had. level", "P_{z}^{had} (GeV)", "Particles", -50., 1000., 0.1, kTRUE, 1);
//   set_hist_atributes(15, "h_had_e", "without cuts, had. level", "E^{had} (GeV)", "Particles", -50., 1000., 0.1, kTRUE, 1);
//   set_hist_atributes(16, "h_part_px_sum_short", "without cuts, part. level", "P_{x, sum}^{part} (GeV)", "Events", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(17, "h_part_py_sum_short", "without cuts, part. level", "P_{y, sum}^{part} (GeV)", "Events", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(18, "h_part_pz_sum", "without cuts, part. level", "P_{z, sum}^{part} (GeV)", "Events", -1000., 300., 0.1, kTRUE, 1);
//   set_hist_atributes(19, "h_part_e_sum", "without cuts, part. level", "E_{sum}^{part} (GeV)", "Events", -1000., 300., 0.1, kTRUE, 1);
//   set_hist_atributes(20, "h_part_px", "without cuts, part. level", "P_{x}^{part} (GeV)", "Particles", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(21, "h_part_py", "without cuts, part. level", "P_{y}^{part} (GeV)", "Particles", -100., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(22, "h_part_pz", "without cuts, part. level", "P_{z}^{part} (GeV)", "Particles", -50., 1000., 0.1, kTRUE, 1);
//   set_hist_atributes(23, "h_part_e", "without cuts, part. level", "E^{part} (GeV)", "Particles", -50., 1000., 0.1, kTRUE, 1);

//   //hadron level
//   set_hist_atributes(24, "h_had_dis_px_sum", "inclusive DIS, had. level", "P_{x, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(25, "h_had_dis_py_sum", "inclusive DIS, had. level", "P_{y, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(26, "h_had_dis_pz_sum", "inclusive DIS, had. level", "P_{z, sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(27, "h_had_dis_e_sum", "inclusive DIS, had. level", "E_{sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(28, "h_had_ijet_px_sum", "inclusive jets, had. level", "P_{x, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(29, "h_had_ijet_py_sum", "inclusive jets, had. level", "P_{y, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(30, "h_had_ijet_pz_sum", "inclusive jets, had. level", "P_{z, sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(31, "h_had_ijet_e_sum", "inclusive jets, had. level", "E_{sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(32, "h_had_idijet_px_sum", "inclusive dijets, had. level", "P_{x, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(33, "h_had_idijet_py_sum", "inclusive dijets, had. level", "P_{y, sum}^{had} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(34, "h_had_idijet_pz_sum", "inclusive dijets, had. level", "P_{z, sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(35, "h_had_idijet_e_sum", "inclusive dijets, had. level", "E_{sum}^{had} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(36, "h_had_idijet_mjj", "inclusive dijets, had. level", "M_{jj} (GeV)", "Events", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(37, "h_had_idijet_pt_breit_mean", "inclusive dijets, had.level", "#bar{P_{T,Breit}} (GeV)", "Events", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(38, "h_had_idijet_eta_prime", "inclusive dijets, had.level", "#eta\'", "Events", 0., 2., 0.1, kTRUE, 1);
//   set_hist_atributes(39, "h_had_idijet_ksi", "inclusive dijets, had.level", "log_{10}#xi", "Events", -2.2, 0., 0.1, kTRUE, 1);

//   set_hist_atributes(40, "h_had_ijet_q2", "inclusive jets, had. level", "Q^{2} (GeV^{2})", "Jets", 0., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(41, "h_had_ijet_pt_breit", "inclusive jets, had.level", "P_{T,Breit} (GeV)", "Jets", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(42, "h_had_ijet_eta_lab", "inclusive jets, had.level", "#eta_{Lab}", "Jets", -1., 2.5, 0.1, kFALSE, 1);
//   set_hist_atributes(43, "h_had_ijet_pt_lab", "inclusive jets, had.level", "P_{T,Lab} (GeV)", "Jets", 0., 120., 0.1, kTRUE, 1);

//   //parton level
//   set_hist_atributes(44, "h_part_dis_px_sum", "inclusive DIS, part. level", "P_{x, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(45, "h_part_dis_py_sum", "inclusive DIS, part. level", "P_{y, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(46, "h_part_dis_pz_sum", "inclusive DIS, part. level", "P_{z, sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(47, "h_part_dis_e_sum", "inclusive DIS, part. level", "E_{sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(48, "h_part_ijet_px_sum", "inclusive jets, part. level", "P_{x, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(49, "h_part_ijet_py_sum", "inclusive jets, part. level", "P_{y, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(50, "h_part_ijet_pz_sum", "inclusive jets, part. level", "P_{z, sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(51, "h_part_ijet_e_sum", "inclusive jets, part. level", "E_{sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(52, "h_part_idijet_px_sum", "inclusive dijets, part. level", "P_{x, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(53, "h_part_idijet_py_sum", "inclusive dijets, part. level", "P_{y, sum}^{part} (GeV)", "Events", -10., 10., 0.1, kTRUE, 1);
//   set_hist_atributes(54, "h_part_idijet_pz_sum", "inclusive dijets, part. level", "P_{z, sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);
//   set_hist_atributes(55, "h_part_idijet_e_sum", "inclusive dijets, part. level", "E_{sum}^{part} (GeV)", "Events", -1000., 100., 0.1, kTRUE, 1);

//   set_hist_atributes(56, "h_part_idijet_mjj", "inclusive dijets, part. level", "M_{jj} (GeV)", "Events", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(57, "h_part_idijet_pt_breit_mean", "inclusive dijets, part.level", "#bar{P_{T,Breit}} (GeV)", "Events", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(58, "h_part_idijet_eta_prime", "inclusive dijets, part.level", "#eta\'", "Events", 0., 2., 0.1, kTRUE, 1);
//   set_hist_atributes(59, "h_part_idijet_ksi", "inclusive dijets, part.level", "log_{10}#xi", "Events", -2.2, 0., 0.1, kTRUE, 1);

//   set_hist_atributes(60, "h_part_ijet_q2", "inclusive jets, part. level", "Q^{2} (GeV^{2})", "Jets", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(61, "h_part_ijet_pt_breit", "inclusive jets, part.level", "P_{T,Breit} (GeV)", "Jets", 0., 120., 0.1, kTRUE, 1);
//   set_hist_atributes(62, "h_part_ijet_eta_lab", "inclusive jets, part.level", "#eta_{Lab}", "Jets", -1., 2.5, 0.1, kFALSE, 1);
//   set_hist_atributes(63, "h_part_ijet_pt_lab", "inclusive jets, part.level", "P_{T,Lab} (GeV)", "Jets", 0., 120., 0.1, kTRUE, 1);

//   set_hist_atributes(64, "h_part_dis_q2", "inclusive DIS, part. level", "Q^{2} (GeV^{2})", "Events", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(65, "h_part_dis_x", "inclusive DIS, part. level", "x", "Events", 1.e-5, 5.e-3, 0.1, kTRUE, 1);
//   set_hist_atributes(66, "h_had_dis_q2", "inclusive DIS, had. level", "Q^{2} (GeV^{2})", "Events", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(67, "h_had_dis_x", "inclusive DIS, had. level", "x", "Events", 1.e-5, 5.e-3, 0.1, kTRUE, 1);

//   set_hist_atributes(68, "h_part_ijet_q2", "inclusive jets, part. level", "Q^{2} (GeV^{2})", "Jets", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(69, "h_part_ijet_x", "inclusive jets, part. level", "x", "Jets", 1.e-5, 5.e-3, 0.1, kTRUE, 1);
//   set_hist_atributes(70, "h_had_ijet_q2", "inclusive jets, had. level", "Q^{2} (GeV^{2})", "Jets", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(71, "h_had_ijet_x", "inclusive jets, had. level", "x", "Jets", 1.e-5, 5.e-3, 0.1, kTRUE, 1);

//   set_hist_atributes(72, "h_part_idijet_q2", "inclusive dijets, part. level", "Q^{2} (GeV^{2})", "Events", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(73, "h_part_idijet_x", "inclusive dijets, part. level", "x", "Events", 1.e-5, 5.e-3, 0.1, kTRUE, 1);
//   set_hist_atributes(74, "h_had_idijet_q2", "inclusive dijets, had. level", "Q^{2} (GeV^{2})", "Events", 0., 100., 1, kTRUE, 1);
//   set_hist_atributes(75, "h_had_idijet_x", "inclusive dijets, had. level", "x", "Events", 1.e-5, 5.e-3, 0.1, kTRUE, 1);

  //  set_hist_atributes(28, "h_had_ijet_px_sum", "inclusive jets", "P_{x, sum}^{had} (GeV)", "Particles", -1000., 100., 0.1, kTRUE, 1);
//  "h_det_et_LAB", "h_det_eta_LAB", "h_det_et_Breit", "h_det_eta_Breit",
// 			     "h_E_el", "h_y_e_calc", "h_Q2_el", "h_Q_el_log",
// 			     "h_had_px_sum_short", "h_had_py_sum_short", "h_had_pz_sum", "h_had_e_sum",
// 			     "h_had_px", "h_had_py", "h_had_pz", "h_had_e",
// 			     "h_part_px_sum_short", "h_part_py_sum_short", "h_part_pz_sum", "h_part_e_sum",
// 			     "h_part_px", "h_part_py", "h_part_pz", "h_part_e"  
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleX(0.15);
    

  //  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/60211h1.root", "read"); 
  //  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/60211orange.root");
  TFile *file1 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_data_h1.root", "read"); 
  TFile *file2 = new TFile("/data/zenith234b/kuprash/analyses/jets_low_q2/0607p_data_orange.root");
  cout << "files attached" << endl;
  

  TH1D *hist1[n_hist];
  TH1D *hist2[n_hist];
  TH1D *hist_ratio[n_hist];
  
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
      
      hist1[i]->Rebin(hist_rebin[i]);
      hist2[i]->Rebin(hist_rebin[i]);
      //      hist2[i]->SetLineWidth(1);
      //      hist[i]->SetLineStyle(i);
    }
  cout << "histos succesfully read" << endl;

//   cout << "before calling hist_maximum: " << endl;
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
    hist_ratio[i] = (TH1D*)hist2[i]->Clone();
    hist_ratio[i]->SetName(s_hist2[i] + "_ratio");
    hist_ratio[i]->Divide(hist1[i]);
  }


  //  Double_t y_max = hist_maximum(hist, n_hist);

  TH2D *h_window[n_hist];
  TH2D *h_window_ratio[n_hist];
  for(Int_t i=0; i<n_hist; i++) {
    Double_t offset = 1.2;
    if(hist_logy[i])
      offset = 12.;
    if(i%4==0)
      offset = 1.5;
    if(i%4==0 && hist_logy[i])
      offset = 1.2e3;
	
    h_window[i] = new TH2D("h_window" + s_hist1[i], "title", 10, hist_xmin[i], hist_xmax[i], 10, hist_ymin[i], offset * hist2[i]->GetMaximum());
    Double_t ratio_min = 0.8;
    Double_t ratio_max = 1.2;
    if(s_hist1[i] == "h_Chivtx") {ratio_min = 0.; ratio_max = 2.;}
    if(s_hist1[i] == "h_dis_sincand") {ratio_min = 0.; ratio_max = 2.;}
    h_window_ratio[i] = new TH2D("h_window" + s_hist1[i] + "_ratio", "", 10, hist_xmin[i], hist_xmax[i], 10, ratio_min, ratio_max);
  }

  //  h_window[4] = new TH2D("h_window4", "", 10, ., 150., 10, 3.e3, 1.2 * hist2[4]->GetMaximum());
  //  h_window[5] = new TH2D("h_window5", "", 10, -50., 50., 10, 0.1, 1.2 * hist2[5]->GetMaximum());
  TCanvas* c1[n_hist/4+1];
  for(Int_t i=0; i<n_hist/4+1; i++)
    {
      TString s; s.Form("c_%i", i);
      c1[i] = new TCanvas(s, s, 800, 600);
      c1[i]->Divide(2, 2);
      make_clean_pads(c1[i], 4, 0, 0);
      for(Int_t j=0; j<4; j++) {
	prepare_canvas_for_cross_sec(c1[i]->GetPad(j+1));
      }
    }

  for(Int_t i=0; i<n_hist; i++) {
    sign_window(c1[i/4]->GetPad(i%4+1)->GetPad(1), h_window[i], s_xaxis[i], s_yaxis[i], s_title[i], "large");
    sign_window(c1[i/4]->GetPad(i%4+1)->GetPad(2), h_window_ratio[i], s_xaxis[i], "ratio", "", "large");
    c1[i/4]->GetPad(i%4+1)->GetPad(1)->SetLogy(hist_logy[i]);
  }
  
  TLegend *leg = new TLegend(0.2, 0.75, 0.65, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist1[0], "Orange Ntuples", "l");
  leg->AddEntry(hist2[0], "Common Ntuples", "p");
  //  for(Int_t i=0; i<n_hist; i++)
  //    {
  //      leg->AddEntry(hist[i], s_legend[i], "l");
  //    }
  //  c1[0]->GetPad(1)->cd();
//   c1[0]->GetPad(1)->SetLogy();
//   c1[0]->GetPad(3)->SetLogy();

//   c1[2]->GetPad(1)->SetLogy();
//   c1[2]->GetPad(2)->SetLogy();
//   c1[2]->GetPad(3)->SetLogy();
//   c1[2]->GetPad(4)->SetLogy();

//   c1[3]->GetPad(1)->SetLogy();
//   c1[3]->GetPad(2)->SetLogy();
//   c1[3]->GetPad(3)->SetLogy();
//   c1[3]->GetPad(4)->SetLogy();

//   c1[4]->GetPad(1)->SetLogy();
//   c1[4]->GetPad(2)->SetLogy();
//   c1[4]->GetPad(3)->SetLogy();
//   c1[4]->GetPad(4)->SetLogy();

//   c1[5]->GetPad(1)->SetLogy();
//   c1[5]->GetPad(2)->SetLogy();
//   c1[5]->GetPad(3)->SetLogy();
//   c1[5]->GetPad(4)->SetLogy();
  //  c1->GetPad(5)->SetLogy();
  //  c1->GetPad(1)->SetLogx();
  //  h_window[0]->GetXaxis()->SetRangeUser(0., 35.);
  TF1 *f = new TF1("f1", "1", -10000., 10000.);
  f->SetLineColor(kBlue);
  f->SetLineWidth(1);

  for(Int_t i=0; i<n_hist; i++){
    c1[i/4]->GetPad(i%4+1)->GetPad(1)->cd();
    h_window[i]->DrawClone();    
    hist1[i]->DrawClone("E X0 HIST SAME");
    //    hist2[i]->SetStats(kTRUE);
    //    hist2[i]->ResetStats();
    hist2[i]->DrawClone("E P X0 SAMES");
    //    gPad->Update();
    //    TPaveStats *st = (TPaveStats*)hist2[i]->FindObject("stats");
    //    st->Draw();
    if(i%4==0)
      leg->DrawClone();
    c1[i/4]->GetPad(i%4+1)->GetPad(2)->cd();
    h_window_ratio[i]->GetYaxis()->SetTitleOffset(0.4);
    h_window_ratio[i]->DrawClone();    
    f->DrawClone("SAME");
    hist_ratio[i]->DrawClone(" SAME");
  }
  for(Int_t i=0; i<n_hist/4+1; i++) {
    TString s; s.Form("data_hanno_common%i.eps", i);
    c1[i]->Print(s);
  }
  //  c1->Print("Q2.png");
  return 0;
}

void set_hist_atributes(Int_t number, TString name, TString xaxis, TString yaxis, Double_t xmin, Double_t xmax, Double_t ymin, Bool_t logy, Int_t rebin)
{
  s_hist1[number] = name;
  s_hist2[number] = name;
  s_xaxis[number] = xaxis;
  s_yaxis[number] = yaxis;
  s_title[number] = "";
  hist_xmin[number] = xmin;
  hist_xmax[number] = xmax;
  hist_ymin[number] = ymin;
  hist_logy[number] = logy;
  hist_rebin[number] = rebin;
}

void set_hist_atributes(Int_t number, TString name, TString title, TString xaxis, TString yaxis, Double_t xmin, Double_t xmax, Double_t ymin, Bool_t logy, Int_t rebin)
{
  s_hist1[number] = name;
  s_hist2[number] = name;
  s_xaxis[number] = xaxis;
  s_yaxis[number] = yaxis;
  s_title[number] = title;
  hist_xmin[number] = xmin;
  hist_xmax[number] = xmax;
  hist_ymin[number] = ymin;
  hist_logy[number] = logy;
  hist_rebin[number] = rebin;
}
