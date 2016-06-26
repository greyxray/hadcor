#define selectorBegin_c
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
//#include "selector.h"
//#ifndef constants_h
#include "constants.h"
//#endif
#include "selector.h"
#include "runinfo.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TSelector.h>
//#include <TLorentzVector.h>
//#include <TObjArray.h>
//#include "Dijets.h"
//#include "KtJet/KtEvent.h"
//#include "KtJet/KtLorentzVector.h"
/*using KtJet::KtLorentzVector;
using KtJet::KtEvent;

extern Bool_t excludedrunlist(Int_t);
extern Bool_t take_run(Int_t);

extern Bool_t poltake(Int_t, Int_t, Int_t, Int_t, RunInfo runinfo[], const int);
extern Double_t polcorr0405e(Double_t );
extern Double_t polcorr0607p(Double_t );
extern Double_t polcorr06e(Double_t );
extern Double_t reweighting(Double_t , TString);
extern void fill_trigger_bits(TH1D* h_dis, Int_t r_Fltw[], Int_t n);
extern Bool_t flt_bit_is_1(UInt_t bit, Int_t r_Fltw[], Int_t n);
*/

/*  extern Double_t        C_scale[number_of_eta_bins];
  extern Double_t        a0_corr[number_of_eta_bins];
  extern Double_t        a1_corr[number_of_eta_bins];*/

void selector::Begin(/*TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  //   TString option = GetOption();
  time_begin.Set();
  cout << "begin time: " << time_begin.AsString() << endl;
  check_cuts = kFALSE;//kTRUE;
  use_ktjetb_prph = kTRUE;//kFALSE;
  root_file_name = "";
  if(Data)
    root_file_name += "data";
  if(!Data)
    root_file_name += mc_type;

  //  if(mc_type != "mc_prph")
  //    {
  root_file_name += period;
  //    }
  if(use_corr) {
    if(take_2ndana_corr)
      root_file_name += "_2ndcorr";
    else
      root_file_name += "_1stcorr";
  }

  //  if(!use_corr)
  //      root_file_name += "_nocorr";

  if(use_clustered)
    root_file_name += "_uc";

  root_file_name += "_parton.root"; //warning
  //  root_file_name += "_inclusive_prph_Matt_trig_thetacut_norew.root";
  //root_file_name = "test_php_signal.root";
  //  root_file_name = "mc_php_res_prph0405e.root";
  cout << "root file name: " << root_file_name << endl;
  f = new TFile(root_file_name, "recreate", "root file with histohrams");
  const Int_t ndirs = 10;
  TDirectory* hist_directory[ndirs];//Inclusive_DIS_dir = f->mkdir("Inclusive_DIS", "inclusive DIS events");
  hist_directory[0]  = f->mkdir("Inclusive_DIS", "inclusive DIS events");
  hist_directory[1]  = f->mkdir("MC_gen_recon_corel", "correlation of generated and reconstructed values");
  hist_directory[2]  = f->mkdir("Photon_with_jet", "photon with accomp. jet events");
  hist_directory[3]  = f->mkdir("Inclusive_Photon", "inclusive photon events");
  hist_directory[4]  = f->mkdir("Cross_Sections_Histograms", "for cross section calculation");
  hist_directory[5]  = f->mkdir("Generated_Level", "generated level histograms");

  hist.Init(hist_directory);
  f->cd();
  for(Int_t i=0; i<number_of_eta_bins; i++) {
    C_scale[i] = 1.; 
  }

  for(Int_t i=0; i<number_of_eta_bins; i++)
    {
      if(use_corr)
	{
	  if(period == "989900") {
	    if(mc_corr_type == "lepto_corr")
	      C_scale[i] = lepto_factor_cells_989900[i];
	    if(mc_corr_type == "ariadne_corr")
	      C_scale[i] = ariadne_factor_cells_989900[i];
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr")
		C_scale[i] = C_scale_lepto_989900[i];
	      if(mc_corr_type == "ariadne_corr")
		C_scale[i] = C_scale_ariadne_989900[i];
	    }
	  }
	  //	  if(period == "040506e") {
	  //	    C_scale[i] = lepto_factor_cells_040506e[i];
	  //	    if(take_2ndana_corr) {
	  //	      C_scale[i] = C_scale_lepto_040506e[i];
	  //	    }
	  //	  }
	  if(period == "0405e") {
	    if(mc_corr_type == "lepto_corr")
	      C_scale[i] = lepto_factor_cells_0405e[i];
	    if(mc_corr_type == "ariadne_corr")
	      C_scale[i] = ariadne_factor_cells_0405e[i];
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr")
		C_scale[i] = C_scale_lepto_0405e[i];
	      if(mc_corr_type == "ariadne_corr")
		C_scale[i] = C_scale_ariadne_0405e[i];
	    }
	  }
	  if(period == "06e") {
	    if(mc_corr_type == "lepto_corr")
	      C_scale[i] = lepto_factor_cells_06e[i];
	    if(mc_corr_type == "ariadne_corr")
	      C_scale[i] = ariadne_factor_cells_06e[i];
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr")
		C_scale[i] = C_scale_lepto_06e[i];
	      if(mc_corr_type == "ariadne_corr")
		C_scale[i] = C_scale_ariadne_06e[i];
	    }
	  }
	  if(period == "0607p" || period == "06p" || period == "07p") {
	    if(mc_corr_type == "lepto_corr")
	      C_scale[i] = lepto_factor_cells_0607p[i];
	    if(mc_corr_type == "ariadne_corr")
	      C_scale[i] = ariadne_factor_cells_0607p[i];
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr")
		C_scale[i] = C_scale_lepto_0607p[i];
	      if(mc_corr_type == "ariadne_corr")
		C_scale[i] = C_scale_ariadne_0607p[i];
	    }
	  }
/*	  if(period == "0607p_test") {
            C_scale[i] = lepto_factor_cells_0607p[i];
            if(take_2ndana_corr) {
              C_scale[i] = C_scale_lepto_0607p_test[i];
            }
          }*/
	} else {
	C_scale[i] = C_scale_no_corr[i];
      }
    }
  cout << "======= C_scale ========" << endl;
  for(Int_t i=0; i<number_of_eta_bins; i++)
    {
      cout << C_scale[i] << endl;
    }
  
  for(Int_t i=0; i<number_of_eta_bins; i++) {
    a0_corr[i] = 0.;
    a1_corr[i] = 1.;
  };
  //a0 and a1 are for data and for detector level of monte carlo
  //q2-reweighting is only for hadron and detector level monte carlo;

  for(Int_t i=0; i<number_of_eta_bins; i++)
    {
      if(use_corr)
	{
	  if(period == "989900") {
	    if(mc_corr_type == "lepto_corr") {
	      a0_corr[i] = lepto_989900_cells[i][0];
	      a1_corr[i] = lepto_989900_cells[i][1];
	    }
	    if(mc_corr_type == "ariadne_corr") {
	      a0_corr[i] = ariadne_989900_cells[i][0];
	      a1_corr[i] = ariadne_989900_cells[i][1];
	    }
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr") {
		a0_corr[i] = a0_corr_lepto_989900[i];
		a1_corr[i] = a1_corr_lepto_989900[i];
	      }
	      if(mc_corr_type == "ariadne_corr") {
		a0_corr[i] = a0_corr_ariadne_989900[i];
		a1_corr[i] = a1_corr_ariadne_989900[i];
	      }
	    }
	  }
	  //	  if(period == "040506e") {
	  //	    a0_corr[i] = lepto_040506e_cells[i][0];
	  //	    a1_corr[i] = lepto_040506e_cells[i][1];
	  //	    if(take_2ndana_corr) {
	  //	      a0_corr[i] = a0_corr_lepto_040506e[i];
	  //	      a1_corr[i] = a1_corr_lepto_040506e[i];
	  //	    }
	  //	  }
	  if(period == "0405e") {
	    if(mc_corr_type == "lepto_corr") {
	      a0_corr[i] = lepto_0405e_cells[i][0];
	      a1_corr[i] = lepto_0405e_cells[i][1];
	    }
	    if(mc_corr_type == "ariadne_corr") {
	      a0_corr[i] = ariadne_0405e_cells[i][0];
	      a1_corr[i] = ariadne_0405e_cells[i][1];
	    }
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr") {
		a0_corr[i] = a0_corr_lepto_0405e[i];
		a1_corr[i] = a1_corr_lepto_0405e[i];
	      }
	      if(mc_corr_type == "ariadne_corr") {
		a0_corr[i] = a0_corr_ariadne_0405e[i];
		a1_corr[i] = a1_corr_ariadne_0405e[i];
	      }
	    }
	  }
	  if(period == "06e") {
	    if(mc_corr_type == "lepto_corr") {
	      a0_corr[i] = lepto_06e_cells[i][0];
	      a1_corr[i] = lepto_06e_cells[i][1];
	    }
	    if(mc_corr_type == "ariadne_corr") {
	      a0_corr[i] = ariadne_06e_cells[i][0];
	      a1_corr[i] = ariadne_06e_cells[i][1];
	    }
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr") {
		a0_corr[i] = a0_corr_lepto_06e[i];
		a1_corr[i] = a1_corr_lepto_06e[i];
	      }
	      if(mc_corr_type == "ariadne_corr") {
		a0_corr[i] = a0_corr_ariadne_06e[i];
		a1_corr[i] = a1_corr_ariadne_06e[i];
	      }
	    }
	  }
	  if(period == "0607p" || period == "06p" || period == "07p") {
	    if(mc_corr_type == "lepto_corr") {
	      a0_corr[i] = lepto_0607p_cells[i][0];
	      a1_corr[i] = lepto_0607p_cells[i][1];
	    }
	    if(mc_corr_type == "ariadne_corr") {
	      a0_corr[i] = ariadne_0607p_cells[i][0];
	      a1_corr[i] = ariadne_0607p_cells[i][1];
	    }
	    if(take_2ndana_corr) {
	      if(mc_corr_type == "lepto_corr") {
		a0_corr[i] = a0_corr_lepto_0607p[i];
		a1_corr[i] = a1_corr_lepto_0607p[i];
	      }
	      if(mc_corr_type == "ariadne_corr") {
		a0_corr[i] = a0_corr_ariadne_0607p[i];
		a1_corr[i] = a1_corr_ariadne_0607p[i];
	      }
	    }
	  }
	  /*          if(period == "0607p_test") {
            a0_corr[i] = lepto_0607p_cells[i][0];
            a1_corr[i] = lepto_0607p_cells[i][1];
            if(take_2ndana_corr) {
              a0_corr[i] = a0_corr_lepto_0607p_test[i];
              a1_corr[i] = a1_corr_lepto_0607p_test[i];
            }
          }*/
	} else {
	a0_corr[i] = a0_no_corr[i];
	a1_corr[i] = a1_no_corr[i];
      }
    }
  cout << "++++++==== a0 ---- a1 ======++++++" << endl;
  for(Int_t i=0; i<number_of_eta_bins; i++)
    {
      cout << a0_corr[i] << "   " << a1_corr[i] << endl;
    }
  /*
  h_events_per_runlumi_eltrigger_l = new TH1D("h_events_per_runlumi_eltrigger_l", "h_events_per_runlumi_eltrigger_l", 31879,30758.0, (62636.0 + 1.0));  
  h_events_per_runlumi_runs = new TH1D("h_events_per_runlumi_runs",  "h_events_per_runlumi_runs", 31879,30758.0, (62636.0 + 1.0));
  h_events_per_runlumi_largebins_lumi = new TH1D("h_events_per_runlumi_largebins_lumi", "h_events_per_runlumi_largebins_lumi",100, 30758.0, 62637.0);

  h_events_per_run = new TH1D("h_events_per_run","h_events_per_run", 31879,30758.0, (62636.0 + 1.0));
  h_dis_events_per_run = new TH1D("h_dis_events_per_run", "h_dis_events_per_run", 31879, 30758.0,	(62636.0 + 1.0));
  h_dis_ev_per_run_without_cleaning = new TH1D("h_dis_ev_per_run_without_cleaning", "h_dis_ev_per_run_without_cleaning",31879, 30758.0, (62636.0 + 1.0) );

  h_Q2_det_crosssec =  new TH1D("Q2 det crossec", "Q2 det crossec", number_Q2bins, Q2_bin);
  h_eta_det_crosssec = new TH1D("eta det crossec", "eta det crossec", number_etabins, eta_bin_crosssec);
  h_m_det_crosssec = new TH1D("m det crossec", "m det crossec", number_mbins, m_bin);
  h_xi_det_crosssec = new TH1D("xi det crossec", "xi det crossec", number_xibins, xi_bin);
  h_xbj_det_crosssec = new TH1D("xbj det crossec", "xbj det crossec", number_xbjbins, xbj_bin);
  h_etmean_det_crosssec = new TH1D("etmean det crossec", "etmean det crossec", number_etmeanbins, etmean_bin);
  tprof_Q2_det_crosssec =  new TProfile("profile Q2 det crossec", "Q2 det crossec", number_Q2bins, Q2_bin);
  tprof_eta_det_crosssec =  new TProfile("profile eta det crossec", "eta det crossec", number_etabins, eta_bin_crosssec);
  tprof_m_det_crosssec =  new TProfile("profile m det crossec", "m det crossec", number_mbins, m_bin);
  tprof_xi_det_crosssec =  new TProfile("profile xi det crossec", "xi det crossec", number_xibins, xi_bin);
  tprof_xbj_det_crosssec =  new TProfile("profile xbj det crossec", "xbj det crossec", number_xbjbins, xbj_bin);
  tprof_etmean_det_crosssec =  new TProfile("profile etmean det crossec", "etmean det crossec", number_etmeanbins, etmean_bin);

  h_Q2_had_crosssec =  new TH1D("Q2 had crossec", "Q2 had crossec", number_Q2bins, Q2_bin);
  h_eta_had_crosssec =  new TH1D("eta had crossec", "eta had crossec", number_etabins, eta_bin_crosssec);
  h_m_had_crosssec =  new TH1D("m had crossec", "m had crossec", number_mbins, m_bin);
  h_xi_had_crosssec =  new TH1D("xi had crossec", "xi had crossec", number_xibins, xi_bin);
  h_xbj_had_crosssec =  new TH1D("xbj had crossec", "xbj had crossec", number_xbjbins, xbj_bin);
  h_etmean_had_crosssec =  new TH1D("etmean had crossec", "etmean had crossec", number_etmeanbins, etmean_bin);

  tprof_Q2_had_crosssec =  new TProfile("profile Q2 had crossec", "Q2 had crossec", number_Q2bins, Q2_bin);
  tprof_eta_had_crosssec =  new TProfile("profile eta had crossec", "eta had crossec", number_etabins, eta_bin_crosssec);
  tprof_m_had_crosssec =  new TProfile("profile m had crossec", "m had crossec", number_mbins, m_bin);
  tprof_xi_had_crosssec =  new TProfile("profile xi had crossec", "xi had crossec", number_xibins, xi_bin);
  tprof_xbj_had_crosssec =  new TProfile("profile xbj had crossec", "xbj had crossec", number_xbjbins, xbj_bin);
  tprof_etmean_had_crosssec =  new TProfile("profile etmean had crossec", "etmean had crossec", number_etmeanbins, etmean_bin);
  h_Q2_part_crosssec =  new TH1D("Q2 part crossec", "Q2 part crossec", number_Q2bins, Q2_bin);
  h_eta_part_crosssec =  new TH1D("eta part crossec", "eta part crossec", number_etabins, eta_bin_crosssec);
  h_m_part_crosssec =  new TH1D("m part crossec", "m part crossec", number_mbins, m_bin);
  h_xi_part_crosssec =  new TH1D("xi part crossec", "xi part crossec", number_xibins, xi_bin);
  h_xbj_part_crosssec =  new TH1D("xbj part crossec", "xbj part crossec", number_xbjbins, xbj_bin);
  h_etmean_part_crosssec =  new TH1D("etmean part crossec", "etmean part crossec", number_etmeanbins, etmean_bin);

  tprof_Q2_part_crosssec =  new TProfile("profile Q2 part crossec", "Q2 part crossec", number_Q2bins, Q2_bin);
  tprof_eta_part_crosssec =  new TProfile("profile eta part crossec", "eta part crossec", number_etabins, eta_bin_crosssec);
  tprof_m_part_crosssec =  new TProfile("profile m part crossec", "m part crossec", number_mbins, m_bin);
  tprof_xi_part_crosssec =  new TProfile("profile xi part crossec", "xi part crossec", number_xibins, xi_bin);
  tprof_xbj_part_crosssec =  new TProfile("profile xbj part crossec", "xbj part crossec", number_xbjbins, xbj_bin);
  tprof_etmean_part_crosssec =  new TProfile("profile etmean part crossec", "etmean part crossec", number_etmeanbins, etmean_bin);

  h_Q2_hd_crosssec =  new TH1D("Q2 hd crossec", "Q2 hd crossec", number_Q2bins, Q2_bin);
  h_eta_hd_crosssec =  new TH1D("eta hd crossec", "eta hd crossec", number_etabins, eta_bin_crosssec);
  h_m_hd_crosssec =  new TH1D("m hd crossec", "m hd crossec", number_mbins, m_bin);
  h_xi_hd_crosssec =  new TH1D("xi hd crossec", "xi hd crossec", number_xibins, xi_bin);
  h_xbj_hd_crosssec =  new TH1D("xbj hd crossec", "xbj hd crossec", number_xbjbins, xbj_bin);
  h_etmean_hd_crosssec =  new TH1D("etmean hd crossec", "etmean hd crossec", number_etmeanbins, etmean_bin);

  h_xi_diff_q2_bins_1 = new TH1D("xi diff q2 bin 1", "xi diff q2 bin 1", xidiffq2bins2_1, xi_diff_q2_bins2_1);
  h_xi_diff_q2_bins_2 = new TH1D("xi diff q2 bin 2", "xi diff q2 bin 2", xidiffq2bins2_2, xi_diff_q2_bins2_2);
  h_xi_diff_q2_bins_3 = new TH1D("xi diff q2 bin 3", "xi diff q2 bin 3", xidiffq2bins2_3, xi_diff_q2_bins2_3);
  h_xi_diff_q2_bins_4 = new TH1D("xi diff q2 bin 4", "xi diff q2 bin 4", xidiffq2bins2_4, xi_diff_q2_bins2_4);
  h_xi_diff_q2_bins_5 = new TH1D("xi diff q2 bin 5", "xi diff q2 bin 5", xidiffq2bins2_5, xi_diff_q2_bins2_5);
  h_xi_diff_q2_bins_6 = new TH1D("xi diff q2 bin 6", "xi diff q2 bin 6", xidiffq2bins2_6, xi_diff_q2_bins2_6);
  tprof_xi_diff_q2_bins_1 =  new TProfile("prof_xi_diff_q2_bins_1", "xi_diff_q2_bins_1", xidiffq2bins2_1, xi_diff_q2_bins2_1);
  tprof_xi_diff_q2_bins_2 =  new TProfile("prof_xi_diff_q2_bins_2", "xi_diff_q2_bins_2", xidiffq2bins2_2, xi_diff_q2_bins2_2);
  tprof_xi_diff_q2_bins_3 =  new TProfile("prof_xi_diff_q2_bins_3", "xi_diff_q2_bins_3", xidiffq2bins2_3, xi_diff_q2_bins2_3);
  tprof_xi_diff_q2_bins_4 =  new TProfile("prof_xi_diff_q2_bins_4", "xi_diff_q2_bins_4", xidiffq2bins2_4, xi_diff_q2_bins2_4);
  tprof_xi_diff_q2_bins_5 =  new TProfile("prof_xi_diff_q2_bins_5", "xi_diff_q2_bins_5", xidiffq2bins2_5, xi_diff_q2_bins2_5);
  tprof_xi_diff_q2_bins_6 =  new TProfile("prof_xi_diff_q2_bins_6", "xi_diff_q2_bins_6", xidiffq2bins2_6, xi_diff_q2_bins2_6);

  h_had_xi_diff_q2_bins_1 = new TH1D("had xi diff q2 bin 1", "had xi diff q2 bin 1", xidiffq2bins2_1, xi_diff_q2_bins2_1);
  h_had_xi_diff_q2_bins_2 = new TH1D("had xi diff q2 bin 2", "had xi diff q2 bin 2", xidiffq2bins2_2, xi_diff_q2_bins2_2);
  h_had_xi_diff_q2_bins_3 = new TH1D("had xi diff q2 bin 3", "had xi diff q2 bin 3", xidiffq2bins2_3, xi_diff_q2_bins2_3);
  h_had_xi_diff_q2_bins_4 = new TH1D("had xi diff q2 bin 4", "had xi diff q2 bin 4", xidiffq2bins2_4, xi_diff_q2_bins2_4);
  h_had_xi_diff_q2_bins_5 = new TH1D("had xi diff q2 bin 5", "had xi diff q2 bin 5", xidiffq2bins2_5, xi_diff_q2_bins2_5);
  h_had_xi_diff_q2_bins_6 = new TH1D("had xi diff q2 bin 6", "had xi diff q2 bin 6", xidiffq2bins2_6, xi_diff_q2_bins2_6);

  h_hd_xi_diff_q2_bins_1 = new TH1D("hd xi diff q2 bin 1", "hd xi diff q2 bin 1", xidiffq2bins2_1, xi_diff_q2_bins2_1);
  h_hd_xi_diff_q2_bins_2 = new TH1D("hd xi diff q2 bin 2", "hd xi diff q2 bin 2", xidiffq2bins2_2, xi_diff_q2_bins2_2);
  h_hd_xi_diff_q2_bins_3 = new TH1D("hd xi diff q2 bin 3", "hd xi diff q2 bin 3", xidiffq2bins2_3, xi_diff_q2_bins2_3);
  h_hd_xi_diff_q2_bins_4 = new TH1D("hd xi diff q2 bin 4", "hd xi diff q2 bin 4", xidiffq2bins2_4, xi_diff_q2_bins2_4);
  h_hd_xi_diff_q2_bins_5 = new TH1D("hd xi diff q2 bin 5", "hd xi diff q2 bin 5", xidiffq2bins2_5, xi_diff_q2_bins2_5);
  h_hd_xi_diff_q2_bins_6 = new TH1D("hd xi diff q2 bin 6", "hd xi diff q2 bin 6", xidiffq2bins2_6, xi_diff_q2_bins2_6);

  h_etmean_diff_q2_bins_1 = new TH1D("etmean diff q2 bin 1", "etmean diff q2 bin 1", etmeandiffq2bins2_1, etmean_diff_q2_bins2_1);
  h_etmean_diff_q2_bins_2 = new TH1D("etmean diff q2 bin 2", "etmean diff q2 bin 2", etmeandiffq2bins2_2, etmean_diff_q2_bins2_2);
  h_etmean_diff_q2_bins_3 = new TH1D("etmean diff q2 bin 3", "etmean diff q2 bin 3", etmeandiffq2bins2_3, etmean_diff_q2_bins2_3);
  h_etmean_diff_q2_bins_4 = new TH1D("etmean diff q2 bin 4", "etmean diff q2 bin 4", etmeandiffq2bins2_4, etmean_diff_q2_bins2_4);
  h_etmean_diff_q2_bins_5 = new TH1D("etmean diff q2 bin 5", "etmean diff q2 bin 5", etmeandiffq2bins2_5, etmean_diff_q2_bins2_5);
  h_etmean_diff_q2_bins_6 = new TH1D("etmean diff q2 bin 6", "etmean diff q2 bin 6", etmeandiffq2bins2_6, etmean_diff_q2_bins2_6);
  tprof_etmean_diff_q2_bins_1 =  new TProfile("prof_etmean_diff_q2_bins_1", "etmean_diff_q2_bins_1", etmeandiffq2bins2_1, etmean_diff_q2_bins2_1);
  tprof_etmean_diff_q2_bins_2 =  new TProfile("prof_etmean_diff_q2_bins_2", "etmean_diff_q2_bins_2", etmeandiffq2bins2_2, etmean_diff_q2_bins2_2);
  tprof_etmean_diff_q2_bins_3 =  new TProfile("prof_etmean_diff_q2_bins_3", "etmean_diff_q2_bins_3", etmeandiffq2bins2_3, etmean_diff_q2_bins2_3);
  tprof_etmean_diff_q2_bins_4 =  new TProfile("prof_etmean_diff_q2_bins_4", "etmean_diff_q2_bins_4", etmeandiffq2bins2_4, etmean_diff_q2_bins2_4);
  tprof_etmean_diff_q2_bins_5 =  new TProfile("prof_etmean_diff_q2_bins_5", "etmean_diff_q2_bins_5", etmeandiffq2bins2_5, etmean_diff_q2_bins2_5);
  tprof_etmean_diff_q2_bins_6 =  new TProfile("prof_etmean_diff_q2_bins_6", "etmean_diff_q2_bins_6", etmeandiffq2bins2_6, etmean_diff_q2_bins2_6);

  h_had_etmean_diff_q2_bins_1 = new TH1D("had etmean diff q2 bin 1", "had etmean diff q2 bin 1", etmeandiffq2bins2_1, etmean_diff_q2_bins2_1);
  h_had_etmean_diff_q2_bins_2 = new TH1D("had etmean diff q2 bin 2", "had etmean diff q2 bin 2", etmeandiffq2bins2_2,      etmean_diff_q2_bins2_2);
  h_had_etmean_diff_q2_bins_3 = new TH1D("had etmean diff q2 bin 3", "had etmean diff q2 bin 3", etmeandiffq2bins2_3,      etmean_diff_q2_bins2_3);
  h_had_etmean_diff_q2_bins_4 = new TH1D("had etmean diff q2 bin 4", "had etmean diff q2 bin 4", etmeandiffq2bins2_4,      etmean_diff_q2_bins2_4);
  h_had_etmean_diff_q2_bins_5 = new TH1D("had etmean diff q2 bin 5", "had etmean diff q2 bin 5", etmeandiffq2bins2_5, etmean_diff_q2_bins2_5);
  h_had_etmean_diff_q2_bins_6 = new TH1D("had etmean diff q2 bin 6", "had etmean diff q2 bin 6", etmeandiffq2bins2_6, etmean_diff_q2_bins2_6);
  h_hd_etmean_diff_q2_bins_1 = new TH1D("hd etmean diff q2 bin 1", "hd etmean diff q2 bin 1", etmeandiffq2bins2_1, etmean_diff_q2_bins2_1);
  h_hd_etmean_diff_q2_bins_2 = new TH1D("hd etmean diff q2 bin 2", "hd etmean diff q2 bin 2", etmeandiffq2bins2_2,         etmean_diff_q2_bins2_2);
  h_hd_etmean_diff_q2_bins_3 = new TH1D("hd etmean diff q2 bin 3", "hd etmean diff q2 bin 3", etmeandiffq2bins2_3,         etmean_diff_q2_bins2_3);
  h_hd_etmean_diff_q2_bins_4 = new TH1D("hd etmean diff q2 bin 4", "hd etmean diff q2 bin 4", etmeandiffq2bins2_4,         etmean_diff_q2_bins2_4);
  h_hd_etmean_diff_q2_bins_5 = new TH1D("hd etmean diff q2 bin 5", "hd etmean diff q2 bin 5", etmeandiffq2bins2_5, etmean_diff_q2_bins2_5);
  h_hd_etmean_diff_q2_bins_6 = new TH1D("hd etmean diff q2 bin 6", "hd etmean diff q2 bin 6", etmeandiffq2bins2_6, etmean_diff_q2_bins2_6);

  //trijets cross_sections
  h_trijet_Q2_det_crosssec = new TH1D("trijet Q2 det crossec", "trijet Q2 det crossec", number_Q2bins_trijet, Q2_bin_trijet);
  h_trijet_eta_det_crosssec = new TH1D("trijet eta det crossec", "trijet eta det crossec", number_etabins, eta_bin_crosssec);
  h_trijet_m_det_crosssec = new TH1D("trijet m det crossec", "trijet m det crossec", number_mbins, m_bin);
  h_trijet_xi_det_crosssec = new TH1D("trijet xi det crossec", "trijet xi det crossec", number_xibins, xi_bin);
  h_trijet_xbj_det_crosssec = new TH1D("trijet xbj det crossec", "trijet xbj det crossec", number_xbjbins, xbj_bin);
  h_trijet_etmean_det_crosssec = new TH1D("trijet etmean det crossec", "trijet etmean det crossec", number_etmeanbins, etmean_bin);

  tprof_trijet_Q2_det_crosssec = new TProfile("profile trijet Q2 det crossec", "trijet Q2 det crossec", number_Q2bins_trijet, Q2_bin_trijet);
  tprof_trijet_eta_det_crosssec = new TProfile("profile trijet eta det crossec", "trijet eta det crossec", number_etabins, eta_bin_crosssec);
  tprof_trijet_m_det_crosssec = new TProfile("profile trijet m det crossec", "trijet m det crossec", number_mbins, m_bin);
  tprof_trijet_xi_det_crosssec = new TProfile("profile trijet xi det crossec", "trijet xi det crossec", number_xibins, xi_bin);
  tprof_trijet_xbj_det_crosssec = new TProfile("profile trijet xbj det crossec", "trijet xbj det crossec", number_xbjbins, xbj_bin);
  tprof_trijet_etmean_det_crosssec = new TProfile("profile trijet etmean det crossec", "trijet etmean det crossec", number_etmeanbins, etmean_bin);

  h_trijet_Q2_had_crosssec = new TH1D("trijet Q2 had crossec", "trijet Q2 had crossec", number_Q2bins_trijet, Q2_bin_trijet);
  h_trijet_eta_had_crosssec = new TH1D("trijet eta had crossec", "trijet eta had crossec", number_etabins, eta_bin_crosssec);
  h_trijet_m_had_crosssec = new TH1D("trijet m had crossec", "trijet m had crossec", number_mbins, m_bin);
  h_trijet_xi_had_crosssec = new TH1D("trijet xi had crossec", "trijet xi had crossec", number_xibins, xi_bin);
  h_trijet_xbj_had_crosssec = new TH1D("trijet xbj had crossec", "trijet xbj had crossec", number_xbjbins, xbj_bin);
  h_trijet_etmean_had_crosssec = new TH1D("trijet etmean had crossec", "trijet etmean had crossec", number_etmeanbins, etmean_bin);

  tprof_trijet_Q2_had_crosssec = new TProfile("profile trijet Q2 had crossec", "trijet Q2 had crossec", number_Q2bins_trijet, Q2_bin_trijet);
  tprof_trijet_eta_had_crosssec = new TProfile("profile trijet eta had crossec", "trijet eta had crossec", number_etabins, eta_bin_crosssec);
  tprof_trijet_m_had_crosssec = new TProfile("profile trijet m had crossec", "trijet m had crossec", number_mbins, m_bin);
  tprof_trijet_xi_had_crosssec = new TProfile("profile trijet xi had crossec", "trijet xi had crossec", number_xibins, xi_bin);
  tprof_trijet_xbj_had_crosssec = new TProfile("profile trijet xbj had crossec", "trijet xbj had crossec", number_xbjbins, xbj_bin);
  tprof_trijet_etmean_had_crosssec = new TProfile("profile trijet etmean had crossec", "trijet etmean had crossec", number_etmeanbins, etmean_bin);

  h_trijet_Q2_part_crosssec = new TH1D("trijet Q2 part crossec", "trijet Q2 part crossec", number_Q2bins_trijet, Q2_bin_trijet);
  h_trijet_eta_part_crosssec = new TH1D("trijet eta part crossec", "trijet eta part crossec", number_etabins, eta_bin_crosssec);
  h_trijet_m_part_crosssec = new TH1D("trijet m part crossec", "trijet m part crossec", number_mbins, m_bin);
  h_trijet_xi_part_crosssec = new TH1D("trijet xi part crossec", "trijet xi part crossec", number_xibins, xi_bin);
  h_trijet_xbj_part_crosssec = new TH1D("trijet xbj part crossec", "trijet xbj part crossec", number_xbjbins, xbj_bin);
  h_trijet_etmean_part_crosssec = new TH1D("trijet etmean part crossec", "trijet etmean part crossec", number_etmeanbins, etmean_bin);

  tprof_trijet_Q2_part_crosssec = new TProfile("profile trijet Q2 part crossec", "trijet Q2 part crossec", number_Q2bins_trijet, Q2_bin_trijet);
  tprof_trijet_eta_part_crosssec = new TProfile("profile trijet eta part crossec", "trijet eta part crossec", number_etabins, eta_bin_crosssec);
  tprof_trijet_m_part_crosssec = new TProfile("profile trijet m part crossec", "trijet m part crossec", number_mbins, m_bin);
  tprof_trijet_xi_part_crosssec = new TProfile("profile trijet xi part crossec", "trijet xi part crossec", number_xibins, xi_bin);
  tprof_trijet_xbj_part_crosssec = new TProfile("profile trijet xbj part crossec", "trijet xbj part crossec", number_xbjbins, xbj_bin);
  tprof_trijet_etmean_part_crosssec = new TProfile("profile trijet etmean part crossec", "trijet etmean part crossec", number_etmeanbins, etmean_bin);

  h_trijet_Q2_hd_crosssec = new TH1D("trijet Q2 hd crossec", "trijet Q2 hd crossec", number_Q2bins_trijet, Q2_bin_trijet);
  h_trijet_eta_hd_crosssec = new TH1D("trijet eta hd crossec", "trijet eta hd crossec", number_etabins, eta_bin_crosssec);
  h_trijet_m_hd_crosssec = new TH1D("trijet m hd crossec", "trijet m hd crossec", number_mbins, m_bin);
  h_trijet_xi_hd_crosssec = new TH1D("trijet xi hd crossec", "trijet xi hd crossec", number_xibins, xi_bin);
  h_trijet_xbj_hd_crosssec = new TH1D("trijet xbj hd crossec", "trijet xbj hd crossec", number_xbjbins, xbj_bin);
  h_trijet_etmean_hd_crosssec = new TH1D("trijet etmean hd crossec", "trijet etmean hd crossec", number_etmeanbins, etmean_bin);

  h_pJetsEtBreit = new TH1D("pJetsEtBreit", "pJetsEtBreit", 160, 0., 80.);
  h_hJetsEtBreit = new TH1D("hJetsEtBreit", "hJetsEtBreit", 160, 0., 80.);
  h_JetsEtBreit = new TH1D("JetsEtBreit", "JetsEtBreit", 20, 0., 100.);
  h_Zvtx = new TH1D("Zvtx", "Zvtx", 100, -50., 50.);
  h_E_da = new TH1D("E_da", "E_da", 100, 0., 40.);
  h_E_el = new TH1D("E_el", "E_el", 100, 0., 40.);
  h_yjb_da = new TH1D("yjb_da", "yjb_da", 100, 0.1, 0.7);
  h_Q2 = new TH1D("h_Q2", "h_Q2", 100, 125., 20000.);
  h_pt_sqrtet = new TH1D("h_pt_sqrtet", "h_pt_sqrtet", 100, 0., 2.5);
  h_hac_cell = new TH1D("h_hac_cell", "h_hac_cell", 100, 0., 80.);
  h_emc_cell = new TH1D("h_emc_cell", "h_emc_cell", 100, 0., 80.);

  //dijet observables: detector level
  h_q2_da = new TH1D("q2_da", "q2_da", 35, 125., 20000.);
  h_M_jj = new TH1D("M_jj", "M_jj", 20, 20., 85.);
  h_log_ksi = new TH1D("log10ksi", "log10ksi", 20, -2.2, 0.);
  h_x_Bj = new TH1D("x_Bj", "x_Bj", 20, 0., 0.1);
  h_jet_Et_mean_Breit = new TH1D("jet_Et_mean_Breit", "jet_Et_mean_Breit", 22, 8., 45.);
  h_Eta_str = new TH1D("Eta_str", "Eta_str", 20, 0., 2.0);

  h_jet1_TBreit = new TH1D("jet1_TBreit", "jet1_TBreit", 160, 0., 80.);
  h_jet2_TBreit = new TH1D("jet2_TBreit", "jet2_TBreit", 160, 0., 80.);

  h_jet1_etaBreit = new TH1D("jet1_etaBreit", "jet1_etaBreit", 100, -2.25, 2.75);
  h_jet2_etaBreit = new TH1D("jet2_etaBreit", "jet2_etaBreit", 100, -2.25, 2.75);

  //dijet observables: hadron level
  h_had_q2_da = new TH1D("had_q2_da", "had_q2_da", 35, 125., 20000.);
  h_had_M_jj = new TH1D("had_M_jj", "had_M_jj", 20, 20., 85.);
  h_had_log_ksi = new TH1D("had_log10ksi", "had_log10ksi", 20, -2.2, 0.);
  h_had_x_Bj = new TH1D("had_x_Bj", "had_x_Bj", 20, 0., 0.1);
  h_had_jet_Et_mean_Breit = new TH1D("had_jet_Et_mean_Breit", "had_jet_Et_mean_Breit", 22, 8., 45.);
  h_had_Eta_str = new TH1D("had_Eta_str", "had_Eta_str", 20, 0., 2.0);

  h_had_jet1_etaBreit = new TH1D("had_jet1_etaBreit", "had_jet1_etaBreit", 100, -2.25, 2.75);
  h_had_jet2_etaBreit = new TH1D("had_jet2_etaBreit", "had_jet2_etaBreit", 100, -2.25, 2.75);
  h_had_jet1_TBreit = new TH1D("had_jet1_TBreit", "had_jet1_TBreit", 160, 0., 80.);
  h_had_jet2_TBreit = new TH1D("had_jet2_TBreit", "had_jet2_TBreit", 160, 0., 80.);
  //dijet observables: parton level
  h_part_q2_da = new TH1D("part_q2_da", "part_q2_da", 35, 125., 20000.);
  h_part_M_jj = new TH1D("part_M_jj", "part_M_jj", 20, 20., 85.);
  h_part_log_ksi = new TH1D("part_log10ksi", "part_log10ksi", 20, -2.2, 0.);
  h_part_x_Bj = new TH1D("part_x_Bj", "part_x_Bj", 20, 0., 0.1);
  h_part_jet_Et_mean_Breit = new TH1D("part_jet_Et_mean_Breit", "part_jet_Et_mean_Breit", 22, 8., 45.);
  h_part_Eta_str = new TH1D("part_Eta_str", "part_Eta_str", 20, 0., 2.0);

  h_part_jet1_etaBreit = new TH1D("part_jet1_etaBreit", "part_jet1_etaBreit", 100, -2.25, 2.75);
  h_part_jet2_etaBreit = new TH1D("part_jet2_etaBreit", "part_jet2_etaBreit", 100, -2.25, 2.75);
  h_part_jet1_TBreit = new TH1D("part_jet1_TBreit", "part_jet1_TBreit", 160, 0., 80.);
  h_part_jet2_TBreit = new TH1D("part_jet2_TBreit", "part_jet2_TBreit", 160, 0., 80.);

  h_Eel_Eda = new TH1D("Eel/Eda", "Eel/Eda", 40, 0., 2.);
  h_theta = new TH1D("Theta", "Theta", 40, 1.5, 3.);


  // forward & backward jets
  h_2jets_control_eta_b = new TH1D("h_2jets_control_eta_b", "h_2jets_control_eta_b", 25, -1.0, 2.5);
  h_2jets_control_eta_f = new TH1D("h_2jets_control_eta_f", "h_2jets_control_eta_f", 25, -1.0, 2.5);

  h_2jets_control_eta_breit_b = new TH1D("h_2jets_control_eta_breit_b", "h_2jets_control_eta_breit_b", 25, -2.0, 3.5);
  h_2jets_control_eta_breit_f = new TH1D("h_2jets_control_eta_breit_f", "h_2jets_control_eta_breit_f", 25, -2.0, 3.5);
  h_2jets_control_eta = new TH1D("h_2jets_control_eta", "h_2jets_control_eta", 25, -1.0, 2.5);
  h_2jets_control_eta_breit = new TH1D("h_2jets_control_eta_breit", "h_2jets_control_eta_breit", 25, -2.0, 3.5);

  //dijets per run
  h_dijets_per_run = new TH1D("h_dijets_per_run","h_dijets_per_run", 31879,30758.0, (62636.0 + 1.0));

  //    TH1D h_xvtx("xvtx", "xvtx", 40, -5., 5.);

  //electron position in RCAL sipos
  h_elec_position = new TH2D("h_elec_position", "elec_position", 1000, -250., 250., 1000, -250., 250.);
  h_elec_position->Sumw2();
  h_elec_position_x = new TH1D("h_elec_position_x", "elec_position_x", 1000, -250., 250.);
  h_elec_position_y = new TH1D("h_elec_position_y", "elec_position_y", 1000, -250., 250);
  //q2 rekonstructed and real
  h_q2_comp = new TH2D("h_q2_comp", "h_q2_comp", 1000, 0., 20000., 1000, 0., 20000.);
  h_q2_deviation = new TH1D("h_q2_deviation", "h_q2_deviation", 100, -1., 1.);
  h_elec_position_x->Sumw2();
  h_elec_position_y->Sumw2();

  //Fltw bits
  h_fltw = new TH1D("h_fltw", "h_fltw", 64, 0., 64.); h_fltw->Sumw2();
  h_fltw_p = new TH1D("h_fltw_p", "h_fltw_p", 64, 0., 64.); h_fltw_p->Sumw2();
  h_fltw_dis = new TH1D("h_fltw_dis", "h_fltw_dis", 64, 0., 64.); h_fltw_dis->Sumw2();
  h_fltw_dis_p = new TH1D("h_fltw_dis_p", "h_fltw_dis_p", 64, 0., 64.); h_fltw_dis_p->Sumw2();

  //jets phase space
  h_et1et2 = new TH2D("h_et1et2", "h_et1_et2", 160, 0., 80., 160, 0., 80.); h_et1et2->Sumw2();
  h_et1et2_mcut = new TH2D("h_et1et2_mcut", "h_et1et2", 160, 0., 80., 160, 0., 80.); h_et1et2_mcut->Sumw2();

  h_et1et2_breit = new TH2D("h_et1et2_breit", "h_et1et2_breit", 160, 0., 80., 160, 0., 80.); h_et1et2_breit->Sumw2();
  h_et1et2_mcut_breit = new TH2D("h_et1et2_mcut_breit", "h_et1et2_mcut_breit", 160, 0., 80., 160, 0., 80.); h_et1et2_mcut_breit->Sumw2();

  //track veto correction
  h_yda_veto_bit30 = new TH1D("h_yda_veto_bit30", "yda veto bit30", 15, 0.2, 0.6);
  h_yda_veto_bit30->Sumw2();
  h_yda_bit30 = new TH1D("h_yda_bit30", "yda bit30", 15, 0.2, 0.6);
  h_yda_bit30->Sumw2();
  */

  if(Data)
    pcele_mode = 2;
  else
    pcele_mode = 20;
  bpres_conerad = 1.0;

  //  counters
  take_det_dijet = 0;
  cutC_1 = 0;
  cutC_2 = 0;
  cutC_3 = 0;
  cutPS_1 = 0;
  cutPS_2 = 0;
  cutEC_1 = 0;
  cutEC_2 = 0;
  cutEC_3 = 0;
  cutEC_4 = 0;
  cutEC_5 = 0;
  cutEC_6_1 = 0;
  cutEC_6_2 = 0;
  cutEC_6_3 = 0;
  cutEC_6_4 = 0;
  cutElCand = 0;
  cutCTDacc = 0;
  cutJPS_1 = 0;
  cutJPS_23 = 0;
  cutJC_1 = 0;
  cutJC_2 = 0;
  cutJC_3 = 0;
  totalSelectedJetsCount = 0;
  localSelectedJetsCount = 0;
  goodCellsCount = 0;
  selectedEventsCount = 0;
  totalNumberOfEntries = 0;
  processedEventsCount = 0;
  cut_DIS03 = 0;
  cut_Sincand = 0;
  cut_Siecorr = 0;
  cut_q2_bos = 0;
  cut_y_bos = 0;
  cut_part_Et_jet = 0;
  cut_part_eta = 0;
  cut_hadr_Et_jet = 0;
  cut_hadr_eta = 0;
  cut_assymetric = 0;
  cut_track_momentum = 0;
  cut_electron_track_dca = 0;
  cut_electron_distance_module_edge = 0;
  cut_Mjj = 0;
  partonDijetEvents = 0;
  hadronDijetEvents = 0;
  detector_dijets = 0;
  cut_poltake = 0;
  cut_QEDC = 0;
  cut_yx2 = 0;
  cut_p_Mjj = 0;
  cut_h_Mjj = 0;
  cut_ngoodtrks = 0;
  cut_gammahad = 0;
  cut_RCAL = 0;
  cut_triggers = 0;
  count_nat = 0;
  v_uncorr_prompt_photon = new TLorentzVector();
  v_corr_prompt_photon = new TLorentzVector();
  v_prompt_photon_jet = new TLorentzVector();
  v_true_prompt_photon = new TLorentzVector();
  v_true_jet_cont_prompt_photon = new TLorentzVector();
  v_accomp_corr_jet = new TLorentzVector();
  v_accomp_uncorr_jet = new TLorentzVector();
  v_true_acc_jet = new TLorentzVector();
  v_true_parton_acc_jet = new TLorentzVector();
  v_true_scattered_electron = new TLorentzVector();
  glob_fmax = -999.;
  glob_deltaz = -999.;
  hadron_level_jet_cont_photon_index = -1;
  hadron_level_ephoton_over_ejet = -1;
  if(!Data) getBpreReadout();
  if(Data) getRunsBPRE();

  ET_JET_CUT = 2.5;//warning
  inclusive_prph_analysis = kFALSE; //warning;
  fCalib = TFile::Open("calib/photcalib.v3.5.root");
  had_et = -999.;
  had_eta = -999.;
  had_et_jet = -999.;
  had_eta_jet = -999.;
  had_x = -999.;
  had_Q2 = -999.;

  part_et = -999.;
  part_eta = -999.;
  part_et_jet = -999.;
  part_eta_jet = -999.;
  part_x = -999.;
  part_Q2 = -999.;
}

void selector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  //   TString option = GetOption();

}
