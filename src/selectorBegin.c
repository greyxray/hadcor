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

  root_file_name += "_parton";
  root_file_name +=  "_" + TString::Itoa(filenum, 10);
  root_file_name += ".root"; //warning
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
