//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 19 10:10:56 2010 by ROOT version 5.20/00
// from TChain orange/
//////////////////////////////////////////////////////////

//#ifndef selector_h
//#define selector_h
#ifndef selector_cxx
#define selector_cxx
#include <iostream>
#include <string>
using namespace std;

//#include <TROOT.h>
#include <TChain.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TString.h>
#include <TFile.h>
#include <TDatime.h>
#include <TRandom3.h>
#include <TEventList.h>
#include "constants.h"
#include "hist.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
using KtJet::KtLorentzVector;
using KtJet::KtEvent;

class selector {
 public :

  static const Float_t m_Lepton_Energy = 27.52;
  static const Float_t m_Proton_Energy = 920.;
  static const  Double_t E_e = 27.52;
  static const  Double_t E_p = 920.;
  static const Int_t maxNjetsInEvent = 50;
  static const Int_t maxNofJets = 250;

  bool en_mom_conservation;//per event flag
  bool check_en_mom_conservation;//per process at all
  bool check_en_mom_conservation_on_parton_level;//per process at all
  static const Double_t px_cons = - 0.35;
  static const Double_t py_cons = - 0.22;
  static const Double_t pz_cons = 892.5;
  static const Double_t E_cons = 947.5;

    bool nodebugmode;
   Double_t delta_phi(Double_t phi1, Double_t phi2);
  Double_t jet_en_corr(Double_t eta, Double_t et, TString period, TString mc_type);
  Double_t jet_en_corr(Double_t eta, Double_t et, Double_t* C_scale, Double_t* a0_corr, Double_t* a1_corr, Double_t* a2_corr);
  Double_t stretch_calib(Double_t value,       // Value of the distribution
			 string method,       // "e5" for ELEC5 or "zu" for Zufo
			 string distribution, // Fmax or deltaZ
			 Double_t eta,
			 Double_t et);
  Double_t BPRES_mips(TVector3 v_photon, Double_t radius);
  void  getBpreReadout(); // A.Volynets
  Double_t doBpreReadout(Double_t mips); // A.Volynets
  void getRunsBPRE(); // A.Volynets
  Double_t Land[100];          // BPRE smearing constants
  map<int,int>   BadBpre;
  void fill_trigger_bits(TH1D* h[], Int_t trigger[]);
  void fill_trigger_bits_general(TH1D* h[], Int_t trigger[], Int_t n);
  Bool_t SelectHadronLevel(Bool_t take_det_event);
  Bool_t SelectPartonLevel(Bool_t take_det_event, Bool_t take_had_event);
  Bool_t SelectPrPhKtJetB(Int_t i_true_rad_photon,  Int_t electron_in_zufo);
  Double_t q2_reweighting(Float_t Mc_q2, TString mc_type);
  Double_t etaph_reweighting(Double_t eta, TString mc_type);
    Double_t eph_reweighting(Double_t e, TString mc_type);
  Double_t q2_reweighting_separately(Float_t Mc_q2, TString mc_type, TString period);
  /* THIS ALL WAS NOT PRESENT HERE
    Double_t y_reweighting(Float_t Mc_y, TString mc_type);
    Double_t y_reweighting_binbybin(Float_t Mc_y, TString mc_type);
    Double_t eph_reweighting(Double_t e, TString mc_type);
    Double_t etph_bg_reweighting(Double_t et);
    Double_t eph_bg_reweighting(Double_t e);
    Double_t etaph_reweighting(Double_t eta, TString mc_type);
    Double_t etph_reweighting(Double_t et, TString mc_type);
  */
  //THIS IS SOMETHING NEW - BLOCK
    vector<KtJet::KtLorentzVector> hadron_level_jets;
    Int_t           hadron_level_jet_cont_photon_index;
    Double_t           hadron_level_ephoton_over_ejet;
    Double_t        had_et;
    Double_t        had_eta;
    Double_t        had_et_jet;
    Double_t        had_eta_jet;
    Double_t        had_x;
    Double_t        had_Q2;
    Double_t         had_xgamma;
    Double_t        had_xp;
    Double_t        had_dphi;
    Double_t        had_deta;
    Double_t        had_dphi_e_ph;
    Double_t        had_deta_e_ph;
    
    Double_t        part_et;
    Double_t        part_eta;
    Double_t        part_et_jet;
    Double_t        part_eta_jet;
    Double_t        part_x;
    Double_t        part_Q2;
    Double_t        part_xgamma;
    Double_t        part_xp;
    Double_t        part_dphi;
    Double_t        part_deta;
    Double_t        part_dphi_e_ph;
    Double_t        part_deta_e_ph;
  //--------------------------------
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TEventList     *event_list; 
  Hist            hist;
  Bool_t          Data;
  Bool_t          use_corr;
  Bool_t          take_2ndana_corr;
  Bool_t          use_clustered;
  Bool_t          check_cuts;
  Bool_t          use_ktjetb_prph;
  Bool_t          here_is_prph;
  Bool_t          inclusive_prph_analysis;
  Double_t        jet_energy_frac;
  Double_t        cell_energy_frac;
  Double_t        photon_corr_et;
  Double_t        max_et_candidate;
  Double_t        max_et_candidate_uncorr;
  Double_t        max_et_accomp_jet_had_level;// THIS IS SOMETHING NEW
  Double_t        max_et_accomp_jet_had_level_eta;// THIS IS SOMETHING NEW
  Int_t           index_phot_vector;// NOT PRESENTED HERE but I added
  Double_t        glob_fmax;
  Double_t        glob_deltaz;
  Int_t           det_cross_sec_bin_et;
  Int_t           det_cross_sec_bin_eta;
  Int_t           det_cross_sec_bin_q2;
  Int_t           det_cross_sec_bin_x;
  Int_t           det_cross_sec_bin_et_jet;
  Int_t           det_cross_sec_bin_eta_jet;
  Int_t           det_cross_sec_bin_xgamma;
  Int_t           det_cross_sec_bin_xp;
  Int_t           det_cross_sec_bin_dphi;
  Int_t           det_cross_sec_bin_deta;
  Int_t           det_cross_sec_bin_dphi_e_ph;
  Int_t           det_cross_sec_bin_deta_e_ph;

  Double_t        deltaRtrack;// NOT PRESENTED HERE 
  TString         period;
  TString         mc_type;
  TString         mc_corr_type;
  TString         root_file_name;
  TString         treeName;
  TFile           *f;
  TFile           *fCalib;
  TTree           *selectedTree;
  TDatime         time_begin;
  TDatime         time_end;
  TRandom3        *random;
  TLorentzVector  *v_uncorr_prompt_photon;
  TLorentzVector  *v_corr_prompt_photon;
  TLorentzVector  *v_true_prompt_photon;
  TLorentzVector  *v_true_jet_cont_prompt_photon;
  TLorentzVector  *v_prompt_photon_jet;
  TLorentzVector  *v_accomp_corr_jet;
  TLorentzVector  *v_accomp_uncorr_jet;
  TLorentzVector  *v_true_acc_jet;
  TLorentzVector  *v_true_parton_acc_jet;// THIS IS SOMETHING NEW is filled in SelectHadronLevel
  TLorentzVector *v_true_scattered_electron;// THIS IS SOMETHING NEW is filled in SelectHadronLevel
  Int_t           radiated_candidate_number;
  Int_t           max_et_candidate_number;
  Int_t           candidate_jet_number;
  //  Int_t       
  Double_t yjb;
  Double_t Empz;
  Double_t x_gamma;
  Double_t wtx;

  //cut values
  Double_t ET_JET_CUT;
  
  Int_t pcele_mode;
  Double_t bpres_conerad;
  //variables to write in a tree:
  Int_t w_Eventnr, w_Runnr;
  Int_t w_Fltw[2];
  Int_t w_Fltpsw[2];
  Int_t w_Sltw[6];
  Int_t w_Sltupw[6];
  Int_t w_Tltw[15];
  Int_t w_Dstw[4];
  Float_t         w_bosene;
  Float_t w_q2da, w_empz, w_yjb, w_eda, w_xbj, w_q2_true, w_y_true, w_x_true;
  Float_t r_q2da, r_empz, r_yjb, r_eda, r_xbj, r_q2_true, r_y_true, r_x_true;
  Int_t w_NofJets, w_NofHadronJets, w_NofPartonJets;
  Int_t r_NofJets, r_NofHadronJets, r_NofPartonJets;
  Float_t w_jets_Lab[maxNofJets][4], w_hadron_jets_Lab[maxNofJets][4], w_parton_jets_Lab[maxNofJets][4];
  Float_t r_jets_Lab[maxNofJets][4], r_hadron_jets_Lab[maxNofJets][4], r_parton_jets_Lab[maxNofJets][4];
  Float_t w_Scattered_Electron_DA[4], w_Scattered_electron_mc_true[4], w_Scattered_Electron[4];
  Float_t r_Scattered_Electron_DA[4], r_Scattered_electron_mc_true[4], r_Scattered_Electron[4];
  Float_t w_weight, r_weight;
  Int_t w_take_det_dijet_event;// = 0;
  Int_t w_take_hadron_dijet_event;// = 0;
  Int_t w_take_parton_dijet_event;// = 0;
  Int_t r_take_det_dijet_event;// = 0;
  Int_t r_take_hadron_dijet_event;// = 0;
  Int_t r_take_parton_dijet_event;// = 0;
  Double_t        C_scale[number_of_eta_bins];
  Double_t        a0_corr[number_of_eta_bins];
  Double_t        a1_corr[number_of_eta_bins];
  
  //declaration of histograms
  TH1D* h_events_per_runlumi_eltrigger_l;
  TH1D* h_events_per_runlumi_runs;
  TH1D* h_events_per_runlumi_largebins_lumi;

  TH1D* h_events_per_run;
  TH1D* h_dis_events_per_run;
  TH1D* h_dis_ev_per_run_without_cleaning;
  
  //cross section histos
  TH1D* h_Q2_det_crosssec;
  TH1D* h_eta_det_crosssec;
  TH1D* h_m_det_crosssec;
  TH1D* h_xi_det_crosssec;
  TH1D* h_xbj_det_crosssec;
  TH1D* h_etmean_det_crosssec;
  TProfile* tprof_Q2_det_crosssec;
  TProfile* tprof_eta_det_crosssec;
  TProfile* tprof_m_det_crosssec;
  TProfile* tprof_xi_det_crosssec;
  TProfile* tprof_xbj_det_crosssec;
  TProfile* tprof_etmean_det_crosssec;
 
  TH1D* h_Q2_had_crosssec;
  TH1D* h_eta_had_crosssec;
  TH1D* h_m_had_crosssec;
  TH1D* h_xi_had_crosssec;
  TH1D* h_xbj_had_crosssec;
  TH1D* h_etmean_had_crosssec;
  TProfile* tprof_Q2_had_crosssec;
  TProfile* tprof_eta_had_crosssec;
  TProfile* tprof_m_had_crosssec;
  TProfile* tprof_xi_had_crosssec;
  TProfile* tprof_xbj_had_crosssec;
  TProfile* tprof_etmean_had_crosssec;
 
  TH1D* h_Q2_part_crosssec;
  TH1D* h_eta_part_crosssec;
  TH1D* h_m_part_crosssec;
  TH1D* h_xi_part_crosssec;
  TH1D* h_xbj_part_crosssec;
  TH1D* h_etmean_part_crosssec;
  TProfile* tprof_Q2_part_crosssec;
  TProfile* tprof_eta_part_crosssec;
  TProfile* tprof_m_part_crosssec;
  TProfile* tprof_xi_part_crosssec;
  TProfile* tprof_xbj_part_crosssec;
  TProfile* tprof_etmean_part_crosssec;
 
  TH1D* h_Q2_hd_crosssec;
  TH1D* h_eta_hd_crosssec;
  TH1D* h_m_hd_crosssec;
  TH1D* h_xi_hd_crosssec;
  TH1D* h_xbj_hd_crosssec;
  TH1D* h_etmean_hd_crosssec;
 
  //double differential cross sections
  TH1D* h_xi_diff_q2_bins_1;
  TH1D* h_xi_diff_q2_bins_2;
  TH1D* h_xi_diff_q2_bins_3;
  TH1D* h_xi_diff_q2_bins_4;
  TH1D* h_xi_diff_q2_bins_5;
  TH1D* h_xi_diff_q2_bins_6;
  TProfile* tprof_xi_diff_q2_bins_1;
  TProfile* tprof_xi_diff_q2_bins_2;
  TProfile* tprof_xi_diff_q2_bins_3;
  TProfile* tprof_xi_diff_q2_bins_4;
  TProfile* tprof_xi_diff_q2_bins_5;
  TProfile* tprof_xi_diff_q2_bins_6;
 
  TH1D* h_had_xi_diff_q2_bins_1;
  TH1D* h_had_xi_diff_q2_bins_2;
  TH1D* h_had_xi_diff_q2_bins_3;
  TH1D* h_had_xi_diff_q2_bins_4;
  TH1D* h_had_xi_diff_q2_bins_5;
  TH1D* h_had_xi_diff_q2_bins_6;
 
  TH1D* h_hd_xi_diff_q2_bins_1;
  TH1D* h_hd_xi_diff_q2_bins_2;
  TH1D* h_hd_xi_diff_q2_bins_3;
  TH1D* h_hd_xi_diff_q2_bins_4;
  TH1D* h_hd_xi_diff_q2_bins_5;
  TH1D* h_hd_xi_diff_q2_bins_6;
 
  TH1D* h_etmean_diff_q2_bins_1;
  TH1D* h_etmean_diff_q2_bins_2;
  TH1D* h_etmean_diff_q2_bins_3;
  TH1D* h_etmean_diff_q2_bins_4;
  TH1D* h_etmean_diff_q2_bins_5;
  TH1D* h_etmean_diff_q2_bins_6;
  TProfile* tprof_etmean_diff_q2_bins_1;
  TProfile* tprof_etmean_diff_q2_bins_2;
  TProfile* tprof_etmean_diff_q2_bins_3;
  TProfile* tprof_etmean_diff_q2_bins_4;
  TProfile* tprof_etmean_diff_q2_bins_5;
  TProfile* tprof_etmean_diff_q2_bins_6;
  
  TH1D* h_had_etmean_diff_q2_bins_1;
  TH1D* h_had_etmean_diff_q2_bins_2;
  TH1D* h_had_etmean_diff_q2_bins_3;
  TH1D* h_had_etmean_diff_q2_bins_4;
  TH1D* h_had_etmean_diff_q2_bins_5;
  TH1D* h_had_etmean_diff_q2_bins_6;
  TH1D* h_hd_etmean_diff_q2_bins_1;
  TH1D* h_hd_etmean_diff_q2_bins_2;
  TH1D* h_hd_etmean_diff_q2_bins_3;
  TH1D* h_hd_etmean_diff_q2_bins_4;
  TH1D* h_hd_etmean_diff_q2_bins_5;
  TH1D* h_hd_etmean_diff_q2_bins_6;
  
  //trijets cross_sections
  TH1D* h_trijet_Q2_det_crosssec;
  TH1D* h_trijet_eta_det_crosssec;
  TH1D* h_trijet_m_det_crosssec;
  TH1D* h_trijet_xi_det_crosssec;
  TH1D* h_trijet_xbj_det_crosssec;
  TH1D* h_trijet_etmean_det_crosssec;
 
  TProfile* tprof_trijet_Q2_det_crosssec;
  TProfile* tprof_trijet_eta_det_crosssec;
  TProfile* tprof_trijet_m_det_crosssec;
  TProfile* tprof_trijet_xi_det_crosssec;
  TProfile* tprof_trijet_xbj_det_crosssec;
  TProfile* tprof_trijet_etmean_det_crosssec;
 
  TH1D* h_trijet_Q2_had_crosssec;
  TH1D* h_trijet_eta_had_crosssec;
  TH1D* h_trijet_m_had_crosssec;
  TH1D* h_trijet_xi_had_crosssec;
  TH1D* h_trijet_xbj_had_crosssec;
  TH1D* h_trijet_etmean_had_crosssec;
 
  TProfile* tprof_trijet_Q2_had_crosssec;
  TProfile* tprof_trijet_eta_had_crosssec;
  TProfile* tprof_trijet_m_had_crosssec;
  TProfile* tprof_trijet_xi_had_crosssec;
  TProfile* tprof_trijet_xbj_had_crosssec;
  TProfile* tprof_trijet_etmean_had_crosssec;
 
  TH1D* h_trijet_Q2_part_crosssec;
  TH1D* h_trijet_eta_part_crosssec;
  TH1D* h_trijet_m_part_crosssec;
  TH1D* h_trijet_xi_part_crosssec;
  TH1D* h_trijet_xbj_part_crosssec;
  TH1D* h_trijet_etmean_part_crosssec;
 
  TProfile* tprof_trijet_Q2_part_crosssec;
  TProfile* tprof_trijet_eta_part_crosssec;
  TProfile* tprof_trijet_m_part_crosssec;
  TProfile* tprof_trijet_xi_part_crosssec;
  TProfile* tprof_trijet_xbj_part_crosssec;
  TProfile* tprof_trijet_etmean_part_crosssec;
 
  TH1D* h_trijet_Q2_hd_crosssec;
  TH1D* h_trijet_eta_hd_crosssec;
  TH1D* h_trijet_m_hd_crosssec;
  TH1D* h_trijet_xi_hd_crosssec;
  TH1D* h_trijet_xbj_hd_crosssec;
  TH1D* h_trijet_etmean_hd_crosssec;
 
  TH1D* h_pJetsEtBreit;
  TH1D* h_hJetsEtBreit;
  TH1D* h_JetsEtBreit;
  TH1D* h_Zvtx;
  TH1D* h_E_da;
  TH1D* h_E_el;
  TH1D* h_yjb_da;
  TH1D* h_Q2;
  TH1D* h_pt_sqrtet;
  TH1D* h_hac_cell;
  TH1D* h_emc_cell;

  //dijet observables: detector level
  TH1D* h_q2_da;
  TH1D* h_M_jj;
  TH1D* h_log_ksi;
  TH1D* h_x_Bj;
  TH1D* h_jet_Et_mean_Breit;
  TH1D* h_Eta_str;
 
  TH1D* h_jet1_TBreit;
  TH1D* h_jet2_TBreit;
 
  TH1D* h_jet1_etaBreit;
  TH1D* h_jet2_etaBreit;
   
  //dijet observables: hadron level
  TH1D* h_had_q2_da;
  TH1D* h_had_M_jj;
  TH1D* h_had_log_ksi;
  TH1D* h_had_x_Bj;
  TH1D* h_had_jet_Et_mean_Breit;
  TH1D* h_had_Eta_str;
 
  TH1D* h_had_jet1_etaBreit;
  TH1D* h_had_jet2_etaBreit;
  TH1D* h_had_jet1_TBreit;
  TH1D* h_had_jet2_TBreit;


  //dijet observables: parton level
  TH1D* h_part_q2_da;
  TH1D* h_part_M_jj;
  TH1D* h_part_log_ksi;
  TH1D* h_part_x_Bj;
  TH1D* h_part_jet_Et_mean_Breit;
  TH1D* h_part_Eta_str;
 
  TH1D* h_part_jet1_etaBreit;
  TH1D* h_part_jet2_etaBreit;
  TH1D* h_part_jet1_TBreit;
  TH1D* h_part_jet2_TBreit;
 
  TH1D* h_Eel_Eda;
  TH1D* h_theta;
   
  // forward & backward jets
  TH1D* h_2jets_control_eta_b;
  TH1D* h_2jets_control_eta_f;
  TH1D* h_2jets_control_eta_breit_b;
  TH1D* h_2jets_control_eta_breit_f;
  TH1D* h_2jets_control_eta;
  TH1D* h_2jets_control_eta_breit;
   
  //dijets per run
  TH1D* h_dijets_per_run;
  
  //TH1D h_xvtx("xvtx", "xvtx", 40, -5., 5.);
  
  //electron position in RCAL sipos
  TH2D* h_elec_position;
  //h_elec_position->Sumw2();
  TH1D* h_elec_position_x;
  TH1D* h_elec_position_y;
  //h_elec_position_x->Sumw2();
  //h_elec_position_y->Sumw2();

  //q2 rekonstructed and real
  TH2D* h_q2_comp;
  TH1D* h_q2_deviation;
  
  //Fltw bits
  TH1D* h_fltw;
  TH1D* h_fltw_p;
  TH1D* h_fltw_dis;
  TH1D* h_fltw_dis_p;

  //jets phase space
  TH2D* h_et1et2;
  TH2D* h_et1et2_mcut;
  TH2D* h_et1et2_breit;
  TH2D* h_et1et2_mcut_breit;

  //for track veto correction
  TH1D* h_yda_veto_bit30;
  TH1D* h_yda_bit30;

  //Subjet variables
  TH1D* h_deltaR;

  //counters
  Long64_t phot_count;
  Long_t take_det_dijet;
  Long_t cutC_1;
  Long_t cutC_2;
  Long_t cutC_3;
  Long_t cutPS_1;
  Long_t cutPS_2;
  Long_t cutEC_1;
  Long_t cutEC_2;
  Long_t cutEC_3;
  Long_t cutEC_4;
  Long_t cutEC_5;
  Long_t cutEC_6_1;
  Long_t cutEC_6_2;
  Long_t cutEC_6_3;
  Long_t cutEC_6_4;
  Long_t cutElCand;
  Long_t cutCTDacc;
  Long_t cutJPS_1;
  Long_t cutJPS_23;
  Long_t cutJC_1;
  Long_t cutJC_2;
  Long_t cutJC_3;
  Long_t totalSelectedJetsCount;
  Long_t localSelectedJetsCount;
  Long_t goodCellsCount;
  Long_t selectedEventsCount;
  Long_t totalNumberOfEntries;
  Long_t processedEventsCount;
  Long_t cut_DIS03;
  Long_t cut_Sincand;
  Long_t cut_Siecorr;
  Long_t cut_q2_bos;
  Long_t cut_y_bos;
  Long_t cut_part_Et_jet;
  Long_t cut_part_eta;
  Long_t cut_hadr_Et_jet;
  Long_t cut_hadr_eta;
  Long_t cut_assymetric;
  Long_t cut_track_momentum;
  Long_t cut_electron_track_dca;
  Long_t cut_electron_distance_module_edge;
  Long_t cut_Mjj;
  Long_t partonDijetEvents;
  Long_t hadronDijetEvents;
  Long_t detector_dijets;
  Long_t cut_poltake;
  Long_t cut_QEDC;
  Long_t cut_yx2;
  Long_t cut_p_Mjj;
  Long_t cut_h_Mjj;
  Long_t cut_ngoodtrks;
  Long_t cut_gammahad;
  Long_t cut_RCAL;
  Long_t cut_triggers;
  Long_t count_nat;

  Float_t         myKzufodeltaz[250]; //to divide by 5.45 //IN THE WRONG PLACE
   
   // Declaration of leaf types
   Int_t           Runnr;
   Int_t           Eventnr;
   Int_t         Evtake;
   Int_t         Evtake_iwant;
   Int_t          Mbtake;
   Int_t           Fltw[2];
   Int_t           Fltpsw[2];
   Int_t           Fltfl;
   Int_t           Gslt_global;
   Int_t           Sltw[6];
   Int_t           Sltupw[6];
   Int_t           Tltw[15];
   Int_t           Dstw[4];
   Int_t           Fltpsfcw[2];
   Float_t         Cal_px;
   Float_t         Cal_py;
   Float_t         Cal_pz;
   Float_t         Cal_e;
   Float_t         Cal_et;
   Float_t         Cal_empz;
   Float_t         Cal_pt;
   Float_t         Cal_phi;
   Float_t         Remc[6];
   Float_t         Bemc[6];
   Float_t         Femc[6];
   Float_t         Rhac[6];
   Float_t         Bhac[6];
   Float_t         Fhac[6];
   Float_t         Bhac2[6];
   Float_t         Fhac2[6];
   Int_t           Nfemc;
   Int_t           Nfhac1;
   Int_t           Nfhac2;
   Int_t           Nbemc;
   Int_t           Nbhac1;
   Int_t           Nbhac2;
   Int_t           Nremc;
   Int_t           Nrhac;
   Float_t         Cal_tf;
   Float_t         Cal_tb;
   Float_t         Cal_tr;
   Float_t         Cal_tg;
   Float_t         Cal_tu;
   Float_t         Cal_td;
   Float_t         Cal_tf_e;
   Float_t         Cal_tb_e;
   Float_t         Cal_tr_e;
   Float_t         Cal_tg_e;
   Float_t         Cal_tu_e;
   Float_t         Cal_td_e;
   Int_t        Cal_tf_n;
   Int_t        Cal_tb_n;
   Int_t        Cal_tr_n;
   Int_t        Cal_tg_n;
   Int_t           Cal_tu_n;
   Int_t           Cal_td_n;
   Float_t         Etamax_ce;
   Float_t         Etamax_ce4;
   Float_t         Cal_et10;
   Float_t         Mtrknoe_pi;
   Float_t         Mtrknoe_k;
   Float_t         E_hk;
   Float_t         Unmen_pi;
   Float_t         Unmen_k;
   Int_t           Sparkf;
  // THIS IS SOMETHING NEW BLOCK 
      //Block: DisTrue
       Int_t         nlepton;
       Int_t         nradpho;
       Int_t         nboson;
       Int_t         nquark;
       Int_t         ngluon;
       Int_t        idscatlep;
       Int_t        idradpho;
       Int_t        idboson;
       Int_t        idquark;
       Int_t        idgluon;
       Int_t         dolepton;
       Int_t         doradpho;
       Int_t         doboson;
       Int_t         doquark;
       Int_t         dogluon;
       Float_t         plepton[5];
       Float_t         pradpho[5];
       Float_t         pboson[5];
       Float_t         pquark[5];
       Float_t         pgluon[5];
       Int_t           nqg;
       Int_t           quarkprt;
       Int_t        idqg[255];   //[nqg]
       Int_t         doqg[255];   //[nqg]
       Int_t           prtqg[255];   //[nqg]
       Float_t         pqg[255][5];   //[nqg]
   //----------------------------
   Int_t           E5ncand;
   Int_t           E5error;
   Float_t         E5prob[4];   //[E5ncand]
   Float_t         E5pos[4][3];   //[E5ncand]
   Float_t         E5calpos[4][3];   //[E5ncand]
   Float_t         E5calene[4];   //[E5ncand]
   Float_t         E5ein[4];   //[E5ncand]
   Float_t         E5enin[4];   //[E5ncand]
   Float_t         E5ecorr[4][3];   //[E5ncand]
   Float_t         E5th[4];   //[E5ncand]
   Float_t         E5ph[4];   //[E5ncand]
   Float_t         E5pt[4];   //[E5ncand]
   Int_t         E5xdet[4][3];   //[E5ncand]
   Int_t         E5ydet[4][3];   //[E5ncand]
   Int_t         E5trknr[4];   //[E5ncand]
   Int_t         E5nrsl[4];   //[E5ncand]
   Float_t         E5dca[4];   //[E5ncand]
   Float_t         E5dcabeam[4];   //[E5ncand]
   Float_t         E5trkp[4];   //[E5ncand]
   Float_t         E5trkth[4];   //[E5ncand]
   Float_t         E5trkph[4];   //[E5ncand]
   Float_t         E5trkq[4];   //[E5ncand]
   Float_t         E5trkdme[4];   //[E5ncand]
   Float_t         E5trkdce[4];   //[E5ncand]
   Float_t         E5trkpos[4][3];   //[E5ncand]
   Float_t         E5xel[4];   //[E5ncand]
   Float_t         E5yel[4];   //[E5ncand]
   Float_t         E5q2el[4];   //[E5ncand]
   Float_t         E5xda[4];   //[E5ncand]
   Float_t         E5yda[4];   //[E5ncand]
   Float_t         E5q2da[4];   //[E5ncand]
   Float_t         E5xda_cell[4];   //[E5ncand]
   Float_t         E5yda_cell[4];   //[E5ncand]
   Float_t         E5q2da_cell[4];   //[E5ncand]
   Float_t         E5xjb[4];   //[E5ncand]
   Float_t         E5yjb[4];   //[E5ncand]
   Float_t         E5q2jb[4];   //[E5ncand]
   Float_t         E5xjb_cell[4];   //[E5ncand]
   Float_t         E5yjb_cell[4];   //[E5ncand]
   Float_t         E5q2jb_cell[4];   //[E5ncand]
   Int_t         E5ncell[4];   //[E5ncand]
   Int_t        E5cellptr[4];   //[E5ncand]
   Float_t         E5femc[4];   //[E5ncand]
   Float_t         E5fcell[4][7];   //[E5ncand]
   Float_t         E5fmodu[4][5];   //[E5ncand]
   Float_t         E5imbmod[4][5];   //[E5ncand]
   Float_t         E5subprob[4][7];   //[E5ncand]
   Float_t         E5dphi[4];   //[E5ncand]
   Float_t         E5dtheta[4];   //[E5ncand]
   Float_t         E5calprob[4];   //[E5ncand]
   Int_t         E5calprobrank[4];   //[E5ncand]
   Float_t         E5fmaxbemc[4];   //[E5ncand]
   Float_t         E5fmaxremc[4];   //[E5ncand]
   Float_t         E5fmaxfemc[4];   //[E5ncand]
   Float_t         E5deltaz[4];   //[E5ncand]
   Float_t         E5deltax[4];   //[E5ncand]
   Float_t         E5deltay[4];   //[E5ncand]


   Float_t         Mc_q2;   
   Float_t         Mc_x;
   Float_t         Mc_y;

   //Block: FMCKin
    Int_t           Npart;
    Int_t           idlepton;
    Int_t           idphoton;
    Int_t           Part_prt[511];   //[Npart]
    Int_t           Part_id[511];   //[Npart]
    Float_t         Part_p[511][4];   //[Npart]
    Int_t           Part_motherid[511];   //[Npart]
    Int_t           Part_motherprt[511];   //[Npart]

  //Block: FMCZufo
       Int_t           Fmc_ezisl;
       Int_t           Fmc_gzisl;
       Int_t           Fmc_ezufo;
       Int_t           Fmc_gzufo;
       Int_t           Fmc_zisl[511];   //[Npart]
       Int_t           Fmc_zufo[511];   //[Npart]

   // THIS IS SOMETHING NEW BLOCK 
       Int_t           fmce5true[4];   //[E5ncand]
       Int_t           fmce5npart[4];   //[E5ncand]
       Int_t           fmce5partlist[4][5];   //[E5ncand]
       Float_t         fmce5enelist[4][5];   //[E5ncand]

       Float_t         Mc_ez;
       Float_t         Mc_esum;
       Float_t         Mc_etcone;
       Float_t         Mc_ercal;
       Float_t         Mc_el;
       Float_t         Mc_ep;

       Float_t         Mc_mu;
       Float_t         Mc_pt;
       Float_t         Mc_xpro;
       Float_t         Mc_xgam;
       Float_t         Mc_xpom;
       Float_t         Mc_beta;
       Float_t         Mc_t;
       Int_t           Mc_idl;
       Int_t           Mc_pidl;
       Int_t           Mc_pidp;
       Int_t           Mc_idfsl;
       Float_t         Mc_pfsl[4];
       Float_t         Mc_pfsph[4];
       Float_t         Mc_vtx[3];
       Int_t           Mc_ichnn;
       Float_t         Mc_q2_cr;
       Float_t         Mc_x_cr;
       Float_t         Mcvtx[3];
   //----------------------------
   Int_t           Sincand;
   Int_t           Sierror;
   Float_t         Siprob[4];   //[Sincand]
   Float_t         Sipos[4][3];   //[Sincand]
   Float_t         Sicalpos[4][3];   //[Sincand]
   Float_t         Sicalene[4];   //[Sincand]
   Float_t         Siein[4];   //[Sincand]
   Float_t         Sienin[4];   //[Sincand]
   Float_t         Siecorr[4][3];   //[Sincand]
   Float_t         Sith[4];   //[Sincand]
   Float_t         Siph[4];   //[Sincand]
   Float_t         Sipt[4];   //[Sincand]
   Int_t         Sixdet[4][3];   //[Sincand]
   Int_t         Siydet[4][3];   //[Sincand]
   Int_t         Sitrknr[4];   //[Sincand]
   Int_t         Sinrsl[4];   //[Sincand]
   Float_t         Sidca[4];   //[Sincand]
   Float_t         Sitrkp[4];   //[Sincand]
   Float_t         Sitrkth[4];   //[Sincand]
   Float_t         Sitrkph[4];   //[Sincand]
   Float_t         Sitrkq[4];   //[Sincand]
   Float_t         Sitrkdme[4];   //[Sincand]
   Float_t         Sitrkpos[4][3];   //[Sincand]
   Int_t           Sisrtf[4];   //[Sincand]
   Int_t         Sisrtquad[4];   //[Sincand]
   Int_t           Sihesf[4];   //[Sincand]
   Int_t         Sicorrcode[4];   //[Sincand]
   Float_t         Sisrtpos[4][2];   //[Sincand]
   Float_t         Sisrtene[4];   //[Sincand]
   Float_t         Sihespos[4][2];   //[Sincand]
   Float_t         Sihesene[4];   //[Sincand]
   Float_t         Sihesr[4];   //[Sincand]
   Float_t         Siprsene[4][3];   //[Sincand]
   Float_t         Siccet[4];   //[Sincand]
   Float_t         Siccempz[4];   //[Sincand]
   Float_t         Sietamax[4];   //[Sincand]
   Float_t         Sicehmom[4][4];   //[Sincand]
   Float_t         Sizuhmom[4][4];   //[Sincand]
   Float_t         Sicchmom[4][4];   //[Sincand]
   Float_t         Sixel[4];   //[Sincand]
   Float_t         Siyel[4];   //[Sincand]
   Float_t         Siq2el[4];   //[Sincand]
   Float_t         Sixda[4];   //[Sincand]
   Float_t         Siyda[4];   //[Sincand]
   Float_t         Siq2da[4];   //[Sincand]
   Float_t         Sixda_cell[4];   //[Sincand]
   Float_t         Siyda_cell[4];   //[Sincand]
   Float_t         Siq2da_cell[4];   //[Sincand]
   Float_t         Sixjb[4];   //[Sincand]
   Float_t         Siyjb[4];   //[Sincand]
   Float_t         Siq2jb[4];   //[Sincand]
   Float_t         Sixjb_cell[4];   //[Sincand]
   Float_t         Siyjb_cell[4];   //[Sincand]
   Float_t         Siq2jb_cell[4];   //[Sincand]

   // THIS IS SOMETHING NEW BLOCK 
       Int_t           fmcsitrue[4];   //[Sincand]
       Int_t           fmcsinpart[4];   //[Sincand]
       Int_t           fmcsipartlist[4][5];   //[Sincand]
       Float_t         fmcsienelist[4][5];   //[Sincand]
   //----------------------------

   Int_t           Nbpchn;
   Float_t         Bpmip[432];   //[Nbpchn]
   Float_t         Bpxyz[432][3];   //[Nbpchn]
   Int_t        Bpchn[432];   //[Nbpchn]
   Float_t         Bptim[432];   //[Nbpchn]
   Int_t           Ntrkvtx;
   Float_t         Xvtx;
   Float_t         Yvtx;
   Float_t         Zvtx;
   Float_t         Chivtx;
   Int_t           Nsecvtx;
   Float_t         Xsecvtx[40];   //[Nsecvtx]
   Float_t         Ysecvtx[40];   //[Nsecvtx]
   Float_t         Zsecvtx[40];   //[Nsecvtx]
   Float_t         Chisecvtx[40];   //[Nsecvtx]
   Float_t         Fetatr;
   Float_t         Betatr;
   Float_t         Pt_tr;
   Float_t         Empz_tr_pi;
   Float_t         Et_tr;
   Float_t         E_tr_pi;
   Float_t         phi_tr;
   Float_t         zvtx_fcal;
   Int_t        fcal_nrgoodcells;
   Float_t         fcal_vzerr;
   Float_t         V_h_px_zu;
   Float_t         V_h_py_zu;
   Float_t         V_h_pz_zu;
   Float_t         V_h_e_zu;
   Float_t         Etamax_zu;
   Float_t         Etamax_zu4;
   Float_t         Fgap;
   Float_t         Bgap;
   Int_t           Nzufos;
   //Nmu = 1000 - number of muon cadidate
   Int_t        Tufo[1000][4];   //[Nzufos]
   Int_t         zufo_bsp[1000];   //[Nzufos]
   Float_t         Zufo[1000][4];   //[Nzufos]
   Float_t         cufo[1000];   //[Nzufos]
   //Block: ZUFO_CAL
   Float_t         Zufoecal[1000];   //[Nzufos]
   Float_t         Zufoeemc[1000];   //[Nzufos]

   Float_t         zufopos[1000][3];   //[Nzufos]
   //trk_ntracks=400 total number of tracks in the tracking block [0,300]
   Int_t           Trk_ntracks;
   Int_t         trk_type[400];   //[Trk_ntracks]
   Int_t         ntrack_type[4];
   Int_t         def_trk_type;
   Int_t         trk_id[400];   //[Trk_ntracks]
   Int_t           trk_id2[400];   //[Trk_ntracks]
   Float_t         Trk_px[400];   //[Trk_ntracks]
   Float_t         Trk_py[400];   //[Trk_ntracks]
   Float_t         Trk_pz[400];   //[Trk_ntracks]
   Float_t         trk_charge[400];   //[Trk_ntracks]
   Char_t          trk_vtx[400];   //[Trk_ntracks]
   Int_t         Trk_prim_vtx[400];   //[Trk_ntracks]
   Int_t         trk_sec_vtx[400];   //[Trk_ntracks]
   Int_t         trk_vxid[400];   //[Trk_ntracks]
   Float_t         trk_endpos[400][3];   //[Trk_ntracks]
   Int_t         trk_nzbyt[400];   //[Trk_ntracks]
   Int_t         trk_naxial[400];   //[Trk_ntracks]
   Int_t         trk_nstereo[400];   //[Trk_ntracks]
   Int_t         trk_ndof[400];   //[Trk_ntracks]
   Int_t         trk_layinner[400];   //[Trk_ntracks]
   Int_t         trk_layouter[400];   //[Trk_ntracks]
   Float_t         trk_dedxctd[400];   //[Trk_ntracks]
   Float_t         trk_dedxctdcr[400];   //[Trk_ntracks]
   Short_t         trk_dedxctderr[400];   //[Trk_ntracks]
   Int_t         trk_dedxctdnh[400];   //[Trk_ntracks]
   Int_t         trk_dedxctdsaturated[400];   //[Trk_ntracks]
   Float_t         trk_chi2[400];   //[Trk_ntracks]
   Float_t         trk_vchi2[400];   //[Trk_ntracks]
   Float_t         trk_imppar[400];   //[Trk_ntracks]
   Float_t         trk_imperr[400];   //[Trk_ntracks]
   Float_t         trk_pca[400][3];   //[Trk_ntracks]
   Int_t           filter;

    //Block: QCDHAD
      Int_t           Nfmckin;
      Int_t           Idfmckin[511];   //[Nfmckin]
      Float_t         Ppfmckin[511][4];   //[Nfmckin]

    //Block: QCDBOSON
      Float_t         bospx;
      Float_t         bospy;
      Float_t         bospz;
      Float_t         bosene;

    //Block: QCDPAR
      Int_t           Nppart;
      Int_t           Idpart[511];   //[Nppart]
      Float_t         Ppart[511][4];   //[Nppart]

  // Not in ntuples08

       Int_t           npart_my;
       Int_t           Part_jetid[511];   //[npart_my]
       Int_t           Part_isthep[511];   //[npart_my]
       Int_t           Part_charge[511];   //[npart_my]
       Int_t           Photn;  //(!)
       Int_t           photprphjetindex;
       Int_t           Photid[100];   //[Photn] //(!)


  // THIS IS SOMETHING NEW BLOCK Not in ntuples08
       Float_t         Photp[100][4];   //[Photn]
       Int_t           photmothprt[100];   //[Photn]
       Int_t           photmothid[100];   //[Photn]
       Int_t           photisthep[100];   //[Photn]
  //----------------------------
  // Not in ntuples08

       Int_t           Knjets;
       Float_t         Kpjets[250][4];   //[Knjets]
       Float_t         Kpjetet[250];   //[Knjets]
       Float_t         Kpjetpt[250];   //[Knjets]
       Float_t         Kpjeteta[250];   //[Knjets]
       Float_t         Kpjetphi[250];   //[Knjets]
       Int_t         Kpjetnzu[250];   //[Knjets]
       Float_t         Kpjetemcfrac[250];   //[Knjets]
       Int_t         Kpjetnisl[250];   //[Knjets]
       Float_t         Kpjetfmax[250];   //[Knjets]
       Float_t         Kpjetdeltaz[250];   //[Knjets]
       Int_t           Ktrnjets;
       Float_t         Ktrjets[250][4];   //[Ktrnjets]
       Float_t         Ktrjetet[250];   //[Ktrnjets]
       Float_t         Ktrjetpt[250];   //[Ktrnjets]
       Float_t         Ktrjeteta[250];   //[Ktrnjets]
       Float_t         Ktrjetphi[250];   //[Ktrnjets]
       Int_t           Ktrjetntrac[250];   //[Ktrnjets]
       Int_t           Ktrjetchar[250];   //[Ktrnjets]
       Int_t           Knzufos;
       Float_t         Kzufos[1000][4];   //[Knzufos]
       Float_t         Kzufoet[1000];   //[Knzufos]
       Float_t         Kzufopt[1000];   //[Knzufos]
       Float_t         Kzufoeta[1000];   //[Knzufos]
       Float_t         Kzufophi[1000];   //[Knzufos]
       Float_t         Kzufoemcfrac[1000];   //[Knzufos]
       Float_t         Kzufofmax[1000];   //[Knzufos]
       Float_t         Kzufodeltaz[1000];   //[Knzufos]
       Int_t         Kzufotype[1000];   //[Knzufos]
       Int_t         Kzufoidjet[1000];   //[Knzufos]
       Int_t           Kzufoncells[1000];   //[Knzufos]

       Int_t           dattyp;
   //----------------------------

  // Not in ntuples08
       Int_t           bit3_tlt4;// Not in ntuples08
       Int_t           tlt4[32];// Not in ntuples08
   //----------------------------
   // THIS IS SOMETHING NEW BLOCK 
       Float_t         Siein03;
       Float_t         Sienin03;
       Int_t           mc_fla[8];
       Float_t         mc_juscal;
       Float_t         mc_jux1ge;
       Float_t         mc_jux2ge;
       Int_t           mc_np;
       Float_t         mc_pa[30][10];
       Float_t         mc_eparton[10];
       Int_t           mc_efla;
       Float_t         mc_ygen;
       Float_t         mc_xgen;
       Float_t         mc_q2gen;
       Float_t         mc_gaparton[10];
       Float_t         mc_yapp;
       Float_t         mc_xapp;
       Float_t         mc_q2app;
       Float_t         mc_xxt;
       Float_t         mc_yyt;
       Float_t         mc_zzt;
       Float_t         mc_gparton[10];
       Float_t         mc_pparton[8][10];
       Int_t           Mcebeam;
       Int_t           Mcpbeam;
       Float_t         Mcelec;
       Float_t         Mcethe;
       Float_t         Mcephi;
       Float_t         Mcq2;
       Float_t         Mcxbj;
       Float_t         Mcybj;
       Float_t         Mcgam[4];
       Float_t         vtxpos[3];


      Int_t Fmck_nstor;
      Int_t Fmck_id[511];
      Int_t Fmck_prt[511];
      Float_t Fmck_px[511];
      Float_t Fmck_py[511];
      Float_t Fmck_pz[511]; 
      Float_t Fmck_e[511];  
      Float_t Fmck_m[511];  
      Int_t Fmck_isthep[511];  
      Int_t Fmck_daug[511];  
   //----------------------------
   

   // List of branches

    TBranch        *b_Fmck_nstor;
    TBranch        *b_Fmck_id;
    TBranch        *b_Fmck_prt;
    TBranch        *b_Fmck_px;
    TBranch        *b_Fmck_py;
    TBranch        *b_Fmck_pz;
    TBranch        *b_Fmck_e;
    TBranch        *b_Fmck_m;
    TBranch        *b_Fmck_isthep;
    TBranch        *b_Fmck_daug;

   TBranch        *b_Runnr;   //!
   TBranch        *b_Eventnr;   //!
   TBranch        *b_Evtake;   //!
   TBranch        *b_Evtake_iwant;   //!
   TBranch        *b_Mbtake;   //!
   TBranch        *b_Fltw;   //!
   TBranch        *b_Fltpsw;   //!
   TBranch        *b_Fltfl;   //!
   TBranch        *b_Gslt_global;   //!
   TBranch        *b_Sltw;   //!
   TBranch        *b_Sltupw;   //!
   TBranch        *b_Tltw;   //!
   TBranch        *b_Dstw;   //!
   TBranch        *b_Fltpsfcw;   //!
   TBranch        *b_Cal_px;   //!
   TBranch        *b_Cal_py;   //!
   TBranch        *b_Cal_pz;   //!
   TBranch        *b_Cal_e;   //!
   TBranch        *b_Cal_et;   //!
   TBranch        *b_Cal_empz;   //!
   TBranch        *b_Cal_pt;   //!
   TBranch        *b_Cal_phi;   //!
   TBranch        *b_Remc;   //!
   TBranch        *b_Bemc;   //!
   TBranch        *b_Femc;   //!
   TBranch        *b_Rhac;   //!
   TBranch        *b_Bhac;   //!
   TBranch        *b_Fhac;   //!
   TBranch        *b_Bhac2;   //!
   TBranch        *b_Fhac2;   //!
   TBranch        *b_Nfemc;   //!
   TBranch        *b_Nfhac1;   //!
   TBranch        *b_Nfhac2;   //!
   TBranch        *b_Nbemc;   //!
   TBranch        *b_Nbhac1;   //!
   TBranch        *b_Nbhac2;   //!
   TBranch        *b_Nremc;   //!
   TBranch        *b_Nrhac;   //!
   TBranch        *b_Cal_tf;   //!
   TBranch        *b_Cal_tb;   //!
   TBranch        *b_Cal_tr;   //!
   TBranch        *b_Cal_tg;   //!
   TBranch        *b_Cal_tu;   //!
   TBranch        *b_Cal_td;   //!
   TBranch        *b_Cal_tf_e;   //!
   TBranch        *b_Cal_tb_e;   //!
   TBranch        *b_Cal_tr_e;   //!
   TBranch        *b_Cal_tg_e;   //!
   TBranch        *b_Cal_tu_e;   //!
   TBranch        *b_Cal_td_e;   //!
   TBranch        *b_Cal_tf_n;   //!
   TBranch        *b_Cal_tb_n;   //!
   TBranch        *b_Cal_tr_n;   //!
   TBranch        *b_Cal_tg_n;   //!
   TBranch        *b_Cal_tu_n;   //!
   TBranch        *b_Cal_td_n;   //!
   TBranch        *b_Etamax_ce;   //!
   TBranch        *b_Etamax_ce4;   //!
   TBranch        *b_Cal_et10;   //!
   TBranch        *b_Mtrknoe_pi;   //!
   TBranch        *b_Mtrknoe_k;   //!
   TBranch        *b_E_hk;   //!
   TBranch        *b_Unmen_pi;   //!
   TBranch        *b_Unmen_k;   //!
   TBranch        *b_Sparkf;   //!
   // THIS IS SOMETHING NEW BLOCK 
       TBranch        *b_nlepton;   //!
       TBranch        *b_nradpho;   //!
       TBranch        *b_nboson;   //!
       TBranch        *b_nquark;   //!
       TBranch        *b_ngluon;   //!
       TBranch        *b_idscatlep;   //!
       TBranch        *b_idradpho;   //!
       TBranch        *b_idboson;   //!
       TBranch        *b_idquark;   //!
       TBranch        *b_idgluon;   //!
       TBranch        *b_dolepton;   //!
       TBranch        *b_doradpho;   //!
       TBranch        *b_doboson;   //!
       TBranch        *b_doquark;   //!
       TBranch        *b_dogluon;   //!
       TBranch        *b_plepton;   //!
       TBranch        *b_pradpho;   //!
       TBranch        *b_pboson;   //!
       TBranch        *b_pquark;   //!
       TBranch        *b_pgluon;   //!
       TBranch        *b_nqg;   //!
       TBranch        *b_quarkprt;   //!
       TBranch        *b_idqg;   //!
       TBranch        *b_doqg;   //!
       TBranch        *b_prtqg;   //!
       TBranch        *b_pqg;   //!
    //----------------------------
   TBranch        *b_E5ncand;   //!
   TBranch        *b_E5error;   //!
   TBranch        *b_E5prob;   //!
   TBranch        *b_E5pos;   //!
   TBranch        *b_E5calpos;   //!
   TBranch        *b_E5calene;   //!
   TBranch        *b_E5ein;   //!
   TBranch        *b_E5enin;   //!
   TBranch        *b_E5ecorr;   //!
   TBranch        *b_E5th;   //!
   TBranch        *b_E5ph;   //!
   TBranch        *b_E5pt;   //!
   TBranch        *b_E5xdet;   //!
   TBranch        *b_E5ydet;   //!
   TBranch        *b_E5trknr;   //!
   TBranch        *b_E5nrsl;   //!
   TBranch        *b_E5dca;   //!
   TBranch        *b_E5dcabeam;   //!
   TBranch        *b_E5trkp;   //!
   TBranch        *b_E5trkth;   //!
   TBranch        *b_E5trkph;   //!
   TBranch        *b_E5trkq;   //!
   TBranch        *b_E5trkdme;   //!
   TBranch        *b_E5trkdce;   //!
   TBranch        *b_E5trkpos;   //!
   TBranch        *b_E5xel;   //!
   TBranch        *b_E5yel;   //!
   TBranch        *b_E5q2el;   //!
   TBranch        *b_E5xda;   //!
   TBranch        *b_E5yda;   //!
   TBranch        *b_E5q2da;   //!
   TBranch        *b_E5xda_cell;   //!
   TBranch        *b_E5yda_cell;   //!
   TBranch        *b_E5q2da_cell;   //!
   TBranch        *b_E5xjb;   //!
   TBranch        *b_E5yjb;   //!
   TBranch        *b_E5q2jb;   //!
   TBranch        *b_E5xjb_cell;   //!
   TBranch        *b_E5yjb_cell;   //!
   TBranch        *b_E5q2jb_cell;   //!
   TBranch        *b_E5ncell;   //!
   TBranch        *b_E5cellptr;   //!
   TBranch        *b_E5femc;   //!
   TBranch        *b_E5fcell;   //!
   TBranch        *b_E5fmodu;   //!
   TBranch        *b_E5imbmod;   //!
   TBranch        *b_E5subprob;   //!
   TBranch        *b_E5dphi;   //!
   TBranch        *b_E5dtheta;   //!
   TBranch        *b_E5calprob;   //!
   TBranch        *b_E5calprobrank;   //!
   TBranch        *b_E5fmaxbemc;   //!
   TBranch        *b_E5fmaxremc;   //!
   TBranch        *b_E5fmaxfemc;   //!
   TBranch        *b_E5deltaz;   //!
   TBranch        *b_E5deltax;   //!
   TBranch        *b_E5deltay;   //!

  //Block: FMCZufo
       TBranch        *b_Fmc_ezisl;   //!
       TBranch        *b_Fmc_gzisl;   //!
       TBranch        *b_Fmc_ezufo;   //!
       TBranch        *b_Fmc_gzufo;   //!
       TBranch        *b_Fmc_zisl;   //!
       TBranch        *b_Fmc_zufo;   //!

   // THIS IS SOMETHING NEW BLOCK 
       TBranch        *b_fmce5true;   //!
       TBranch        *b_fmce5npart;   //!
       TBranch        *b_fmce5partlist;   //!
       TBranch        *b_fmce5enelist;   //!
       TBranch        *b_npart;   //!
       TBranch        *b_idlepton;   //!
       TBranch        *b_idphoton;   //!
       TBranch        *b_Part_id;   //!
       TBranch        *b_Part_motherid;   //!
       TBranch        *b_Part_motherprt;   //!
       TBranch        *b_Part_prt;   //!
       TBranch        *b_Part_p;   //!
       TBranch        *b_Mc_ez;   //!
       TBranch        *b_Mc_esum;   //!
       TBranch        *b_Mc_etcone;   //!
       TBranch        *b_Mc_ercal;   //!
       TBranch        *b_Mc_el;   //!
       TBranch        *b_Mc_ep;   //!
       TBranch        *b_Mc_x;   //!
       TBranch        *b_Mc_y;   //!
       TBranch        *b_Mc_q2;   //!
       TBranch        *b_Mc_mu;   //!
       TBranch        *b_Mc_pt;   //!
       TBranch        *b_Mc_xpro;   //!
       TBranch        *b_Mc_xgam;   //!
       TBranch        *b_Mc_xpom;   //!
       TBranch        *b_Mc_beta;   //!
       TBranch        *b_Mc_t;   //!
       TBranch        *b_Mc_idl;   //!
       TBranch        *b_Mc_pidl;   //!
       TBranch        *b_Mc_pidp;   //!
       TBranch        *b_Mc_idfsl;   //!
       TBranch        *b_Mc_pfsl;   //!
       TBranch        *b_Mc_pfsph;   //!
       TBranch        *b_Mc_vtx;   //!
       TBranch        *b_Mc_ichnn;   //!
       TBranch        *b_Mc_q2_cr;   //!
       TBranch        *b_Mc_x_cr;   //!
       TBranch        *b_Mcvtx;   //!
    //----------------------------
   TBranch        *b_Sincand;   //!
   TBranch        *b_Sierror;   //!
   TBranch        *b_Siprob;   //!
   TBranch        *b_Sipos;   //!
   TBranch        *b_Sicalpos;   //!
   TBranch        *b_Sicalene;   //!
   TBranch        *b_Siein;   //!
   TBranch        *b_Sienin;   //!
   TBranch        *b_Siecorr;   //!
   TBranch        *b_Sith;   //!
   TBranch        *b_Siph;   //!
   TBranch        *b_Sipt;   //!
   TBranch        *b_Sixdet;   //!
   TBranch        *b_Siydet;   //!
   TBranch        *b_Sitrknr;   //!
   TBranch        *b_Sinrsl;   //!
   TBranch        *b_Sidca;   //!
   TBranch        *b_Sitrkp;   //!
   TBranch        *b_Sitrkth;   //!
   TBranch        *b_Sitrkph;   //!
   TBranch        *b_Sitrkq;   //!
   TBranch        *b_Sitrkdme;   //!
   TBranch        *b_Sitrkpos;   //!
   TBranch        *b_Sisrtf;   //!
   TBranch        *b_Sisrtquad;   //!
   TBranch        *b_Sihesf;   //!
   TBranch        *b_Sicorrcode;   //!
   TBranch        *b_Sisrtpos;   //!
   TBranch        *b_Sisrtene;   //!
   TBranch        *b_Sihespos;   //!
   TBranch        *b_Sihesene;   //!
   TBranch        *b_Sihesr;   //!
   TBranch        *b_Siprsene;   //!
   TBranch        *b_Siccet;   //!
   TBranch        *b_Siccempz;   //!
   TBranch        *b_Sietamax;   //!
   TBranch        *b_Sicehmom;   //!
   TBranch        *b_Sizuhmom;   //!
   TBranch        *b_Sicchmom;   //!
   TBranch        *b_Sixel;   //!
   TBranch        *b_Siyel;   //!
   TBranch        *b_Siq2el;   //!
   TBranch        *b_Sixda;   //!
   TBranch        *b_Siyda;   //!
   TBranch        *b_Siq2da;   //!
   TBranch        *b_Sixda_cell;   //!
   TBranch        *b_Siyda_cell;   //!
   TBranch        *b_Siq2da_cell;   //!
   TBranch        *b_Sixjb;   //!
   TBranch        *b_Siyjb;   //!
   TBranch        *b_Siq2jb;   //!
   TBranch        *b_Sixjb_cell;   //!
   TBranch        *b_Siyjb_cell;   //!
   TBranch        *b_Siq2jb_cell;   //!
   // THIS IS SOMETHING NEW BLOCK 
       TBranch        *b_fmcsitrue;   //!
       TBranch        *b_fmcsinpart;   //!
       TBranch        *b_fmcsipartlist;   //!
       TBranch        *b_fmcsienelist;   //!
       //----------------------------
   TBranch        *b_nbpchn;   //!
   TBranch        *b_Bpmip;   //!
   TBranch        *b_Bpxyz;   //!
   TBranch        *b_Bpchn;   //!
   TBranch        *b_Bptim;   //!
   TBranch        *b_Ntrkvtx;   //!
   TBranch        *b_Xvtx;   //!
   TBranch        *b_Yvtx;   //!
   TBranch        *b_Zvtx;   //!
   TBranch        *b_Chivtx;   //!
   TBranch        *b_Nsecvtx;   //!
   TBranch        *b_Xsecvtx;   //!
   TBranch        *b_Ysecvtx;   //!
   TBranch        *b_Zsecvtx;   //!
   TBranch        *b_Chisecvtx;   //!
   TBranch        *b_Fetatr;   //!
   TBranch        *b_Betatr;   //!
   TBranch        *b_Pt_tr;   //!
   TBranch        *b_Empz_tr_pi;   //!
   TBranch        *b_Et_tr;   //!
   TBranch        *b_E_tr_pi;   //!
   TBranch        *b_phi_tr;   //!
   TBranch        *b_zvtx_fcal;   //!
   TBranch        *b_fcal_nrgoodcells;   //!
   TBranch        *b_fcal_vzerr;   //!
   TBranch        *b_V_h_px_zu;   //!
   TBranch        *b_V_h_py_zu;   //!
   TBranch        *b_V_h_pz_zu;   //!
   TBranch        *b_V_h_e_zu;   //!
   TBranch        *b_Etamax_zu;   //!
   TBranch        *b_Etamax_zu4;   //!
   TBranch        *b_Fgap;   //!
   TBranch        *b_Bgap;   //!
   TBranch        *b_Nzufos;   //!
   TBranch        *b_tufo;   //!
   TBranch        *b_zufo_bsp;   //!
   TBranch        *b_zufo;   //!
   TBranch        *b_cufo;   //!
   TBranch        *b_zufoecal;   //!
   TBranch        *b_zufoeemc;   //!
   TBranch        *b_zufopos;   //!
   TBranch        *b_trk_ntracks;   //!
   TBranch        *b_trk_type;   //!
   TBranch        *b_ntrack_type;   //!
   TBranch        *b_def_trk_type;   //!
   TBranch        *b_trk_id;   //!
   TBranch        *b_trk_id2;   //!
   TBranch        *b_trk_px;   //!
   TBranch        *b_trk_py;   //!
   TBranch        *b_trk_pz;   //!
   TBranch        *b_trk_charge;   //!
   TBranch        *b_trk_vtx;   //!
   TBranch        *b_Trk_prim_vtx;   //!
   TBranch        *b_trk_sec_vtx;   //!
   TBranch        *b_trk_vxid;   //!
   TBranch        *b_trk_endpos;   //!
   TBranch        *b_trk_nzbyt;   //!
   TBranch        *b_trk_naxial;   //!
   TBranch        *b_trk_nstereo;   //!
   TBranch        *b_trk_ndof;   //!
   TBranch        *b_trk_layinner;   //!
   TBranch        *b_trk_layouter;   //!
   TBranch        *b_trk_dedxctd;   //!
   TBranch        *b_trk_dedxctdcr;   //!
   TBranch        *b_trk_dedxctderr;   //!
   TBranch        *b_trk_dedxctdnh;   //!
   TBranch        *b_trk_dedxctdsaturated;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_vchi2;   //!
   TBranch        *b_trk_imppar;   //!
   TBranch        *b_trk_imperr;   //!
   TBranch        *b_trk_pca;   //!
   TBranch        *b_filter;   //!
   // THIS IS SOMETHING NEW BLOCK 
       TBranch        *b_npart_my;   //!
       TBranch        *b_Part_jetid;   //!
       TBranch        *b_Part_isthep;   //!
       TBranch        *b_Part_charge;   //!
       TBranch        *b_Nfmckin;   //!
       TBranch        *b_Idfmckin;   //!
       TBranch        *b_Ppfmckin;   //!
       TBranch        *b_photn;   //!
       TBranch        *b_photprphjetindex;   //!
       TBranch        *b_photid;   //!
       TBranch        *b_photp;   //!
       TBranch        *b_photmothprt;   //!
       TBranch        *b_photmothid;   //!
       TBranch        *b_photisthep;   //!

       TBranch        *b_Knjets;   //!
       TBranch        *b_Kpjets;   //!
       TBranch        *b_Kpjetet;   //!
       TBranch        *b_Kpjetpt;   //!
       TBranch        *b_Kpjeteta;   //!
       TBranch        *b_Kpjetphi;   //!
       TBranch        *b_Kpjetnzu;   //!
       TBranch        *b_Kpjetemcfrac;   //!
       TBranch        *b_Kpjetnisl;   //!
       TBranch        *b_Kpjetfmax;   //!
       TBranch        *b_Kpjetdeltaz;   //!
   TBranch        *b_Ktrnjets;   //!
   TBranch        *b_Ktrjets;   //!
       TBranch        *b_Ktrjetet;   //!
       TBranch        *b_Ktrjetpt;   //!
       TBranch        *b_Ktrjeteta;   //!
       TBranch        *b_Ktrjetphi;   //!
       TBranch        *b_Ktrjetntrac;   //!
       TBranch        *b_Ktrjetchar;   //!
       TBranch        *b_Knzufos;   //!
       TBranch        *b_Kzufos;   //!
   TBranch        *b_Kzufoet;   //!
   TBranch        *b_Kzufopt;   //!
   TBranch        *b_Kzufoeta;   //!
   TBranch        *b_Kzufophi;   //!
   TBranch        *b_Kzufoemcfrac;   //!
   TBranch        *b_Kzufofmax;   //!
   TBranch        *b_Kzufodeltaz;   //!
   TBranch        *b_Kzufotype;   //!
   TBranch        *b_Kzufoidjet;   //!
   TBranch        *b_Kzufoncells;   //!
       TBranch        *b_dattyp;   //!
       TBranch        *b_bospx;   //!
       TBranch        *b_bospy;   //!
       TBranch        *b_bospz;   //!
       TBranch        *b_bosene;   //!
       TBranch        *b_bit3_tlt4;   //!
       TBranch        *b_tlt4;   //!
       TBranch        *b_Nppart;   //!
       TBranch        *b_Idpart;   //!
       TBranch        *b_Ppart;   //!
       TBranch        *b_Siein03;   //!
       TBranch        *b_Sienin03;   //!
       TBranch        *b_mc_fla;   //!
       TBranch        *b_mc_juscal;   //!
       TBranch        *b_mc_jux1ge;   //!
       TBranch        *b_mc_jux2ge;   //!
       TBranch        *b_mc_np;   //!
       TBranch        *b_mc_pa;   //!
       TBranch        *b_mc_eparton;   //!
       TBranch        *b_mc_efla;   //!
       TBranch        *b_mc_ygen;   //!
       TBranch        *b_mc_xgen;   //!
       TBranch        *b_mc_q2gen;   //!
       TBranch        *b_mc_gaparton;   //!
       TBranch        *b_mc_yapp;   //!
       TBranch        *b_mc_xapp;   //!
       TBranch        *b_mc_q2app;   //!
       TBranch        *b_mc_xxt;   //!
       TBranch        *b_mc_yyt;   //!
       TBranch        *b_mc_zzt;   //!
       TBranch        *b_mc_gparton;   //!
       TBranch        *b_mc_pparton;   //!
       TBranch        *b_Mcebeam;   //!
       TBranch        *b_Mcpbeam;   //!
       TBranch        *b_Mcelec;   //!
       TBranch        *b_Mcethe;   //!
       TBranch        *b_Mcephi;   //!
       TBranch        *b_Mcq2;   //!
       TBranch        *b_Mcxbj;   //!
       TBranch        *b_Mcybj;   //!
       TBranch        *b_Mcgam;   //!
       TBranch        *b_vtxpos;   //!
   //----------------------------
   //   selector(TTree * /*tree*/ =0) { }
   virtual ~selector() {cout <<"destructor" << endl; }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(/*TTree *tree*/);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree, TString run_period, Bool_t b_Data, TString mc_type, TString mc_corr_type, Bool_t b_usecorr, Bool_t b_use2ndcorr, Bool_t b_use_clustered);
  virtual Bool_t  Notify();
  //   virtual Bool_t  Process(Long64_t entry);
  virtual Bool_t Process();
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  //   virtual void    SetOption(const char *option) { fOption = option; }
  //   virtual void    SetObject(TObject *obj) { fObject = obj; }
  //   virtual void    SetInputList(TList *input) { fInput = input; }
  //   virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
   
  //   ClassDef(selector,0);
};

#endif

#ifdef selector_c
//#ifndef selector_init
//#define selector_init
void selector::Init(TTree *tree, TString run_period, Bool_t b_Data, TString s_mc_type, TString s_mc_corr_type, Bool_t b_usecorr, Bool_t b_use2ndcorr, Bool_t b_use_clustered)
{
  wtx = 1;
  check_en_mom_conservation = true;
  check_en_mom_conservation_on_parton_level = false;
  en_mom_conservation = true;
  nodebugmode = kFALSE;
  event_list = new TEventList;
  period = run_period;
  Data = b_Data;
  use_corr = b_usecorr;
  take_2ndana_corr = b_use2ndcorr;
  use_clustered = b_use_clustered;
  mc_type = s_mc_type;
  mc_corr_type = s_mc_corr_type;
  cout <<"initializing....\nperiod = " << period << "\nData = " << Data << "\nmc_corr_type = " << mc_corr_type
       << "\nuse_corr = " << use_corr 
       << "\ntake_2ndana_corr = " << take_2ndana_corr << "\nuse_clustered = " << use_clustered << endl;
  if(!Data)
    cout << "!Data: mc_type = " << mc_type << "\nmc_corr_type = " << mc_corr_type << endl;
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  cout << "in selector.h: added new chain with " << tree->GetEntries() << " entries " << endl;
  fChain->SetMakeClass(1);

  //I added
        //Block: QCDHAD
         fChain->SetBranchAddress("Nfmckin", &Nfmckin, &b_Nfmckin);
         fChain->SetBranchAddress("Idfmckin", Idfmckin, &b_Idfmckin);
         fChain->SetBranchAddress("Ppfmckin", Ppfmckin, &b_Ppfmckin);

   fChain->SetBranchAddress("Runnr", &Runnr, &b_Runnr);
   fChain->SetBranchAddress("Eventnr", &Eventnr, &b_Eventnr);
   fChain->SetBranchAddress("Evtake", &Evtake, &b_Evtake);
   fChain->SetBranchAddress("Evtake_iwant", &Evtake_iwant, &b_Evtake_iwant);
   fChain->SetBranchAddress("Mbtake", &Mbtake, &b_Mbtake);
   fChain->SetBranchAddress("Fltw", Fltw, &b_Fltw);
   fChain->SetBranchAddress("Fltpsw", Fltpsw, &b_Fltpsw);
   fChain->SetBranchAddress("Fltfl", &Fltfl, &b_Fltfl);
   fChain->SetBranchAddress("Gslt_global", &Gslt_global, &b_Gslt_global);
   fChain->SetBranchAddress("Sltw", Sltw, &b_Sltw);
   fChain->SetBranchAddress("Sltupw", Sltupw, &b_Sltupw);
   fChain->SetBranchAddress("Tltw", Tltw, &b_Tltw);
   fChain->SetBranchAddress("Dstw", Dstw, &b_Dstw);
   fChain->SetBranchAddress("Fltpsfcw", Fltpsfcw, &b_Fltpsfcw);
   fChain->SetBranchAddress("Cal_px", &Cal_px, &b_Cal_px);
   fChain->SetBranchAddress("Cal_py", &Cal_py, &b_Cal_py);
   fChain->SetBranchAddress("Cal_pz", &Cal_pz, &b_Cal_pz);
   fChain->SetBranchAddress("Cal_e", &Cal_e, &b_Cal_e);
   fChain->SetBranchAddress("Cal_et", &Cal_et, &b_Cal_et);
   fChain->SetBranchAddress("Cal_empz", &Cal_empz, &b_Cal_empz);
   fChain->SetBranchAddress("Cal_pt", &Cal_pt, &b_Cal_pt);
   fChain->SetBranchAddress("Cal_phi", &Cal_phi, &b_Cal_phi);
   fChain->SetBranchAddress("Remc", Remc, &b_Remc);
   fChain->SetBranchAddress("Bemc", Bemc, &b_Bemc);
   fChain->SetBranchAddress("Femc", Femc, &b_Femc);
   fChain->SetBranchAddress("Rhac", Rhac, &b_Rhac);
   fChain->SetBranchAddress("Bhac", Bhac, &b_Bhac);
   fChain->SetBranchAddress("Fhac", Fhac, &b_Fhac);
   fChain->SetBranchAddress("Bhac2", Bhac2, &b_Bhac2);
   fChain->SetBranchAddress("Fhac2", Fhac2, &b_Fhac2);
   fChain->SetBranchAddress("Nfemc", &Nfemc, &b_Nfemc);
   fChain->SetBranchAddress("Nfhac1", &Nfhac1, &b_Nfhac1);
   fChain->SetBranchAddress("Nfhac2", &Nfhac2, &b_Nfhac2);
   fChain->SetBranchAddress("Nbemc", &Nbemc, &b_Nbemc);
   fChain->SetBranchAddress("Nbhac1", &Nbhac1, &b_Nbhac1);
   fChain->SetBranchAddress("Nbhac2", &Nbhac2, &b_Nbhac2);
   fChain->SetBranchAddress("Nremc", &Nremc, &b_Nremc);
   fChain->SetBranchAddress("Nrhac", &Nrhac, &b_Nrhac);
   fChain->SetBranchAddress("Cal_tf", &Cal_tf, &b_Cal_tf);
   fChain->SetBranchAddress("Cal_tb", &Cal_tb, &b_Cal_tb);
   fChain->SetBranchAddress("Cal_tr", &Cal_tr, &b_Cal_tr);
   fChain->SetBranchAddress("Cal_tg", &Cal_tg, &b_Cal_tg);
   fChain->SetBranchAddress("Cal_tu", &Cal_tu, &b_Cal_tu);
   fChain->SetBranchAddress("Cal_td", &Cal_td, &b_Cal_td);
   fChain->SetBranchAddress("Cal_tf_e", &Cal_tf_e, &b_Cal_tf_e);
   fChain->SetBranchAddress("Cal_tb_e", &Cal_tb_e, &b_Cal_tb_e);
   fChain->SetBranchAddress("Cal_tr_e", &Cal_tr_e, &b_Cal_tr_e);
   fChain->SetBranchAddress("Cal_tg_e", &Cal_tg_e, &b_Cal_tg_e);
   fChain->SetBranchAddress("Cal_tu_e", &Cal_tu_e, &b_Cal_tu_e);
   fChain->SetBranchAddress("Cal_td_e", &Cal_td_e, &b_Cal_td_e);
   fChain->SetBranchAddress("Cal_tf_n", &Cal_tf_n, &b_Cal_tf_n);
   fChain->SetBranchAddress("Cal_tb_n", &Cal_tb_n, &b_Cal_tb_n);
   fChain->SetBranchAddress("Cal_tr_n", &Cal_tr_n, &b_Cal_tr_n);
   fChain->SetBranchAddress("Cal_tg_n", &Cal_tg_n, &b_Cal_tg_n);
   fChain->SetBranchAddress("Cal_tu_n", &Cal_tu_n, &b_Cal_tu_n);
   fChain->SetBranchAddress("Cal_td_n", &Cal_td_n, &b_Cal_td_n);
   fChain->SetBranchAddress("Etamax_ce", &Etamax_ce, &b_Etamax_ce);
   fChain->SetBranchAddress("Etamax_ce4", &Etamax_ce4, &b_Etamax_ce4);
   fChain->SetBranchAddress("Cal_et10", &Cal_et10, &b_Cal_et10);
   fChain->SetBranchAddress("Mtrknoe_pi", &Mtrknoe_pi, &b_Mtrknoe_pi);
   fChain->SetBranchAddress("Mtrknoe_k", &Mtrknoe_k, &b_Mtrknoe_k);
   fChain->SetBranchAddress("E_hk", &E_hk, &b_E_hk);
   fChain->SetBranchAddress("Unmen_pi", &Unmen_pi, &b_Unmen_pi);
   fChain->SetBranchAddress("Unmen_k", &Unmen_k, &b_Unmen_k);
   fChain->SetBranchAddress("Sparkf", &Sparkf, &b_Sparkf);

   // THIS IS SOMETHING NEW BLOCK 
      //Error in <TTree::SetBranchStatus>: unknown branch -> 
        /*
         fChain->SetBranchAddress("nlepton", &nlepton, &b_nlepton);
         fChain->SetBranchAddress("nradpho", &nradpho, &b_nradpho);
         fChain->SetBranchAddress("nboson", &nboson, &b_nboson);
         fChain->SetBranchAddress("nquark", &nquark, &b_nquark);
         fChain->SetBranchAddress("ngluon", &ngluon, &b_ngluon);
         fChain->SetBranchAddress("idscatlep", &idscatlep, &b_idscatlep);
         fChain->SetBranchAddress("idradpho", &idradpho, &b_idradpho);
         fChain->SetBranchAddress("idboson", &idboson, &b_idboson);
         fChain->SetBranchAddress("idquark", &idquark, &b_idquark);
         fChain->SetBranchAddress("idgluon", &idgluon, &b_idgluon);
         fChain->SetBranchAddress("dolepton", &dolepton, &b_dolepton);
         fChain->SetBranchAddress("doradpho", &doradpho, &b_doradpho);
         fChain->SetBranchAddress("doboson", &doboson, &b_doboson);
         fChain->SetBranchAddress("doquark", &doquark, &b_doquark);
         fChain->SetBranchAddress("dogluon", &dogluon, &b_dogluon);
         fChain->SetBranchAddress("plepton", plepton, &b_plepton);
         fChain->SetBranchAddress("pradpho", pradpho, &b_pradpho);
         fChain->SetBranchAddress("pboson", pboson, &b_pboson);
         fChain->SetBranchAddress("pquark", pquark, &b_pquark);
         fChain->SetBranchAddress("pgluon", pgluon, &b_pgluon);
         fChain->SetBranchAddress("nqg", &nqg, &b_nqg);
         fChain->SetBranchAddress("quarkprt", &quarkprt, &b_quarkprt);
         fChain->SetBranchAddress("idqg", idqg, &b_idqg);
         fChain->SetBranchAddress("doqg", doqg, &b_doqg);
         fChain->SetBranchAddress("prtqg", prtqg, &b_prtqg);
         fChain->SetBranchAddress("pqg", pqg, &b_pqg);
         */
  //---------------------------

   fChain->SetBranchAddress("E5ncand", &E5ncand, &b_E5ncand);
   fChain->SetBranchAddress("E5error", &E5error, &b_E5error);
   fChain->SetBranchAddress("E5prob", E5prob, &b_E5prob);
   fChain->SetBranchAddress("E5pos", E5pos, &b_E5pos);
   fChain->SetBranchAddress("E5calpos", E5calpos, &b_E5calpos);
   fChain->SetBranchAddress("E5calene", E5calene, &b_E5calene);
   fChain->SetBranchAddress("E5ein", E5ein, &b_E5ein);
   fChain->SetBranchAddress("E5enin", E5enin, &b_E5enin);
   fChain->SetBranchAddress("E5ecorr", E5ecorr, &b_E5ecorr);
   fChain->SetBranchAddress("E5th", E5th, &b_E5th);
   fChain->SetBranchAddress("E5ph", E5ph, &b_E5ph);
   fChain->SetBranchAddress("E5pt", E5pt, &b_E5pt);
   fChain->SetBranchAddress("E5xdet", E5xdet, &b_E5xdet);
   fChain->SetBranchAddress("E5ydet", E5ydet, &b_E5ydet);
   fChain->SetBranchAddress("E5trknr", E5trknr, &b_E5trknr);
   fChain->SetBranchAddress("E5nrsl", E5nrsl, &b_E5nrsl);
   fChain->SetBranchAddress("E5dca", E5dca, &b_E5dca);
   fChain->SetBranchAddress("E5dcabeam", E5dcabeam, &b_E5dcabeam);
   fChain->SetBranchAddress("E5trkp", E5trkp, &b_E5trkp);
   fChain->SetBranchAddress("E5trkth", E5trkth, &b_E5trkth);
   fChain->SetBranchAddress("E5trkph", E5trkph, &b_E5trkph);
   fChain->SetBranchAddress("E5trkq", E5trkq, &b_E5trkq);
   fChain->SetBranchAddress("E5trkdme", E5trkdme, &b_E5trkdme);
   fChain->SetBranchAddress("E5trkdce", E5trkdce, &b_E5trkdce);
   fChain->SetBranchAddress("E5trkpos", E5trkpos, &b_E5trkpos);
   fChain->SetBranchAddress("E5xel", E5xel, &b_E5xel);
   fChain->SetBranchAddress("E5yel", E5yel, &b_E5yel);
   fChain->SetBranchAddress("E5q2el", E5q2el, &b_E5q2el);
   fChain->SetBranchAddress("E5xda", E5xda, &b_E5xda);
   fChain->SetBranchAddress("E5yda", E5yda, &b_E5yda);
   fChain->SetBranchAddress("E5q2da", E5q2da, &b_E5q2da);
   fChain->SetBranchAddress("E5xda_cell", E5xda_cell, &b_E5xda_cell);
   fChain->SetBranchAddress("E5yda_cell", E5yda_cell, &b_E5yda_cell);
   fChain->SetBranchAddress("E5q2da_cell", E5q2da_cell, &b_E5q2da_cell);
   fChain->SetBranchAddress("E5xjb", E5xjb, &b_E5xjb);
   fChain->SetBranchAddress("E5yjb", E5yjb, &b_E5yjb);
   fChain->SetBranchAddress("E5q2jb", E5q2jb, &b_E5q2jb);
   fChain->SetBranchAddress("E5xjb_cell", E5xjb_cell, &b_E5xjb_cell);
   fChain->SetBranchAddress("E5yjb_cell", E5yjb_cell, &b_E5yjb_cell);
   fChain->SetBranchAddress("E5q2jb_cell", E5q2jb_cell, &b_E5q2jb_cell);
   fChain->SetBranchAddress("E5ncell", E5ncell, &b_E5ncell);
   fChain->SetBranchAddress("E5cellptr", E5cellptr, &b_E5cellptr);
   fChain->SetBranchAddress("E5femc", E5femc, &b_E5femc);
   fChain->SetBranchAddress("E5fcell", E5fcell, &b_E5fcell);
   fChain->SetBranchAddress("E5fmodu", E5fmodu, &b_E5fmodu);
   fChain->SetBranchAddress("E5imbmod", E5imbmod, &b_E5imbmod);
   fChain->SetBranchAddress("E5subprob", E5subprob, &b_E5subprob);
   fChain->SetBranchAddress("E5dphi", E5dphi, &b_E5dphi);
   fChain->SetBranchAddress("E5dtheta", E5dtheta, &b_E5dtheta);
   fChain->SetBranchAddress("E5calprob", E5calprob, &b_E5calprob);
   fChain->SetBranchAddress("E5calprobrank", E5calprobrank, &b_E5calprobrank);
   fChain->SetBranchAddress("E5fmaxbemc", E5fmaxbemc, &b_E5fmaxbemc);
   fChain->SetBranchAddress("E5fmaxremc", E5fmaxremc, &b_E5fmaxremc);
   fChain->SetBranchAddress("E5fmaxfemc", E5fmaxfemc, &b_E5fmaxfemc);
   fChain->SetBranchAddress("E5deltaz", E5deltaz, &b_E5deltaz);
   fChain->SetBranchAddress("E5deltax", E5deltax, &b_E5deltax);
   fChain->SetBranchAddress("E5deltay", E5deltay, &b_E5deltay);

   // THIS IS SOMETHING NEW BLOCK 
  //---------------------------
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmce5true", fmce5true, &b_fmce5true);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmce5npart", fmce5npart, &b_fmce5npart);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmce5partlist", fmce5partlist, &b_fmce5partlist);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmce5enelist", fmce5enelist, &b_fmce5enelist);
   //Block: FMCZufo
      fChain->SetBranchAddress("Fmc_ezisl", &Fmc_ezisl, &b_Fmc_ezisl);
      fChain->SetBranchAddress("Fmc_gzisl", &Fmc_gzisl, &b_Fmc_gzisl);
      fChain->SetBranchAddress("Fmc_ezufo", &Fmc_ezufo, &b_Fmc_ezufo);
      fChain->SetBranchAddress("Fmc_gzufo", &Fmc_gzufo, &b_Fmc_gzufo);
      fChain->SetBranchAddress("Fmc_zisl", Fmc_zisl, &b_Fmc_zisl);
      fChain->SetBranchAddress("Fmc_zufo", Fmc_zufo, &b_Fmc_zufo);

       fChain->SetBranchAddress("Mc_ez", &Mc_ez, &b_Mc_ez);
       fChain->SetBranchAddress("Mc_esum", &Mc_esum, &b_Mc_esum);
       fChain->SetBranchAddress("Mc_etcone", &Mc_etcone, &b_Mc_etcone);
       fChain->SetBranchAddress("Mc_ercal", &Mc_ercal, &b_Mc_ercal);
       fChain->SetBranchAddress("Mc_el", &Mc_el, &b_Mc_el);
       fChain->SetBranchAddress("Mc_ep", &Mc_ep, &b_Mc_ep);
       fChain->SetBranchAddress("Mc_x", &Mc_x, &b_Mc_x);
       fChain->SetBranchAddress("Mc_y", &Mc_y, &b_Mc_y);
       fChain->SetBranchAddress("Mc_q2", &Mc_q2, &b_Mc_q2);
       fChain->SetBranchAddress("Mc_mu", &Mc_mu, &b_Mc_mu);
       fChain->SetBranchAddress("Mc_pt", &Mc_pt, &b_Mc_pt);
       fChain->SetBranchAddress("Mc_xpro", &Mc_xpro, &b_Mc_xpro);
       fChain->SetBranchAddress("Mc_xgam", &Mc_xgam, &b_Mc_xgam);
       fChain->SetBranchAddress("Mc_xpom", &Mc_xpom, &b_Mc_xpom);
       fChain->SetBranchAddress("Mc_beta", &Mc_beta, &b_Mc_beta);
       fChain->SetBranchAddress("Mc_t", &Mc_t, &b_Mc_t);
       fChain->SetBranchAddress("Mc_idl", &Mc_idl, &b_Mc_idl);
       fChain->SetBranchAddress("Mc_pidl", &Mc_pidl, &b_Mc_pidl);
       fChain->SetBranchAddress("Mc_pidp", &Mc_pidp, &b_Mc_pidp);
       fChain->SetBranchAddress("Mc_idfsl", &Mc_idfsl, &b_Mc_idfsl);
       fChain->SetBranchAddress("Mc_pfsl", Mc_pfsl, &b_Mc_pfsl);
       fChain->SetBranchAddress("Mc_pfsph", Mc_pfsph, &b_Mc_pfsph);
       fChain->SetBranchAddress("Mc_vtx", Mc_vtx, &b_Mc_vtx);
       fChain->SetBranchAddress("Mc_ichnn", &Mc_ichnn, &b_Mc_ichnn);
       fChain->SetBranchAddress("Mc_q2_cr", &Mc_q2_cr, &b_Mc_q2_cr);
       fChain->SetBranchAddress("Mc_x_cr", &Mc_x_cr, &b_Mc_x_cr);
       fChain->SetBranchAddress("Mcvtx", Mcvtx, &b_Mcvtx);
   // THIS IS SOMETHING NEW BLOCK 
  //---------------------------

   fChain->SetBranchAddress("Sincand", &Sincand, &b_Sincand);
   fChain->SetBranchAddress("Sierror", &Sierror, &b_Sierror);
   fChain->SetBranchAddress("Siprob", Siprob, &b_Siprob);
   fChain->SetBranchAddress("Sipos", Sipos, &b_Sipos);
   fChain->SetBranchAddress("Sicalpos", Sicalpos, &b_Sicalpos);
   fChain->SetBranchAddress("Sicalene", Sicalene, &b_Sicalene);
   fChain->SetBranchAddress("Siein", Siein, &b_Siein);
   fChain->SetBranchAddress("Sienin", Sienin, &b_Sienin);
   fChain->SetBranchAddress("Siecorr", Siecorr, &b_Siecorr);
   fChain->SetBranchAddress("Sith", Sith, &b_Sith);
   fChain->SetBranchAddress("Siph", Siph, &b_Siph);
   fChain->SetBranchAddress("Sipt", Sipt, &b_Sipt);
   fChain->SetBranchAddress("Sixdet", Sixdet, &b_Sixdet);
   fChain->SetBranchAddress("Siydet", Siydet, &b_Siydet);
   fChain->SetBranchAddress("Sitrknr", Sitrknr, &b_Sitrknr);
   fChain->SetBranchAddress("Sinrsl", Sinrsl, &b_Sinrsl);
   fChain->SetBranchAddress("Sidca", Sidca, &b_Sidca);
   fChain->SetBranchAddress("Sitrkp", Sitrkp, &b_Sitrkp);
   fChain->SetBranchAddress("Sitrkth", Sitrkth, &b_Sitrkth);
   fChain->SetBranchAddress("Sitrkph", Sitrkph, &b_Sitrkph);
   fChain->SetBranchAddress("Sitrkq", Sitrkq, &b_Sitrkq);
   fChain->SetBranchAddress("Sitrkdme", Sitrkdme, &b_Sitrkdme);
   fChain->SetBranchAddress("Sitrkpos", Sitrkpos, &b_Sitrkpos);
   fChain->SetBranchAddress("Sisrtf", Sisrtf, &b_Sisrtf);
   fChain->SetBranchAddress("Sisrtquad", Sisrtquad, &b_Sisrtquad);
   fChain->SetBranchAddress("Sihesf", Sihesf, &b_Sihesf);
   fChain->SetBranchAddress("Sicorrcode", Sicorrcode, &b_Sicorrcode);
   fChain->SetBranchAddress("Sisrtpos", Sisrtpos, &b_Sisrtpos);
   fChain->SetBranchAddress("Sisrtene", Sisrtene, &b_Sisrtene);
   fChain->SetBranchAddress("Sihespos", Sihespos, &b_Sihespos);
   fChain->SetBranchAddress("Sihesene", Sihesene, &b_Sihesene);
   fChain->SetBranchAddress("Sihesr", Sihesr, &b_Sihesr);
   fChain->SetBranchAddress("Siprsene", Siprsene, &b_Siprsene);
   fChain->SetBranchAddress("Siccet", Siccet, &b_Siccet);
   fChain->SetBranchAddress("Siccempz", Siccempz, &b_Siccempz);
   fChain->SetBranchAddress("Sietamax", Sietamax, &b_Sietamax);
   fChain->SetBranchAddress("Sicehmom", Sicehmom, &b_Sicehmom);
   fChain->SetBranchAddress("Sizuhmom", Sizuhmom, &b_Sizuhmom);
   fChain->SetBranchAddress("Sicchmom", Sicchmom, &b_Sicchmom);
   fChain->SetBranchAddress("Sixel", Sixel, &b_Sixel);
   fChain->SetBranchAddress("Siyel", Siyel, &b_Siyel);
   fChain->SetBranchAddress("Siq2el", Siq2el, &b_Siq2el);
   fChain->SetBranchAddress("Sixda", Sixda, &b_Sixda);
   fChain->SetBranchAddress("Siyda", Siyda, &b_Siyda);
   fChain->SetBranchAddress("Siq2da", Siq2da, &b_Siq2da);
   fChain->SetBranchAddress("Sixda_cell", Sixda_cell, &b_Sixda_cell);
   fChain->SetBranchAddress("Siyda_cell", Siyda_cell, &b_Siyda_cell);
   fChain->SetBranchAddress("Siq2da_cell", Siq2da_cell, &b_Siq2da_cell);
   fChain->SetBranchAddress("Sixjb", Sixjb, &b_Sixjb);
   fChain->SetBranchAddress("Siyjb", Siyjb, &b_Siyjb);
   fChain->SetBranchAddress("Siq2jb", Siq2jb, &b_Siq2jb);
   fChain->SetBranchAddress("Sixjb_cell", Sixjb_cell, &b_Sixjb_cell);
   fChain->SetBranchAddress("Siyjb_cell", Siyjb_cell, &b_Siyjb_cell);
   fChain->SetBranchAddress("Siq2jb_cell", Siq2jb_cell, &b_Siq2jb_cell);
   
   // THIS IS SOMETHING NEW BLOCK 
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmcsitrue", fmcsitrue, &b_fmcsitrue);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmcsinpart", fmcsinpart, &b_fmcsinpart);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmcsipartlist", fmcsipartlist, &b_fmcsipartlist);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fmcsienelist", fmcsienelist, &b_fmcsienelist);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("phi_tr", &phi_tr, &b_phi_tr);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("zvtx_fcal", &zvtx_fcal, &b_zvtx_fcal);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fcal_nrgoodcells", &fcal_nrgoodcells, &b_fcal_nrgoodcells);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("fcal_vzerr", &fcal_vzerr, &b_fcal_vzerr);
       fChain->SetBranchAddress("V_h_px_zu", &V_h_px_zu, &b_V_h_px_zu);
       fChain->SetBranchAddress("V_h_py_zu", &V_h_py_zu, &b_V_h_py_zu);
       fChain->SetBranchAddress("V_h_pz_zu", &V_h_pz_zu, &b_V_h_pz_zu);
       fChain->SetBranchAddress("V_h_e_zu", &V_h_e_zu, &b_V_h_e_zu);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("npart_my", &npart_my, &b_npart_my);// NOT IN NTUPLES08
       fChain->SetBranchAddress("Part_jetid", Part_jetid, &b_Part_jetid);
       fChain->SetBranchAddress("Part_isthep", Part_isthep, &b_Part_isthep);
       fChain->SetBranchAddress("Part_charge", Part_charge, &b_Part_charge);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photn", &Photn, &b_photn);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photprphjetindex", &photprphjetindex, &b_photprphjetindex);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photid", Photid, &b_photid);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photp", Photp, &b_photp);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photmothprt", photmothprt, &b_photmothprt);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photmothid", photmothid, &b_photmothid);
       //Error in <TTree::SetBranchStatus>: unknown branch -> fChain->SetBranchAddress("photisthep", photisthep, &b_photisthep);

       // NOT IN NTUPLES08
       fChain->SetBranchAddress("Ktrnjets", &Ktrnjets, &b_Ktrnjets);
       fChain->SetBranchAddress("Ktrjets", Ktrjets, &b_Ktrjets);
       fChain->SetBranchAddress("Ktrjetet", Ktrjetet, &b_Ktrjetet);
       fChain->SetBranchAddress("Ktrjetpt", Ktrjetpt, &b_Ktrjetpt);
       fChain->SetBranchAddress("Ktrjeteta", Ktrjeteta, &b_Ktrjeteta);
       fChain->SetBranchAddress("Ktrjetphi", Ktrjetphi, &b_Ktrjetphi);
       fChain->SetBranchAddress("Ktrjetntrac", Ktrjetntrac, &b_Ktrjetntrac);
       fChain->SetBranchAddress("Ktrjetchar", Ktrjetchar, &b_Ktrjetchar);
       //---------------------------
  //---------------------------




   fChain->SetBranchAddress("Nbpchn", &Nbpchn, &b_nbpchn);
   fChain->SetBranchAddress("Bpmip", Bpmip, &b_Bpmip);
   fChain->SetBranchAddress("Bpxyz", Bpxyz, &b_Bpxyz);
   fChain->SetBranchAddress("Bpchn", Bpchn, &b_Bpchn);
   fChain->SetBranchAddress("Bptim", Bptim, &b_Bptim);
   fChain->SetBranchAddress("Ntrkvtx", &Ntrkvtx, &b_Ntrkvtx);
   fChain->SetBranchAddress("Xvtx", &Xvtx, &b_Xvtx);
   fChain->SetBranchAddress("Yvtx", &Yvtx, &b_Yvtx);
   fChain->SetBranchAddress("Zvtx", &Zvtx, &b_Zvtx);
   fChain->SetBranchAddress("Chivtx", &Chivtx, &b_Chivtx);
   fChain->SetBranchAddress("Nsecvtx", &Nsecvtx, &b_Nsecvtx);
   fChain->SetBranchAddress("Xsecvtx", Xsecvtx, &b_Xsecvtx);
   fChain->SetBranchAddress("Ysecvtx", Ysecvtx, &b_Ysecvtx);
   fChain->SetBranchAddress("Zsecvtx", Zsecvtx, &b_Zsecvtx);
   fChain->SetBranchAddress("Chisecvtx", Chisecvtx, &b_Chisecvtx);
   fChain->SetBranchAddress("Fetatr", &Fetatr, &b_Fetatr);
   fChain->SetBranchAddress("Betatr", &Betatr, &b_Betatr);
   fChain->SetBranchAddress("Pt_tr", &Pt_tr, &b_Pt_tr);
   fChain->SetBranchAddress("Empz_tr_pi", &Empz_tr_pi, &b_Empz_tr_pi);
   fChain->SetBranchAddress("Et_tr", &Et_tr, &b_Et_tr);
   fChain->SetBranchAddress("E_tr_pi", &E_tr_pi, &b_E_tr_pi);
   
   fChain->SetBranchAddress("Etamax_zu", &Etamax_zu, &b_Etamax_zu);
   fChain->SetBranchAddress("Etamax_zu4", &Etamax_zu4, &b_Etamax_zu4);
   fChain->SetBranchAddress("Fgap", &Fgap, &b_Fgap);
   fChain->SetBranchAddress("Bgap", &Bgap, &b_Bgap);
   fChain->SetBranchAddress("Nzufos", &Nzufos, &b_Nzufos);
  
      fChain->SetBranchAddress("Tufo", Tufo, &b_tufo);
      fChain->SetBranchAddress("Zufo_bsp", zufo_bsp, &b_zufo_bsp);
      fChain->SetBranchAddress("Zufo", Zufo, &b_zufo);
      fChain->SetBranchAddress("Cufo", cufo, &b_cufo);
     fChain->SetBranchAddress("Zufoecal", Zufoecal, &b_zufoecal);
     fChain->SetBranchAddress("Zufoeemc", Zufoeemc, &b_zufoeemc);
     fChain->SetBranchAddress("Zufopos", zufopos, &b_zufopos);
     fChain->SetBranchAddress("Trk_ntracks", &Trk_ntracks, &b_trk_ntracks);
     fChain->SetBranchAddress("Trk_type", trk_type, &b_trk_type);
     fChain->SetBranchAddress("Ntrack_type", ntrack_type, &b_ntrack_type);
     fChain->SetBranchAddress("Def_trk_type", &def_trk_type, &b_def_trk_type);
     fChain->SetBranchAddress("Trk_id", trk_id, &b_trk_id);
     fChain->SetBranchAddress("Trk_id2", trk_id2, &b_trk_id2);
     fChain->SetBranchAddress("Trk_px", Trk_px, &b_trk_px);
     fChain->SetBranchAddress("Trk_py", Trk_py, &b_trk_py);
     fChain->SetBranchAddress("Trk_pz", Trk_pz, &b_trk_pz);
     fChain->SetBranchAddress("Trk_charge", trk_charge, &b_trk_charge);
     fChain->SetBranchAddress("Trk_vtx", trk_vtx, &b_trk_vtx);
     fChain->SetBranchAddress("Trk_prim_vtx", Trk_prim_vtx, &b_Trk_prim_vtx);
     fChain->SetBranchAddress("Trk_sec_vtx", trk_sec_vtx, &b_trk_sec_vtx);
     fChain->SetBranchAddress("Trk_vxid", trk_vxid, &b_trk_vxid);
     fChain->SetBranchAddress("Trk_endpos", trk_endpos, &b_trk_endpos);
     fChain->SetBranchAddress("Trk_nzbyt", trk_nzbyt, &b_trk_nzbyt);
     fChain->SetBranchAddress("Trk_naxial", trk_naxial, &b_trk_naxial);
     fChain->SetBranchAddress("Trk_nstereo", trk_nstereo, &b_trk_nstereo);
     fChain->SetBranchAddress("Trk_ndof", trk_ndof, &b_trk_ndof);
     fChain->SetBranchAddress("Trk_layinner", trk_layinner, &b_trk_layinner);
     fChain->SetBranchAddress("Trk_layouter", trk_layouter, &b_trk_layouter);
     fChain->SetBranchAddress("Trk_dedxctd", trk_dedxctd, &b_trk_dedxctd);
     fChain->SetBranchAddress("Trk_dedxctdcr", trk_dedxctdcr, &b_trk_dedxctdcr);
     fChain->SetBranchAddress("Trk_dedxctderr", trk_dedxctderr, &b_trk_dedxctderr);
     fChain->SetBranchAddress("Trk_dedxctdnh", trk_dedxctdnh, &b_trk_dedxctdnh);
     fChain->SetBranchAddress("Trk_dedxctdsaturated", trk_dedxctdsaturated, &b_trk_dedxctdsaturated);
     fChain->SetBranchAddress("Trk_chi2", trk_chi2, &b_trk_chi2);
     fChain->SetBranchAddress("Trk_vchi2", trk_vchi2, &b_trk_vchi2);
     fChain->SetBranchAddress("Trk_imppar", trk_imppar, &b_trk_imppar);
     fChain->SetBranchAddress("Trk_imperr", trk_imperr, &b_trk_imperr);
     fChain->SetBranchAddress("Trk_pca", trk_pca, &b_trk_pca);
     fChain->SetBranchAddress("Filter", &filter, &b_filter);
   
   fChain->SetBranchAddress("Knjets", &Knjets, &b_Knjets);
   fChain->SetBranchAddress("Kpjets", Kpjets, &b_Kpjets);
   fChain->SetBranchAddress("Kpjetet", Kpjetet, &b_Kpjetet);
   fChain->SetBranchAddress("Kpjetpt", Kpjetpt, &b_Kpjetpt);
   fChain->SetBranchAddress("Kpjeteta", Kpjeteta, &b_Kpjeteta);
   fChain->SetBranchAddress("Kpjetphi", Kpjetphi, &b_Kpjetphi);
   fChain->SetBranchAddress("Kpjetnzu", Kpjetnzu, &b_Kpjetnzu);
   fChain->SetBranchAddress("Kpjetemcfrac", Kpjetemcfrac, &b_Kpjetemcfrac);
   fChain->SetBranchAddress("Kpjetnisl", Kpjetnisl, &b_Kpjetnisl);
   fChain->SetBranchAddress("Kpjetfmax", Kpjetfmax, &b_Kpjetfmax);
   fChain->SetBranchAddress("Kpjetdeltaz", Kpjetdeltaz, &b_Kpjetdeltaz);

   fChain->SetBranchAddress("Knzufos", &Knzufos, &b_Knzufos);
   fChain->SetBranchAddress("Kzufos", Kzufos, &b_Kzufos);
   fChain->SetBranchAddress("Kzufoet", Kzufoet, &b_Kzufoet);
   fChain->SetBranchAddress("Kzufopt", Kzufopt, &b_Kzufopt);
   fChain->SetBranchAddress("Kzufoeta", Kzufoeta, &b_Kzufoeta);
   fChain->SetBranchAddress("Kzufophi", Kzufophi, &b_Kzufophi);
   fChain->SetBranchAddress("Kzufoemcfrac", Kzufoemcfrac, &b_Kzufoemcfrac);
   fChain->SetBranchAddress("Kzufofmax", Kzufofmax, &b_Kzufofmax);
   fChain->SetBranchAddress("Kzufodeltaz", Kzufodeltaz, &b_Kzufodeltaz);
   fChain->SetBranchAddress("Kzufotype", Kzufotype, &b_Kzufotype);
   fChain->SetBranchAddress("Kzufoidjet", Kzufoidjet, &b_Kzufoidjet);
   fChain->SetBranchAddress("Kzufoncells", Kzufoncells, &b_Kzufoncells);
      //Block: QCDBOSON
        fChain->SetBranchAddress("Bospx", &bospx, &b_bospx);
        fChain->SetBranchAddress("Bospy", &bospy, &b_bospy);
        fChain->SetBranchAddress("Bospz", &bospz, &b_bospz);
        fChain->SetBranchAddress("Bosene", &bosene, &b_bosene);
      //Block: QCDPAR
       fChain->SetBranchAddress("Nppart", &Nppart, &b_Nppart);
       fChain->SetBranchAddress("Idpart", Idpart, &b_Idpart);
       fChain->SetBranchAddress("Ppart", Ppart, &b_Ppart);
      //Block: QCDHAD
        fChain->SetBranchAddress("Nfmckin", &Nfmckin, &b_Nfmckin);
        fChain->SetBranchAddress("Idfmckin", Idfmckin, &b_Idfmckin);
        fChain->SetBranchAddress("Ppfmckin", Ppfmckin, &b_Ppfmckin);
      //Block: FMCKin
        fChain->SetBranchAddress("Npart", &Npart, &b_npart);
        fChain->SetBranchAddress("Idlepton", &idlepton, &b_idlepton);
        fChain->SetBranchAddress("Idphoton", &idphoton, &b_idphoton);
        fChain->SetBranchAddress("Part_id", Part_id, &b_Part_id);
        fChain->SetBranchAddress("Part_motherid", Part_motherid, &b_Part_motherid);
        fChain->SetBranchAddress("Part_motherprt", Part_motherprt, &b_Part_motherprt);
        fChain->SetBranchAddress("Part_prt", Part_prt, &b_Part_prt);
        fChain->SetBranchAddress("Part_p", Part_p, &b_Part_p);
      //Block: FMCKIN2
        fChain->SetBranchAddress("Fmck_nstor", &Fmck_nstor, &b_Fmck_nstor);
        fChain->SetBranchAddress("Fmck_id", &Fmck_id, &b_Fmck_id);
        fChain->SetBranchAddress("Fmck_prt", &Fmck_prt, &b_Fmck_prt);
        fChain->SetBranchAddress("Fmck_px", &Fmck_px, &b_Fmck_px);
        fChain->SetBranchAddress("Fmck_py", &Fmck_py, &b_Fmck_py);
        fChain->SetBranchAddress("Fmck_pz", &Fmck_pz, &b_Fmck_pz);
        fChain->SetBranchAddress("Fmck_e", &Fmck_e, &b_Fmck_e);
        fChain->SetBranchAddress("Fmck_m", &Fmck_m, &b_Fmck_m);
        fChain->SetBranchAddress("Fmck_isthep", &Fmck_isthep, &b_Fmck_isthep);
        fChain->SetBranchAddress("Fmck_daug", &Fmck_daug, &b_Fmck_daug);


  // THIS IS SOMETHING NEW BLOCK 
   /*Error in <TTree::SetBranchStatus>: unknown branch -> 
     fChain->SetBranchAddress("dattyp", &dattyp, &b_dattyp);
     fChain->SetBranchAddress("bit3_tlt4", &bit3_tlt4, &b_bit3_tlt4);
     fChain->SetBranchAddress("tlt4", tlt4, &b_tlt4);
     fChain->SetBranchAddress("Siein03", &Siein03, &b_Siein03);
     fChain->SetBranchAddress("Sienin03", &Sienin03, &b_Sienin03);
     fChain->SetBranchAddress("mc_fla", mc_fla, &b_mc_fla);
     fChain->SetBranchAddress("mc_juscal", &mc_juscal, &b_mc_juscal);
     fChain->SetBranchAddress("mc_jux1ge", &mc_jux1ge, &b_mc_jux1ge);
     fChain->SetBranchAddress("mc_jux2ge", &mc_jux2ge, &b_mc_jux2ge);
     fChain->SetBranchAddress("mc_np", &mc_np, &b_mc_np);
     fChain->SetBranchAddress("mc_pa", mc_pa, &b_mc_pa);
     fChain->SetBranchAddress("mc_eparton", mc_eparton, &b_mc_eparton);
     fChain->SetBranchAddress("mc_efla", &mc_efla, &b_mc_efla);
     fChain->SetBranchAddress("mc_ygen", &mc_ygen, &b_mc_ygen);
     fChain->SetBranchAddress("mc_xgen", &mc_xgen, &b_mc_xgen);
     fChain->SetBranchAddress("mc_q2gen", &mc_q2gen, &b_mc_q2gen);
     fChain->SetBranchAddress("mc_gaparton", mc_gaparton, &b_mc_gaparton);
     fChain->SetBranchAddress("mc_yapp", &mc_yapp, &b_mc_yapp);
     fChain->SetBranchAddress("mc_xapp", &mc_xapp, &b_mc_xapp);
     fChain->SetBranchAddress("mc_q2app", &mc_q2app, &b_mc_q2app);
     fChain->SetBranchAddress("mc_xxt", &mc_xxt, &b_mc_xxt);
     fChain->SetBranchAddress("mc_yyt", &mc_yyt, &b_mc_yyt);
     fChain->SetBranchAddress("mc_zzt", &mc_zzt, &b_mc_zzt);
     fChain->SetBranchAddress("mc_gparton", mc_gparton, &b_mc_gparton);
     fChain->SetBranchAddress("mc_pparton", mc_pparton, &b_mc_pparton);
     fChain->SetBranchAddress("Mcebeam", &Mcebeam, &b_Mcebeam);
     fChain->SetBranchAddress("Mcpbeam", &Mcpbeam, &b_Mcpbeam);
     fChain->SetBranchAddress("Mcelec", &Mcelec, &b_Mcelec);
     fChain->SetBranchAddress("Mcethe", &Mcethe, &b_Mcethe);
     fChain->SetBranchAddress("Mcephi", &Mcephi, &b_Mcephi);
     fChain->SetBranchAddress("Mcq2", &Mcq2, &b_Mcq2);
     fChain->SetBranchAddress("Mcxbj", &Mcxbj, &b_Mcxbj);
     fChain->SetBranchAddress("Mcybj", &Mcybj, &b_Mcybj);
     fChain->SetBranchAddress("Mcgam", Mcgam, &b_Mcgam);
     fChain->SetBranchAddress("vtxpos", vtxpos, &b_vtxpos);
   */
   //---------------------------
}

Bool_t selector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef selector_cxx
//#endif


