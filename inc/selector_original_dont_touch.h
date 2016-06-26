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

class selector {
 public :
  Double_t jet_en_corr(Double_t eta, Double_t et, TString period, TString mc_type);
  Double_t jet_en_corr(Double_t eta, Double_t et, Double_t* C_scale, Double_t* a0_corr, Double_t* a1_corr, Double_t* a2_corr);
  Double_t BPRES_mips(TVector3 v_photon, Double_t radius);
  void  getBpreReadout(); // A.Volynets
  Double_t doBpreReadout(Double_t mips); // A.Volynets
  void getRunsBPRE(); // A.Volynets
  Double_t Land[100];          // BPRE smearing constants
  map<int,int>   BadBpre;
  void fill_trigger_bits(TH1D* h[], Int_t trigger[]);
  void fill_trigger_bits_general(TH1D* h[], Int_t trigger[], Int_t n);
  Bool_t SelectHadronLevel(Bool_t take_det_event);
  Bool_t SelectPrPhKtJetB(Int_t i_true_rad_photon);
  Double_t q2_reweighting(Float_t Mc_q2, TString mc_type, TString period);
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
  Double_t        jet_energy_frac;
  Double_t        cell_energy_frac;
  Double_t        photon_corr_et;
  Double_t        max_et_candidate;
  Double_t        max_et_candidate_uncorr;
  Double_t        glob_fmax;
  Double_t        glob_deltaz;
  Int_t           det_cross_sec_bin_et;
  Int_t           det_cross_sec_bin_eta;
  Int_t           det_cross_sec_bin_q2;
  Int_t           det_cross_sec_bin_x;
  Int_t           det_cross_sec_bin_et_jet;
  Int_t           det_cross_sec_bin_eta_jet;
  TString         period;
  TString         mc_type;
  TString         mc_corr_type;
  TString         root_file_name;
  TString         treeName;
  TFile           *f;
  TTree           *selectedTree;
  TDatime         time_begin;
  TDatime         time_end;
  TRandom3        *random;
  TLorentzVector  *v_uncorr_prompt_photon;
  TLorentzVector  *v_corr_prompt_photon;
  TLorentzVector  *v_prompt_photon_jet;
  TLorentzVector  *v_accomp_corr_jet;
  TLorentzVector  *v_accomp_uncorr_jet;
  Int_t           radiated_candidate_number;
  Int_t           max_et_candidate_number;
  Int_t           candidate_jet_number;
  //  Int_t           
  static const Float_t m_Lepton_Energy = 27.52;
  static const Float_t m_Proton_Energy = 920.;
  static const  Double_t E_e = 27.52;
  static const  Double_t E_p = 920.;
  static const Int_t maxNjetsInEvent = 50;
  static const Int_t maxNofJets = 200;
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
  /*
    TH1D* h_cutcount_q2 =  new TH1D("h_cutcount_q2","h_cutcount_q2", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_y_da =  new TH1D("h_cutcount_y_da","h_cutcount_y_da", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_zvtx =  new TH1D("h_cutcount_zvtx","h_cutcount_zvtx", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_sincand =  new TH1D("h_cutcount_sincand","h_cutcount_sincand", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_el_energie =  new TH1D("cutcount_el_energie","cutcount_el_energie", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_distance_to_beam =  new TH1D("distance_to_beam","distance_to_beam", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_enin =  new TH1D("h_cutcount_enin","h_cutcount_enin", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_empz =  new TH1D("h_cutcount_empz","h_cutcount_empz", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_pt_et =  new TH1D("h_cutcount_pt_et","h_cutcount_pt_et", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_chimney =  new TH1D("h_cutcount_chimney","h_cutcount_chimney", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_cutcount_super_crack =  new TH1D("h_cutcount_super_crack","h_cutcount_super_crack", 31879,30758.0, (62636.0 + 1.0));
    TH1D* h_jetsB5_eta_lab = new TH1D("h_jetsB5_eta_lab", "h_jetsB5_eta_lab", 20, -1., 2.5);
  */
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


   // Declaration of leaf types
   Int_t           Runnr;
   Int_t           Eventnr;
   Int_t           Evtake;
   Int_t           Evtake_iwant;
   Int_t           Mbtake;
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
   Int_t           Cal_tf_n;
   Int_t           Cal_tb_n;
   Int_t           Cal_tr_n;
   Int_t           Cal_tg_n;
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
   Int_t           E5ncand;
   Int_t           E5error;
   Float_t         E5prob[4];   //[E5Ncand]
   Float_t         E5pos[4][3];   //[E5Ncand]
   Float_t         E5calpos[4][3];   //[E5Ncand]
   Float_t         E5calene[4];   //[E5Ncand]
   Float_t         E5ein[4];   //[E5Ncand]
   Float_t         E5enin[4];   //[E5Ncand]
   Float_t         E5ecorr[4][3];   //[E5Ncand]
   Float_t         E5th[4];   //[E5Ncand]
   Float_t         E5ph[4];   //[E5Ncand]
   Float_t         E5pt[4];   //[E5Ncand]
   Int_t           E5xdet[4][3];   //[E5Ncand]
   Int_t           E5ydet[4][3];   //[E5Ncand]
   Int_t           E5trknr[4];   //[E5Ncand]
   Int_t           E5nrsl[4];   //[E5Ncand]
   Float_t         E5dca[4];   //[E5Ncand]
   Float_t         E5dcabeam[4];   //[E5Ncand]
   Float_t         E5trkp[4];   //[E5Ncand]
   Float_t         E5trkth[4];   //[E5Ncand]
   Float_t         E5trkph[4];   //[E5Ncand]
   Float_t         E5trkq[4];   //[E5Ncand]
   Float_t         E5trkdme[4];   //[E5Ncand]
   Float_t         E5trkdce[4];   //[E5Ncand]
   Float_t         E5trkpos[4][3];   //[E5Ncand]
   Float_t         E5xel[4];   //[E5Ncand]
   Float_t         E5yel[4];   //[E5Ncand]
   Float_t         E5q2el[4];   //[E5Ncand]
   Float_t         E5xda[4];   //[E5Ncand]
   Float_t         E5yda[4];   //[E5Ncand]
   Float_t         E5q2da[4];   //[E5Ncand]
   Float_t         E5xda_cell[4];   //[E5Ncand]
   Float_t         E5yda_cell[4];   //[E5Ncand]
   Float_t         E5q2da_cell[4];   //[E5Ncand]
   Float_t         E5xjb[4];   //[E5Ncand]
   Float_t         E5yjb[4];   //[E5Ncand]
   Float_t         E5q2jb[4];   //[E5Ncand]
   Float_t         E5xjb_cell[4];   //[E5Ncand]
   Float_t         E5yjb_cell[4];   //[E5Ncand]
   Float_t         E5q2jb_cell[4];   //[E5Ncand]
   Int_t           E5ncell[4];   //[E5Ncand]
   Int_t           E5cellptr[4];   //[E5Ncand]
   Float_t         E5femc[4];   //[E5Ncand]
   Float_t         E5fcell[4][7];   //[E5Ncand]
   Float_t         E5fmodu[4][5];   //[E5Ncand]
   Float_t         E5imbmod[4][5];   //[E5Ncand]
   Float_t         E5subprob[4][7];   //[E5Ncand]
   Float_t         E5dphi[4];   //[E5Ncand]
   Float_t         E5dtheta[4];   //[E5Ncand]
   Float_t         E5calprob[4];   //[E5Ncand]
   Int_t           E5calprobrank[4];   //[E5Ncand]
   Float_t         E5fmaxbemc[4];   //[E5Ncand]
   Float_t         E5fmaxremc[4];   //[E5Ncand]
   Float_t         E5fmaxfemc[4];   //[E5Ncand]
   Float_t         E5deltaz[4];   //[E5Ncand]
   Float_t         E5deltax[4];   //[E5Ncand]
   Float_t         E5deltay[4];   //[E5Ncand]
   Int_t           Sincand;
   Int_t           Sierror;
   Float_t         Siprob[4];   //[SiNcand]
   Float_t         Sipos[4][3];   //[SiNcand]
   Float_t         Sicalpos[4][3];   //[SiNcand]
   Float_t         Sicalene[4];   //[SiNcand]
   Float_t         Siein[4];   //[SiNcand]
   Float_t         Sienin[4];   //[SiNcand]
   Float_t         Siecorr[4][3];   //[SiNcand]
   Float_t         Sith[4];   //[SiNcand]
   Float_t         Siph[4];   //[SiNcand]
   Float_t         Sipt[4];   //[SiNcand]
   Int_t           Sixdet[4][3];   //[SiNcand]
   Int_t           Siydet[4][3];   //[SiNcand]
   Int_t           Sitrknr[4];   //[SiNcand]
   Int_t           Sinrsl[4];   //[SiNcand]
   Float_t         Sidca[4];   //[SiNcand]
   Float_t         Sitrkp[4];   //[SiNcand]
   Float_t         Sitrkth[4];   //[SiNcand]
   Float_t         Sitrkph[4];   //[SiNcand]
   Float_t         Sitrkq[4];   //[SiNcand]
   Float_t         Sitrkdme[4];   //[SiNcand]
   Float_t         Sitrkpos[4][3];   //[SiNcand]
   Int_t           Sisrtf[4];   //[SiNcand]
   Int_t           Sisrtquad[4];   //[SiNcand]
   Int_t           Sihesf[4];   //[SiNcand]
   Int_t           Sicorrcode[4];   //[SiNcand]
   Float_t         Sisrtpos[4][2];   //[SiNcand]
   Float_t         Sisrtene[4];   //[SiNcand]
   Float_t         Sihespos[4][2];   //[SiNcand]
   Float_t         Sihesene[4];   //[SiNcand]
   Float_t         Sihesr[4];   //[SiNcand]
   Float_t         Siprsene[4][3];   //[SiNcand]
   Float_t         Siccet[4];   //[SiNcand]
   Float_t         Siccempz[4];   //[SiNcand]
   Float_t         Sietamax[4];   //[SiNcand]
   Float_t         Sicehmom[4][4];   //[SiNcand]
   Float_t         Sizuhmom[4][4];   //[SiNcand]
   Float_t         Sicchmom[4][4];   //[SiNcand]
   Float_t         Sixel[4];   //[SiNcand]
   Float_t         Siyel[4];   //[SiNcand]
   Float_t         Siq2el[4];   //[SiNcand]
   Float_t         Sixda[4];   //[SiNcand]
   Float_t         Siyda[4];   //[SiNcand]
   Float_t         Siq2da[4];   //[SiNcand]
   Float_t         Sixda_cell[4];   //[SiNcand]
   Float_t         Siyda_cell[4];   //[SiNcand]
   Float_t         Siq2da_cell[4];   //[SiNcand]
   Float_t         Sixjb[4];   //[SiNcand]
   Float_t         Siyjb[4];   //[SiNcand]
   Float_t         Siq2jb[4];   //[SiNcand]
   Float_t         Sixjb_cell[4];   //[SiNcand]
   Float_t         Siyjb_cell[4];   //[SiNcand]
   Float_t         Siq2jb_cell[4];   //[SiNcand]
   Int_t           Nbpchn;
   Float_t         Bpmip[432];   //[nBPchn]
   Float_t         Bpxyz[432][3];   //[nBPchn]
   Int_t           Bpchn[432];   //[nBPchn]
   Float_t         Bptim[432];   //[nBPchn]
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
   Float_t         Phi_tr;
   Float_t         Zvtx_fcal;
   Int_t           Fcal_nrgoodcells;
   Float_t         Fcal_vzerr;
   Float_t         V_h_px_zu;
   Float_t         V_h_py_zu;
   Float_t         V_h_pz_zu;
   Float_t         V_h_e_zu;
   Float_t         Etamax_zu;
   Float_t         Etamax_zu4;
   Float_t         Fgap;
   Float_t         Bgap;
   Int_t           Nzufos;
   Int_t           Tufo[750][4];   //[Nzufos]
   Int_t           Zufo_bsp[750];   //[Nzufos]
   Float_t         Zufo[750][4];   //[Nzufos]
   Float_t         Cufo[750];   //[Nzufos]
   Float_t         Zufoecal[750];   //[Nzufos]
   Float_t         Zufoeemc[750];   //[Nzufos]
   Float_t         Zufopos[750][3];   //[Nzufos]
   Int_t           Trk_ntracks;
   Int_t           Trk_type[1200];   //[trk_ntracks]
   Int_t           Ntrack_type[4];
   Int_t           Def_trk_type;
   Int_t           Trk_id[1200];   //[trk_ntracks]
   Int_t           Trk_id2[1200];   //[trk_ntracks]
   Float_t         Trk_px[1200];   //[trk_ntracks]
   Float_t         Trk_py[1200];   //[trk_ntracks]
   Float_t         Trk_pz[1200];   //[trk_ntracks]
   Float_t         Trk_charge[1200];   //[trk_ntracks]
   Int_t           Trk_vtx[1200];   //[trk_ntracks]
   Int_t           Trk_prim_vtx[1200];   //[trk_ntracks]
   Int_t           Trk_sec_vtx[1200];   //[trk_ntracks]
   Int_t           Trk_vxid[1200];   //[trk_ntracks]
   Float_t         Trk_endpos[1200][3];   //[trk_ntracks]
   Int_t           Trk_nzbyt[1200];   //[trk_ntracks]
   Int_t           Trk_naxial[1200];   //[trk_ntracks]
   Int_t           Trk_nstereo[1200];   //[trk_ntracks]
   Int_t           Trk_ndof[1200];   //[trk_ntracks]
   Int_t           Trk_layinner[1200];   //[trk_ntracks]
   Int_t           Trk_layouter[1200];   //[trk_ntracks]
   Float_t         Trk_dedxctd[1200];   //[trk_ntracks]
   Float_t         Trk_dedxctdcr[1200];   //[trk_ntracks]
   Int_t           Trk_dedxctderr[1200];   //[trk_ntracks]
   Int_t           Trk_dedxctdnh[1200];   //[trk_ntracks]
   Int_t           Trk_dedxctdsaturated[1200];   //[trk_ntracks]
   Float_t         Trk_chi2[1200];   //[trk_ntracks]
   Float_t         Trk_vchi2[1200];   //[trk_ntracks]
   Float_t         Trk_imppar[1200];   //[trk_ntracks]
   Float_t         Trk_imperr[1200];   //[trk_ntracks]
   Float_t         Trk_pca[1200][3];   //[trk_ntracks]
   Int_t           Filter;
   Int_t           Knjets;
   Float_t         Kpjets[50][4];   //[KNjets]
   Float_t         Kpjetet[50];   //[KNjets]
   Float_t         Kpjetpt[50];   //[KNjets]
   Float_t         Kpjeteta[50];   //[KNjets]
   Float_t         Kpjetphi[50];   //[KNjets]
   Int_t           Kpjetnzu[50];   //[KNjets]
   Float_t         Kpjetemcfrac[50];   //[KNjets]
   Int_t           Kpjetnisl[50];   //[KNjets]
   Float_t         Kpjetfmax[50];   //[KNjets]
   Float_t         Kpjetdeltaz[50];   //[KNjets]
   Int_t           Knzufos;
   Float_t         Kzufos[250][4];   //[KNzufos]
   Float_t         Kzufoet[250];   //[KNzufos]
   Float_t         Kzufopt[250];   //[KNzufos]
   Float_t         Kzufoeta[250];   //[KNzufos]
   Float_t         Kzufophi[250];   //[KNzufos]
   Float_t         Kzufoemcfrac[250];   //[KNzufos]
   Float_t         Kzufofmax[250];   //[KNzufos]
   Float_t         Kzufodeltaz[250];   //[KNzufos]
   Float_t         myKzufodeltaz[250]; //to divide by 5.45
   Int_t           Kzufotype[250];   //[KNzufos]
   Int_t           Kzufoidjet[250];   //[KNzufos]
   Int_t           Kzufoncells[250];   //[KNzufos]
   Float_t         Mc_q2;
   Int_t           Photn;
   Int_t           Photid[100];   //[photN]
   Float_t         Photp[100][4];   //[photN]
   Int_t           Npart;
   Int_t           Part_charge[511];
   Int_t           Part_prt[511];   //[nPart]
   Int_t           Part_isthep[511];   //[nPart]
   Int_t           Part_id[511];   //[nPart]
   Float_t         Part_p[511][4];   //[nPart]
   Float_t         Mc_pfsl[4];
   Int_t           Part_jetid[511];   //[nPart_my]
   Int_t           Ktrnjets;
   Float_t         Ktrjets[50][4];   //[KtrNJets]
   Float_t         Mc_x;
   Float_t         Mc_y;


   // List of branches
   TBranch        *b_RUNNR;   //!
   TBranch        *b_EVENTNR;   //!
   TBranch        *b_EVTAKE;   //!
   TBranch        *b_EVTAKE_iwant;   //!
   TBranch        *b_MBTAKE;   //!
   TBranch        *b_FLTW;   //!
   TBranch        *b_FLTPSW;   //!
   TBranch        *b_FLTfl;   //!
   TBranch        *b_GSLT_global;   //!
   TBranch        *b_SLTW;   //!
   TBranch        *b_SLTUPW;   //!
   TBranch        *b_TLTW;   //!
   TBranch        *b_DSTW;   //!
   TBranch        *b_FLTPSFCW;   //!
   TBranch        *b_CAL_px;   //!
   TBranch        *b_CAL_py;   //!
   TBranch        *b_CAL_pz;   //!
   TBranch        *b_CAL_E;   //!
   TBranch        *b_CAL_Et;   //!
   TBranch        *b_CAL_Empz;   //!
   TBranch        *b_CAL_Pt;   //!
   TBranch        *b_CAL_phi;   //!
   TBranch        *b_REmc;   //!
   TBranch        *b_BEmc;   //!
   TBranch        *b_FEmc;   //!
   TBranch        *b_RHac;   //!
   TBranch        *b_BHac;   //!
   TBranch        *b_FHac;   //!
   TBranch        *b_BHac2;   //!
   TBranch        *b_FHac2;   //!
   TBranch        *b_NFEMC;   //!
   TBranch        *b_NFHAC1;   //!
   TBranch        *b_NFHAC2;   //!
   TBranch        *b_NBEMC;   //!
   TBranch        *b_NBHAC1;   //!
   TBranch        *b_NBHAC2;   //!
   TBranch        *b_NREMC;   //!
   TBranch        *b_NRHAC;   //!
   TBranch        *b_CAL_Tf;   //!
   TBranch        *b_CAL_Tb;   //!
   TBranch        *b_CAL_Tr;   //!
   TBranch        *b_CAL_Tg;   //!
   TBranch        *b_CAL_Tu;   //!
   TBranch        *b_CAL_Td;   //!
   TBranch        *b_CAL_Tf_E;   //!
   TBranch        *b_CAL_Tb_E;   //!
   TBranch        *b_CAL_Tr_E;   //!
   TBranch        *b_CAL_Tg_E;   //!
   TBranch        *b_CAL_Tu_E;   //!
   TBranch        *b_CAL_Td_E;   //!
   TBranch        *b_CAL_Tf_n;   //!
   TBranch        *b_CAL_Tb_n;   //!
   TBranch        *b_CAL_Tr_n;   //!
   TBranch        *b_CAL_Tg_n;   //!
   TBranch        *b_CAL_Tu_n;   //!
   TBranch        *b_CAL_Td_n;   //!
   TBranch        *b_Etamax_ce;   //!
   TBranch        *b_Etamax_ce4;   //!
   TBranch        *b_CAL_Et10;   //!
   TBranch        *b_Mtrknoe_pi;   //!
   TBranch        *b_Mtrknoe_K;   //!
   TBranch        *b_E_HK;   //!
   TBranch        *b_UnmEn_pi;   //!
   TBranch        *b_UnmEn_k;   //!
   TBranch        *b_SparkF;   //!
   TBranch        *b_E5Ncand;   //!
   TBranch        *b_E5Error;   //!
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
   TBranch        *b_SiNcand;   //!
   TBranch        *b_SiError;   //!
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
   TBranch        *b_nBPchn;   //!
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
   TBranch        *b_fcal_NrGoodCells;   //!
   TBranch        *b_fcal_VzErr;   //!
   TBranch        *b_V_H_px_zu;   //!
   TBranch        *b_V_H_py_zu;   //!
   TBranch        *b_V_H_pz_zu;   //!
   TBranch        *b_V_H_E_zu;   //!
   TBranch        *b_Etamax_zu;   //!
   TBranch        *b_Etamax_zu4;   //!
   TBranch        *b_Fgap;   //!
   TBranch        *b_Bgap;   //!
   TBranch        *b_Nzufos;   //!
   TBranch        *b_Tufo;   //!
   TBranch        *b_Zufo_bsp;   //!
   TBranch        *b_Zufo;   //!
   TBranch        *b_Cufo;   //!
   TBranch        *b_Zufoecal;   //!
   TBranch        *b_Zufoeemc;   //!
   TBranch        *b_Zufopos;   //!
   TBranch        *b_trk_ntracks;   //!
   TBranch        *b_Trk_type;   //!
   TBranch        *b_ntrack_type;   //!
   TBranch        *b_def_trk_type;   //!
   TBranch        *b_Trk_id;   //!
   TBranch        *b_Trk_id2;   //!
   TBranch        *b_Trk_px;   //!
   TBranch        *b_Trk_py;   //!
   TBranch        *b_Trk_pz;   //!
   TBranch        *b_Trk_charge;   //!
   TBranch        *b_Trk_vtx;   //!
   TBranch        *b_Trk_prim_vtx;   //!
   TBranch        *b_Trk_sec_vtx;   //!
   TBranch        *b_Trk_vxid;   //!
   TBranch        *b_Trk_endpos;   //!
   TBranch        *b_Trk_nzbyt;   //!
   TBranch        *b_Trk_naxial;   //!
   TBranch        *b_Trk_nstereo;   //!
   TBranch        *b_Trk_ndof;   //!
   TBranch        *b_Trk_layinner;   //!
   TBranch        *b_Trk_layouter;   //!
   TBranch        *b_Trk_dedxctd;   //!
   TBranch        *b_Trk_dedxctdcr;   //!
   TBranch        *b_Trk_dedxctderr;   //!
   TBranch        *b_Trk_dedxctdnh;   //!
   TBranch        *b_Trk_dedxctdsaturated;   //!
   TBranch        *b_Trk_chi2;   //!
   TBranch        *b_Trk_vchi2;   //!
   TBranch        *b_Trk_imppar;   //!
   TBranch        *b_Trk_imperr;   //!
   TBranch        *b_Trk_pca;   //!
   TBranch        *b_filter;   //!
   TBranch        *b_KNjets;   //!
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
   TBranch        *b_KNzufos;   //!
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
   TBranch        *b_photN;   //!
   TBranch        *b_Photid;   //!
   TBranch        *b_Photp;   //!
   TBranch        *b_nPart;   //!
   TBranch        *b_Part_prt;   //!
   TBranch        *b_Part_isthep;   //!
   TBranch        *b_Part_charge; //!
   TBranch        *b_Part_id;   //!
   TBranch        *b_Part_p;   //!
   TBranch        *b_MC_Q2;   //!
   TBranch        *b_MC_PFSL;   //!
   TBranch        *b_Part_jetid;   //!
   TBranch        *b_KtrNJets;   //!
   TBranch        *b_Ktrjets;   //!
   TBranch        *b_MC_X;   //!
   TBranch        *b_MC_Y;   //!


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
   fChain->SetBranchAddress("Runnr", &Runnr, &b_RUNNR);
   fChain->SetBranchAddress("Eventnr", &Eventnr, &b_EVENTNR);
   fChain->SetBranchAddress("Evtake", &Evtake, &b_EVTAKE);
   fChain->SetBranchAddress("Evtake_iwant", &Evtake_iwant, &b_EVTAKE_iwant);
   fChain->SetBranchAddress("Mbtake", &Mbtake, &b_MBTAKE);
   fChain->SetBranchAddress("Fltw", Fltw, &b_FLTW);
   fChain->SetBranchAddress("Fltpsw", Fltpsw, &b_FLTPSW);
   fChain->SetBranchAddress("Fltfl", &Fltfl, &b_FLTfl);
   fChain->SetBranchAddress("Gslt_global", &Gslt_global, &b_GSLT_global);
   fChain->SetBranchAddress("Sltw", Sltw, &b_SLTW);
   fChain->SetBranchAddress("Sltupw", Sltupw, &b_SLTUPW);
   fChain->SetBranchAddress("Tltw", Tltw, &b_TLTW);
   fChain->SetBranchAddress("Dstw", Dstw, &b_DSTW);
   fChain->SetBranchAddress("Fltpsfcw", Fltpsfcw, &b_FLTPSFCW);
   fChain->SetBranchAddress("Cal_px", &Cal_px, &b_CAL_px);
   fChain->SetBranchAddress("Cal_py", &Cal_py, &b_CAL_py);
   fChain->SetBranchAddress("Cal_pz", &Cal_pz, &b_CAL_pz);
   fChain->SetBranchAddress("Cal_e", &Cal_e, &b_CAL_E);
   fChain->SetBranchAddress("Cal_et", &Cal_et, &b_CAL_Et);
   fChain->SetBranchAddress("Cal_empz", &Cal_empz, &b_CAL_Empz);
   fChain->SetBranchAddress("Cal_pt", &Cal_pt, &b_CAL_Pt);
   fChain->SetBranchAddress("Cal_phi", &Cal_phi, &b_CAL_phi);
   fChain->SetBranchAddress("Remc", Remc, &b_REmc);
   fChain->SetBranchAddress("Bemc", Bemc, &b_BEmc);
   fChain->SetBranchAddress("Femc", Femc, &b_FEmc);
   fChain->SetBranchAddress("Rhac", Rhac, &b_RHac);
   fChain->SetBranchAddress("Bhac", Bhac, &b_BHac);
   fChain->SetBranchAddress("Fhac", Fhac, &b_FHac);
   fChain->SetBranchAddress("Bhac2", Bhac2, &b_BHac2);
   fChain->SetBranchAddress("Fhac2", Fhac2, &b_FHac2);
   fChain->SetBranchAddress("Nfemc", &Nfemc, &b_NFEMC);
   fChain->SetBranchAddress("Nfhac1", &Nfhac1, &b_NFHAC1);
   fChain->SetBranchAddress("Nfhac2", &Nfhac2, &b_NFHAC2);
   fChain->SetBranchAddress("Nbemc", &Nbemc, &b_NBEMC);
   fChain->SetBranchAddress("Nbhac1", &Nbhac1, &b_NBHAC1);
   fChain->SetBranchAddress("Nbhac2", &Nbhac2, &b_NBHAC2);
   fChain->SetBranchAddress("Nremc", &Nremc, &b_NREMC);
   fChain->SetBranchAddress("Nrhac", &Nrhac, &b_NRHAC);
   fChain->SetBranchAddress("Cal_tf", &Cal_tf, &b_CAL_Tf);
   fChain->SetBranchAddress("Cal_tb", &Cal_tb, &b_CAL_Tb);
   fChain->SetBranchAddress("Cal_tr", &Cal_tr, &b_CAL_Tr);
   fChain->SetBranchAddress("Cal_tg", &Cal_tg, &b_CAL_Tg);
   fChain->SetBranchAddress("Cal_tu", &Cal_tu, &b_CAL_Tu);
   fChain->SetBranchAddress("Cal_td", &Cal_td, &b_CAL_Td);
   fChain->SetBranchAddress("Cal_tf_e", &Cal_tf_e, &b_CAL_Tf_E);
   fChain->SetBranchAddress("Cal_tb_e", &Cal_tb_e, &b_CAL_Tb_E);
   fChain->SetBranchAddress("Cal_tr_e", &Cal_tr_e, &b_CAL_Tr_E);
   fChain->SetBranchAddress("Cal_tg_e", &Cal_tg_e, &b_CAL_Tg_E);
   fChain->SetBranchAddress("Cal_tu_e", &Cal_tu_e, &b_CAL_Tu_E);
   fChain->SetBranchAddress("Cal_td_e", &Cal_td_e, &b_CAL_Td_E);
   fChain->SetBranchAddress("Cal_tf_n", &Cal_tf_n, &b_CAL_Tf_n);
   fChain->SetBranchAddress("Cal_tb_n", &Cal_tb_n, &b_CAL_Tb_n);
   fChain->SetBranchAddress("Cal_tr_n", &Cal_tr_n, &b_CAL_Tr_n);
   fChain->SetBranchAddress("Cal_tg_n", &Cal_tg_n, &b_CAL_Tg_n);
   fChain->SetBranchAddress("Cal_tu_n", &Cal_tu_n, &b_CAL_Tu_n);
   fChain->SetBranchAddress("Cal_td_n", &Cal_td_n, &b_CAL_Td_n);
   fChain->SetBranchAddress("Etamax_ce", &Etamax_ce, &b_Etamax_ce);
   fChain->SetBranchAddress("Etamax_ce4", &Etamax_ce4, &b_Etamax_ce4);
   fChain->SetBranchAddress("Cal_et10", &Cal_et10, &b_CAL_Et10);
   fChain->SetBranchAddress("Mtrknoe_pi", &Mtrknoe_pi, &b_Mtrknoe_pi);
   fChain->SetBranchAddress("Mtrknoe_k", &Mtrknoe_k, &b_Mtrknoe_K);
   fChain->SetBranchAddress("E_hk", &E_hk, &b_E_HK);
   fChain->SetBranchAddress("Unmen_pi", &Unmen_pi, &b_UnmEn_pi);
   fChain->SetBranchAddress("Unmen_k", &Unmen_k, &b_UnmEn_k);
   fChain->SetBranchAddress("Sparkf", &Sparkf, &b_SparkF);
   fChain->SetBranchAddress("E5ncand", &E5ncand, &b_E5Ncand);
   fChain->SetBranchAddress("E5error", &E5error, &b_E5Error);
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
   fChain->SetBranchAddress("Sincand", &Sincand, &b_SiNcand);
   fChain->SetBranchAddress("Sierror", &Sierror, &b_SiError);
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
   fChain->SetBranchAddress("Nbpchn", &Nbpchn, &b_nBPchn);
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
   fChain->SetBranchAddress("Phi_tr", &Phi_tr, &b_phi_tr);
   fChain->SetBranchAddress("Zvtx_fcal", &Zvtx_fcal, &b_zvtx_fcal);
   fChain->SetBranchAddress("Fcal_nrgoodcells", &Fcal_nrgoodcells, &b_fcal_NrGoodCells);
   fChain->SetBranchAddress("Fcal_vzerr", &Fcal_vzerr, &b_fcal_VzErr);
   fChain->SetBranchAddress("V_h_px_zu", &V_h_px_zu, &b_V_H_px_zu);
   fChain->SetBranchAddress("V_h_py_zu", &V_h_py_zu, &b_V_H_py_zu);
   fChain->SetBranchAddress("V_h_pz_zu", &V_h_pz_zu, &b_V_H_pz_zu);
   fChain->SetBranchAddress("V_h_e_zu", &V_h_e_zu, &b_V_H_E_zu);
   fChain->SetBranchAddress("Etamax_zu", &Etamax_zu, &b_Etamax_zu);
   fChain->SetBranchAddress("Etamax_zu4", &Etamax_zu4, &b_Etamax_zu4);
   fChain->SetBranchAddress("Fgap", &Fgap, &b_Fgap);
   fChain->SetBranchAddress("Bgap", &Bgap, &b_Bgap);
   fChain->SetBranchAddress("Nzufos", &Nzufos, &b_Nzufos);
   fChain->SetBranchAddress("Tufo", Tufo, &b_Tufo);
   fChain->SetBranchAddress("Zufo_bsp", Zufo_bsp, &b_Zufo_bsp);
   fChain->SetBranchAddress("Zufo", Zufo, &b_Zufo);
   fChain->SetBranchAddress("Cufo", Cufo, &b_Cufo);
   fChain->SetBranchAddress("Zufoecal", Zufoecal, &b_Zufoecal);
   fChain->SetBranchAddress("Zufoeemc", Zufoeemc, &b_Zufoeemc);
   fChain->SetBranchAddress("Zufopos", Zufopos, &b_Zufopos);
   fChain->SetBranchAddress("Trk_ntracks", &Trk_ntracks, &b_trk_ntracks);
   fChain->SetBranchAddress("Trk_type", Trk_type, &b_Trk_type);
   fChain->SetBranchAddress("Ntrack_type", Ntrack_type, &b_ntrack_type);
   fChain->SetBranchAddress("Def_trk_type", &Def_trk_type, &b_def_trk_type);
   fChain->SetBranchAddress("Trk_id", Trk_id, &b_Trk_id);
   fChain->SetBranchAddress("Trk_id2", Trk_id2, &b_Trk_id2);
   fChain->SetBranchAddress("Trk_px", Trk_px, &b_Trk_px);
   fChain->SetBranchAddress("Trk_py", Trk_py, &b_Trk_py);
   fChain->SetBranchAddress("Trk_pz", Trk_pz, &b_Trk_pz);
   fChain->SetBranchAddress("Trk_charge", Trk_charge, &b_Trk_charge);
   fChain->SetBranchAddress("Trk_vtx", Trk_vtx, &b_Trk_vtx);
   fChain->SetBranchAddress("Trk_prim_vtx", Trk_prim_vtx, &b_Trk_prim_vtx);
   fChain->SetBranchAddress("Trk_sec_vtx", Trk_sec_vtx, &b_Trk_sec_vtx);
   fChain->SetBranchAddress("Trk_vxid", Trk_vxid, &b_Trk_vxid);
   fChain->SetBranchAddress("Trk_endpos", Trk_endpos, &b_Trk_endpos);
   fChain->SetBranchAddress("Trk_nzbyt", Trk_nzbyt, &b_Trk_nzbyt);
   fChain->SetBranchAddress("Trk_naxial", Trk_naxial, &b_Trk_naxial);
   fChain->SetBranchAddress("Trk_nstereo", Trk_nstereo, &b_Trk_nstereo);
   fChain->SetBranchAddress("Trk_ndof", Trk_ndof, &b_Trk_ndof);
   fChain->SetBranchAddress("Trk_layinner", Trk_layinner, &b_Trk_layinner);
   fChain->SetBranchAddress("Trk_layouter", Trk_layouter, &b_Trk_layouter);
   fChain->SetBranchAddress("Trk_dedxctd", Trk_dedxctd, &b_Trk_dedxctd);
   fChain->SetBranchAddress("Trk_dedxctdcr", Trk_dedxctdcr, &b_Trk_dedxctdcr);
   fChain->SetBranchAddress("Trk_dedxctderr", Trk_dedxctderr, &b_Trk_dedxctderr);
   fChain->SetBranchAddress("Trk_dedxctdnh", Trk_dedxctdnh, &b_Trk_dedxctdnh);
   fChain->SetBranchAddress("Trk_dedxctdsaturated", Trk_dedxctdsaturated, &b_Trk_dedxctdsaturated);
   fChain->SetBranchAddress("Trk_chi2", Trk_chi2, &b_Trk_chi2);
   fChain->SetBranchAddress("Trk_vchi2", Trk_vchi2, &b_Trk_vchi2);
   fChain->SetBranchAddress("Trk_imppar", Trk_imppar, &b_Trk_imppar);
   fChain->SetBranchAddress("Trk_imperr", Trk_imperr, &b_Trk_imperr);
   fChain->SetBranchAddress("Trk_pca", Trk_pca, &b_Trk_pca);
   fChain->SetBranchAddress("Filter", &Filter, &b_filter);
   fChain->SetBranchAddress("Knjets", &Knjets, &b_KNjets);
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
   fChain->SetBranchAddress("Knzufos", &Knzufos, &b_KNzufos);
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
   if(!Data) {
     fChain->SetBranchAddress("Photn", &Photn, &b_photN);
     fChain->SetBranchAddress("Photid", Photid, &b_Photid);
     fChain->SetBranchAddress("Photp", Photp, &b_Photp);
     fChain->SetBranchAddress("Npart", &Npart, &b_nPart);
     fChain->SetBranchAddress("Part_prt", Part_prt, &b_Part_prt);
     fChain->SetBranchAddress("Part_isthep", Part_isthep, &b_Part_isthep);
     fChain->SetBranchAddress("Part_charge", Part_charge, &b_Part_charge);
     fChain->SetBranchAddress("Part_id", Part_id, &b_Part_id);
     fChain->SetBranchAddress("Part_p", Part_p, &b_Part_p);
     fChain->SetBranchAddress("Mc_q2", &Mc_q2, &b_MC_Q2);
     fChain->SetBranchAddress("Mc_pfsl", Mc_pfsl, &b_MC_PFSL);
     fChain->SetBranchAddress("Part_jetid", Part_jetid, &b_Part_jetid);
     fChain->SetBranchAddress("Ktrnjets", &Ktrnjets, &b_KtrNJets);
     fChain->SetBranchAddress("Ktrjets", Ktrjets, &b_Ktrjets);
     fChain->SetBranchAddress("Mc_x", &Mc_x, &b_MC_X);
     fChain->SetBranchAddress("Mc_y", &Mc_y, &b_MC_Y);
   }
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


