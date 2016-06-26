#ifndef HIST_h
#define HIST_h
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TDirectory.h>
#include "constants.h"
class Hist
{
 public:
  Hist();
  ~Hist();
  void Init(TDirectory **dir);
  TH1D* zvtx;
  TH1D* dstw[2];
  TH1D* fltw[2];
  TH1D* sltw[6];
  TH1D* tltw[15];
  //for trigger efficiency study
  TH1D* theta_electron_tef[2];
  TH1D* e_electron_tef[2];
  TH1D* q2_tef[2];
  TH1D* e_gamma_tef[2];
  TH1D* eta_gamma_tef[2];
  TH1D* e_jet_tef[2];
  
  //for fmax & deltaz fit in different bins and cross sections:
  TH1D* det_cross_et;
  TH1D* det_cross_eta;
  TH1D* det_cross_Q2;
  TH1D* det_cross_x;
  TH1D* det_cross_et_jet;
  TH1D* det_cross_eta_jet;  

  TH1D* had_cross_et;
  TH1D* had_cross_eta;
  TH1D* had_cross_Q2;
  TH1D* had_cross_x;
  TH1D* had_cross_et_jet;
  TH1D* had_cross_et_jet2;
  TH1D* had_cross_eta_jet;  

  TH1D* part_cross_et;
  TH1D* part_cross_eta;
  TH1D* part_cross_Q2;
  TH1D* part_cross_x;
  TH1D* part_cross_et_jet;
  TH1D* part_cross_et_jet2;
  TH1D* part_cross_eta_jet;  

  TH1D* part_nohad_cross_et;
  TH1D* part_nohad_cross_eta;
  TH1D* part_nohad_cross_Q2;
  TH1D* part_nohad_cross_x;
  TH1D* part_nohad_cross_et_jet;
  TH1D* part_nohad_cross_et_jet2;
  TH1D* part_nohad_cross_eta_jet;  

  TH1D* had_nopart_cross_et;
  TH1D* had_nopart_cross_eta;
  TH1D* had_nopart_cross_Q2;
  TH1D* had_nopart_cross_x;
  TH1D* had_nopart_cross_et_jet;
  TH1D* had_nopart_cross_et_jet2;
  TH1D* had_nopart_cross_eta_jet;  

  TH1D* hd_cross_et;
  TH1D* hd_cross_eta;
  TH1D* hd_cross_Q2;
  TH1D* hd_cross_x;
  TH1D* hd_cross_et_jet;
  TH1D* hd_cross_eta_jet;  

  TProfile* prof_det_cross_et;
  TProfile* prof_det_cross_eta;
  TProfile* prof_det_cross_Q2;
  TProfile* prof_det_cross_x;
  TProfile* prof_det_cross_et_jet;
  TProfile* prof_det_cross_eta_jet;

  TProfile* prof_had_cross_et;
  TProfile* prof_had_cross_eta;
  TProfile* prof_had_cross_Q2;
  TProfile* prof_had_cross_x;
  TProfile* prof_had_cross_et_jet;
  TProfile* prof_had_cross_eta_jet;

  TProfile* prof_hd_cross_et;
  TProfile* prof_hd_cross_eta;
  TProfile* prof_hd_cross_Q2;
  TProfile* prof_hd_cross_x;
  TProfile* prof_hd_cross_et_jet;
  TProfile* prof_hd_cross_eta_jet;

  TH1D* fmax_et[number_etbins];
  TH1D* deltaz_et[number_etbins];
  TH1D* fmax_q2[number_Q2bins];
  TH1D* deltaz_q2[number_Q2bins];
  TH1D* fmax_x[number_xbins];
  TH1D* deltaz_x[number_xbins];
  TH1D* fmax_eta[number_etabins];
  TH1D* deltaz_eta[number_etabins];
  TH1D* fmax_et_jet[number_et_jetbins];
  TH1D* deltaz_et_jet[number_et_jetbins];
  TH1D* fmax_eta_jet[number_eta_jetbins];
  TH1D* deltaz_eta_jet[number_eta_jetbins];

  TH1D* had_prph_e;
  TH1D* had_jet_e;
  TH1D* had_Q2;
  TH1D* had_x;
  TH1D* had_et_jet;
  TH1D* had_eta_jet;

  TH1D* part_prph_e;
  TH1D* part_jet_e;
  TH1D* part_Q2;
  TH1D* part_x;
  TH1D* part_et_jet;
  TH1D* part_eta_jet;

  //fmax and deltaz for different particles;
  TH1D* fmax_exact_part[9];
  TH1D* deltaz_exact_part[9];

  //control distributions for isolated photon candidate
  TH1D* phcand_jet_deltaEta;
  TH1D* phcand_jet_deltaPhi;
  TH1D* phcand_jet_deltaR;
  
  //control distributions photon + jet
  TH1D* fmax;
  TH1D* deltaz;
  TH2D* fmax_deltaz;
  TH1D* phjet_electron_theta;
  TH1D* phjet_electron_phi;
  TH1D* phjet_electron_e;
  TH1D* phjet_electron_e_corr_burkard;
  TH1D* phjet_electron_probability;
  TH1D* prph_energy;
  TH1D* prph_eta;
  TH1D* prph_phi;
  TH1D* prph_jet_energy_frac;
  TH1D* prph_cell_energy_frac;
  TH1D* comp_prph_e;
  TH1D* comp_jet_e;
  TH1D* comp_Q2;
  TH1D* comp_x;
  TH1D* comp_y;
  TH1D* comp_s;
  TH2D* phjet_dphi_deta;
  TH1D* phjet_jet_et;
  TH1D* phjet_jet_eta;
  TH1D* phjet_jet_phi;
  TH1D* phjet_dphi;
  TH1D* phjet_deta;
  TH1D* phjet_dphi_el_ph;
  TH1D* phjet_deta_el_ph;
  TH1D* phjet_dphi_el_jet;
  TH1D* phjet_deta_el_jet;
  TH1D* et_jet_photon_ratio;
  TH2D* h2d_et_jet_et_photon;
  TH1D* et_3system_event_ratio;
  TH2D* h2d_phijet_phigamma;
  TH2D* h2d_dphijet_dphigamma;
  TH2D* h2d_dphijet_dphigamma_sharp;
  TH1D* phjet_cal_empz;
  TH1D* phjet_zufo_empz;
  TH1D* phjet_zvtx;
  TH2D* h2d_jet_npart_et;
  TH1D* phjet_photon_sincand_dr[4];
  TH2D* phjet2d_photon_sincand_dr_nr;
  TH1D* phjet_q2_dz_gt_035;
  TH1D* phjet_q2_dz_lt_035;
  TH1D* phjet_q2_dz_full;


  TH2D *h2d_phot_et_true_det;
  TH2D *h2d_phot_eta_true_det;
  TH2D *h2d_phot_phi_true_det;
  TH2D *h2d_corr_accjet_et_true_det;
  TH2D *h2d_corr_accjet_eta_true_det;
  TH2D *h2d_corr_accjet_phi_true_det;
  TH2D *h2d_uncorr_accjet_et_true_det;
  TH2D *h2d_uncorr_accjet_eta_true_det;
  TH2D *h2d_uncorr_accjet_phi_true_det;


  //control distributions inclusive photon
  TH1D* incl_phot_fmax;
  TH1D* incl_phot_deltaz;
  TH1D* incl_phot_elphot_dphi;
  TH1D* incl_phot_elphot_deta;
  TH1D* incl_phot_phi;
  TH1D* incl_phot_eta;
  //control distributions inclusive DIS
  TH1D* dis_electron_theta;
  TH1D* dis_electron_e;
  TH1D* dis_electron_probability;
  TH1D* dis_Q2;  
  TH2D* dis_Q2_x;
  TH1D* dis_Q2_el;
  TH1D* dis_Q2_da;
  TH1D* dis_Q2_da_cell;
  TH1D* dis_Q2_jb;
  TH1D* dis_Q2_jb_cell;
  TH1D* dis_Q2_true;
  TH1D* dis_x_true;
  TH1D* dis_y_true;
  TH1D* dis_electron_e_true;
  TH1D* dis_electron_theta_true;
  TH1D* dis_siyjb;
  TH1D* dis_siyjb_cell;
  TH1D* dis_handmade_yjb;
  TH1D* dis_elec_zufo_dr_min;
  TH1D* dis_elec_zufo_de_min;
  TH1D* dis_zufos_near_electron;

  TH1D* dis_zvtx;
  TH1D* dis_n_vtx_tracks;
  TH1D* dis_cal_empz;
  TH1D* dis_zufo_empz;
  TH1D* dis_pt_sqrtet;
  TH1D* dis_x_el;
  TH1D* dis_x_el_log;
  TH1D* dis_y_el;
  TH1D* x_gamma;
  TH1D* dis_px_event;
  TH1D* dis_py_event;
  TH1D* dis_pz_event;
  TH1D* dis_pt_event;
  TH1D* dis_e_event;
  TH1D* dis_et_event;
  TH1D* dis_n_zufos;
  TH1D* dis_emc;
  TH1D* dis_had;
  TH1D* dis_cal;
  TH1D* dis_had_tot_ratio;
  TH1D* dis_zufo_type;


  //MC correlation of generated and reconstructed values
  TH2D* mccorel_q2[3];
  TH2D* mccorel_x[3];
  TH2D* mccorel_y[3];

  TH2D* mccorel_q2_y;
  TH2D* mccorel_q2_x;

  TH2D* mccorel_q2_y_noweight;
  TH2D* mccorel_q2_x_noweight;

  //hadron level histograms
  TH1D* had_sum_e;
  TH1D* had_sum_px;
  TH1D* had_sum_py;
  TH1D* had_sum_pz;
  TH1D* had_sum_et;
  TH1D* had_q2;
  TH1D* had_y;

  //parton level histograms
  TH1D* part_nppart;
  TH1D* part_idpart;
  TH1D* part_sum_e;
  TH1D* part_sum_px;
  TH1D* part_sum_py;
  TH1D* part_sum_pz;
  TH1D* part_q2;
  TH1D* part_y;
  TH1D* part_sum_et;

  TH2D* h2d_etevent_part_had;
};

#endif
