#ifndef  constants_h
#define  constants_h
#include <TH1D.h>

//NEW CS VARIABLES PARAMETERS OF HISTOGRAMS
static  const Int_t number_Q2bins = 5;
static  const Int_t number_etabins = 4;
static  const Int_t number_etbins = 4;
static  const Int_t number_xbins = 4;
static  const Int_t number_et_jetbins = 6;//warning
static  const Int_t number_eta_jetbins = 4;

static  const Int_t number_xgamma_bins = 6; // 6: 1 - 0.4 ; 0.4 - 0.6; 0.6 - 0.7; 0.7-0.8; 0.8-0.9; 0.9-1;
static  const Int_t number_xp_bins = 6; //my 5: 0 - 0.002; 0.002 - 0.004; 0.004 - 0.006; 0.006 - 0.008; 0.008 - 0.01; 0.01 -0.03
static  const Int_t number_dphi_bins = 7;//5; //art 5: 0..45; 45..85; 85-120; 120-145; 145-180
static  const Int_t number_deta_bins = 6; //my 7:  -2.2 - -1.5; -1.5 - -0.8; -0.8 - -1.1;-1.1 - 0.6; 0.6 - 1.3; 1.3 - 2 ; 2 - 2.7
static  const Int_t number_dphi_e_ph_bins = 6; //
static  const Int_t number_deta_e_ph_bins = 5;//7; //
/*
  set_hist_atributes(70,  "h_phjet_dphi_el_ph_true_unmerged", "Photon_with_jet", "#Delta#phi_{e,#gamma}, per bin", 0.0, 3.5, 0.0, kFALSE, 25);
  set_hist_atributes(71,  "h_phjet_deta_el_ph_true_unmerged", "Photon_with_jet", "#Delta#eta_{e,#gamma}, per bin",  -4.5, 0., 0.0, kFALSE, 10);

  phjet_dphi_el_ph_true_unmerged = new TH1D("h_phjet_dphi_el_ph_true_unmerged", "#phi_{el} - #phi_{photon}", 200, -0.5, 3.5);//many unmerged bins
  phjet_deta_el_ph_true = new TH1D("h_phjet_deta_el_ph_true", "#eta_{el} - #eta_{photon}", 25, -3.6, -0.6);//merged bins
  phjet_deta_el_ph_true_unmerged = new TH1D("h_phjet_deta_el_ph_true_unmerged", "#eta_{el} - #eta_{photon}", 200, -7., 7.);//many unmerged bins
*/

static  const Int_t number_mbins = 4;
static  const Int_t number_xibins = 5;
static  const Int_t number_xbjbins = 5;
static const Int_t number_Q2bins_trijet = 5;

//*****************new variables***********************//
static Double_t xgamma_bin[number_xgamma_bins + 1] = {0, 0.4, 0.6, 0.7, 0.8, 0.9, 1};
static Double_t xp_bin[number_xp_bins + 1] = {0, 0.005, 0.01, 0.015, 0.025, 0.04, 0.07 };//{0, 0.004, 0.006, 0.008,  0.010 , 0.015 , 0.030};//{0, 0.004, 0.006, 0.008, 0.009, 0.03};//{0,  0.004, 0.006, 0.008, 0.01, 0.03};
static Double_t dphi_bin[number_dphi_bins + 1] = {0, 90, 130, 140, 150, 160, 170, 180};//{0., 45., 85., 120., 145., 180.};
static Double_t deta_bin[number_deta_bins + 1] = {-2.2, -1.5, -0.8, -0.1, 0.6, 1.3, 2.0};//{-2.2,  -0.8, -0.2, 0.4, 1.0, 1.5, 2.7};//{-2.2,  -0.8, -0.1, 0.6, 1.3, 2., 2.7};
static Double_t dphi_e_ph_bin[number_dphi_e_ph_bins + 1] = {0., 45., 80., 110., 135., 155., 180.};//{0., 30., 60., 90., 120., 150., 180.};
static Double_t deta_e_ph_bin[number_deta_e_ph_bins + 1] = {-3.6, -3.0, -2.4, -1.8, -1.2, -0.6 };//{-4.5, -3.6, -2.75, -2.25, -1.8, -1.2, -0.6, 0.};//{-4.5, -2.75, -2.25, -1.8, -1.2, 0.};//{-4.5, -3.0, -2.25, -1.5, -0.75, 0.};

// For xp it might be better to divide the upper bin:  0.010 - 0.015 - 0.030
// For delta-eta(egamma) it may be good to divide the outer bins each into
// two parts.  Could you remake with these changes for us to look at?


static Double_t et_bin[number_etbins+1] = {4., 6., 8., 10., 15.};
  static Double_t et_bin_err[number_etbins] = {1, 1, 1, 2.5};
static Double_t eta_bin_crosssec[number_etabins + 1] = {-0.7, -0.3, 0.1, 0.5, 0.9};
  static Double_t eta_bin_crosssec_err[number_etabins] = {0.2, 0.2, 0.2, 0.2};
//usual binning *bin1.root
//static Double_t Q2_bin[number_Q2bins+1] = {10., 40., 80., 150., 350.};
//static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
//static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 5., 10., 18., 35.};

//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.75, 0., 0.9, 1.8};

//optimal binning : *.root
//static Double_t Q2_bin[number_Q2bins+1] = {10., 17.65, 32.1, 66.1, 350.};
//static Double_t x_bin[number_xbins+1] = {1.5e-4, 6.5e-4, 1.3e-3, 2.8e-3, 0.15};
//static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 3.5, 5., 7.7, 35.};
//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, 0., 0.6, 1.2, 1.8};

//Peter's suggestion : *bin2.root FINAL!!!!
static Double_t Q2_bin[number_Q2bins+1] = {10., 20., 40., 80., 150., 350.};
  static Double_t Q2_bin_err[number_Q2bins] = {5., 10., 20., 35., 100.};
static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
  static Double_t x_bin_err[number_xbins] = {4.e-4, 1.e-3, 3.5e-3, 5.e-3};
static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 4., 6., 8., 10., 15., 35.};
  static Double_t et_jet_bin_err[number_et_jetbins] = {0.75, 1., 1., 1., 2.5, 10.};
static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.7, 0.1, 0.9, 1.8};
  static Double_t eta_jet_bin_err[number_eta_jetbins] = {0.4, 0.4, 0.4, 0.45};

//**********************************************//

static  Double_t y_Spiesberger_pos_nlo_et[number_etbins] = {2.320, 1.054, 0.522, 0.180}; // for paper
//static  Double_t y_Spiesberger_pos_nlo_et[number_etbins] = {2.5651, 1.1818, 0.5787, 0.2034}; // running aem
static  Double_t y_Spiesberger_pos_nlo_eta[number_etabins] = {6.728, 5.305, 5.027, 4.827};
static  Double_t y_Spiesberger_pos_nlo_Q2[number_Q2bins] = {0.2089, 0.1000, 0.0491, 0.0232, 0.00538};
static  Double_t y_Spiesberger_pos_nlo_x[number_xbins] = {3408, 1623, 344.1, 31.43};
static  Double_t y_Spiesberger_pos_nlo_et_jet[number_et_jetbins] = {1.297, 1.0375, 0.7803,  0.5558, 0.2829, 0.0281};
static  Double_t y_Spiesberger_pos_nlo_eta_jet[number_eta_jetbins] = {1.477, 2.954, 3.380, 2.780};

static  Double_t y_Spiesberger_pos_nlo_et_min[number_etbins] = {0.1420000, 0.0400000, 0.0234000, 0.0052000}; // for paper
//static  Double_t y_Spiesberger_pos_nlo_et_min[number_etbins] = {0.0505, 0.0309, 0.0108, 0.0056}; // running aem
static  Double_t y_Spiesberger_pos_nlo_eta_min[number_etabins] = {0.3920000, 0.2540000, 0.2430000, 0.2190000};
static  Double_t y_Spiesberger_pos_nlo_Q2_min[number_Q2bins] = {0.0140000, 0.0058700, 0.0021700, 0.0013300, 0.0002920};
static  Double_t y_Spiesberger_pos_nlo_x_min[number_xbins] = {233.0000000, 84.0000000, 25.0000000, 2.3000000};
static  Double_t y_Spiesberger_pos_nlo_et_jet_min[number_et_jetbins] = {0.1000000, 0.0633000, 0.0374000, 0.0216000, 0.0079000, 0.0003000};
static  Double_t y_Spiesberger_pos_nlo_eta_jet_min[number_eta_jetbins] = {0.1710000, 0.1810000, 0.2120000, 0.1310000};

static  Double_t y_Spiesberger_pos_nlo_et_max[number_etbins] = {0.1210000, 0.0630000, 0.0244000, 0.0078000}; // for paper
//static  Double_t y_Spiesberger_pos_nlo_et_max[number_etbins] = {0.0111, 0.0133, 0.0145, 0.0011}; // running aem
static  Double_t y_Spiesberger_pos_nlo_eta_max[number_etabins] = {0.3900000, 0.3040000, 0.1800000, 0.1840000};
static  Double_t y_Spiesberger_pos_nlo_Q2_max[number_Q2bins] = {0.0086000, 0.0043100, 0.0028300, 0.0012500, 0.0002900};
static  Double_t y_Spiesberger_pos_nlo_x_max[number_xbins] = {108.0000000, 77.0000000, 26.1000000, 2.6400000};
static  Double_t y_Spiesberger_pos_nlo_et_jet_max[number_et_jetbins] = {0.1070000, 0.0540000, 0.0368000, 0.0260000, 0.0096000, 0.0005000};
static  Double_t y_Spiesberger_pos_nlo_eta_jet_max[number_eta_jetbins] = {0.1480000, 0.1840000, 0.1380000, 0.0890000};

static  Double_t y_Spiesberger_ele_nlo_et[number_etbins] = {2.382, 1.114, 0.567, 0.202}; // for paper
//static  Double_t y_Spiesberger_ele_nlo_et[number_etbins] = {2.6250, 1.2254, 0.6290, 0.2246}; // running aem
static  Double_t y_Spiesberger_ele_nlo_eta[number_etabins] = {7.130, 5.626, 5.270, 4.981};
static  Double_t y_Spiesberger_ele_nlo_Q2[number_Q2bins] = {0.2148, 0.1043, 0.0519, 0.0246, 0.00576};
static  Double_t y_Spiesberger_ele_nlo_x[number_xbins] = {3475, 1710, 372.3, 33.99};
static  Double_t y_Spiesberger_ele_nlo_et_jet[number_et_jetbins] = {1.385, 1.122,  0.8381, 0.5747, 0.2866, 0.0260};
static  Double_t y_Spiesberger_ele_nlo_eta_jet[number_eta_jetbins] = {1.544, 3.066, 3.530, 2.994};

static  Double_t y_Spiesberger_ele_nlo_et_min[number_etbins] = {0.1290000, 0.0490000, 0.0247000, 0.0070000}; // for paper
//static  Double_t y_Spiesberger_ele_nlo_et_min[number_etbins] = {0.0303, 0.0110, 0.0175, 0.0019}; // running aem
static  Double_t y_Spiesberger_ele_nlo_eta_min[number_etabins] = {0.4290000, 0.2530000, 0.2160000, 0.1930000};
static  Double_t y_Spiesberger_ele_nlo_Q2_min[number_Q2bins] = {0.0146000, 0.0061100, 0.0022800, 0.0012800, 0.0003270};
static  Double_t y_Spiesberger_ele_nlo_x_min[number_xbins] = {239.0000000, 84.0000000, 29.4000000, 2.4700000};
static  Double_t y_Spiesberger_ele_nlo_et_jet_min[number_et_jetbins] = {0.1030000, 0.0650000, 0.0286000, 0.0160000, 0.0110000, 0.0007200};
static  Double_t y_Spiesberger_ele_nlo_eta_jet_min[number_eta_jetbins] = {0.2050000, 0.2090000, 0.1920000, 0.1550000};

static  Double_t y_Spiesberger_ele_nlo_et_max[number_etbins] = {0.1320000, 0.0540000, 0.0264000, 0.0072000}; // for paper
//static  Double_t y_Spiesberger_ele_nlo_et_max[number_etbins] = {0.0273, 0.0225, 0.0027, 0.0044}; // running aem
static  Double_t y_Spiesberger_ele_nlo_eta_max[number_etabins] = {0.3390000, 0.2860000, 0.2230000, 0.2170000};
static  Double_t y_Spiesberger_ele_nlo_Q2_max[number_Q2bins] = {0.0086000, 0.0045000, 0.0028900, 0.0014800, 0.0003130};
static  Double_t y_Spiesberger_ele_nlo_x_max[number_xbins] = {116.0000000, 74.0000000, 27.0000000, 2.5100000};
static  Double_t y_Spiesberger_ele_nlo_et_jet_max[number_et_jetbins] = {0.1080000, 0.0560000, 0.0363000, 0.0259000, 0.0107000, 0.0003900};
static  Double_t y_Spiesberger_ele_nlo_eta_jet_max[number_eta_jetbins] = {0.1440000, 0.2100000, 0.1570000, 0.0540000};

//
// here come results obtained with running aem
//

// HERAPDF1.0
/* static  Double_t y_Spiesberger_pos_nlo_ra_et[number_etbins] = {2.54605E+000,1.19055E+000,5.92450E-001,2.03460E-001}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_eta[number_etabins] = {7.4335, 5.8875, 5.58975, 5.29525}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_Q2[number_Q2bins] = {0.23303, 0.111245, 0.0543575, 0.0255742857142857, 0.0059355}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_x[number_xbins] = {3833.75, 1804.8, 376.157142857143, 34.477}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_et_jet[number_et_jetbins] = {1.4264, 1.14915, 0.86305, 0.61005, 0.3175, 0.032177}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_eta_jet[number_eta_jetbins] = {1.60525, 3.28575, 3.7325, 3.12911111111111}; */

/* static  Double_t y_Spiesberger_ele_nlo_ra_et[number_etbins] = {2.61190E+000,1.25150E+000,6.33250E-001,2.26620E-001}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_eta[number_etabins] = {7.852, 6.2275, 5.8245, 5.5495}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_Q2[number_Q2bins] = {0.23892, 0.115905, 0.0572375, 0.0272342857142857, 0.0063785}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_x[number_xbins] = {3901.125, 1897.45, 406.542857142857, 37.365}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_et_jet[number_et_jetbins] = {1.52786666666667, 1.23375, 0.9209, 0.6367, 0.32124, 0.029575}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_eta_jet[number_eta_jetbins] = {1.66175, 3.38675, 3.925875, 3.34788888888889}; */

/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_et[number_etbins] = {2.55655E+000,1.19519E+000,5.88445E-001,2.03834E-001}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_eta[number_etabins] = {7.4335, 5.8875, 5.58975, 5.29525}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_Q2[number_Q2bins] = {0.23303, 0.111245, 0.0543575, 0.0255742857142857, 0.0059355}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_x[number_xbins] = {3833.75, 1804.8, 376.157142857143, 34.477}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_et_jet[number_et_jetbins] = {1.4264, 1.14915, 0.86305, 0.61005, 0.3175, 0.032177}; */
/* static  Double_t y_Spiesberger_pos_nlo_ra_fast_eta_jet[number_eta_jetbins] = {1.60525, 3.28575, 3.7325, 3.12911111111111}; */

/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_et[number_etbins] = {2.63028E+000,1.25297E+000,6.36125E-001,2.27486E-001}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_eta[number_etabins] = {7.852, 6.2275, 5.8245, 5.5495}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_Q2[number_Q2bins] = {0.23892, 0.115905, 0.0572375, 0.0272342857142857, 0.0063785}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_x[number_xbins] = {3901.125, 1897.45, 406.542857142857, 37.365}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_et_jet[number_et_jetbins] = {1.52786666666667, 1.23375, 0.9209, 0.6367, 0.32124, 0.029575}; */
/* static  Double_t y_Spiesberger_ele_nlo_ra_fast_eta_jet[number_eta_jetbins] = {1.66175, 3.38675, 3.925875, 3.34788888888889}; */

//CTEQ6.6 (et ONLY!!!)
static  Double_t y_Spiesberger_pos_nlo_ra_et[number_etbins] = {2.47800E+000,1.14940E+000,5.67750E-001,1.96440E-001};
static  Double_t y_Spiesberger_pos_nlo_ra_eta[number_etabins] = {7.4335, 5.8875, 5.58975, 5.29525};
static  Double_t y_Spiesberger_pos_nlo_ra_Q2[number_Q2bins] = {0.23303, 0.111245, 0.0543575, 0.0255742857142857, 0.0059355};
static  Double_t y_Spiesberger_pos_nlo_ra_x[number_xbins] = {3833.75, 1804.8, 376.157142857143, 34.477};
static  Double_t y_Spiesberger_pos_nlo_ra_et_jet[number_et_jetbins] = {1.4264, 1.14915, 0.86305, 0.61005, 0.3175, 0.032177};
static  Double_t y_Spiesberger_pos_nlo_ra_eta_jet[number_eta_jetbins] = {1.60525, 3.28575, 3.7325, 3.12911111111111};

static  Double_t y_Spiesberger_ele_nlo_ra_et[number_etbins] = {2.54720E+000,1.19550E+000,6.09150E-001,2.17460E-001};
static  Double_t y_Spiesberger_ele_nlo_ra_eta[number_etabins] = {7.852, 6.2275, 5.8245, 5.5495};
static  Double_t y_Spiesberger_ele_nlo_ra_Q2[number_Q2bins] = {0.23892, 0.115905, 0.0572375, 0.0272342857142857, 0.0063785};
static  Double_t y_Spiesberger_ele_nlo_ra_x[number_xbins] = {3901.125, 1897.45, 406.542857142857, 37.365};
static  Double_t y_Spiesberger_ele_nlo_ra_et_jet[number_et_jetbins] = {1.52786666666667, 1.23375, 0.9209, 0.6367, 0.32124, 0.029575};
static  Double_t y_Spiesberger_ele_nlo_ra_eta_jet[number_eta_jetbins] = {1.66175, 3.38675, 3.925875, 3.34788888888889};

static  Double_t y_Spiesberger_pos_nlo_ra_fast_et[number_etbins] = {2.48224E+000,1.15623E+000,5.67670E-001,1.95962E-001};
static  Double_t y_Spiesberger_pos_nlo_ra_fast_eta[number_etabins] = {7.4335, 5.8875, 5.58975, 5.29525};
static  Double_t y_Spiesberger_pos_nlo_ra_fast_Q2[number_Q2bins] = {0.23303, 0.111245, 0.0543575, 0.0255742857142857, 0.0059355};
static  Double_t y_Spiesberger_pos_nlo_ra_fast_x[number_xbins] = {3833.75, 1804.8, 376.157142857143, 34.477};
static  Double_t y_Spiesberger_pos_nlo_ra_fast_et_jet[number_et_jetbins] = {1.4264, 1.14915, 0.86305, 0.61005, 0.3175, 0.032177};
static  Double_t y_Spiesberger_pos_nlo_ra_fast_eta_jet[number_eta_jetbins] = {1.60525, 3.28575, 3.7325, 3.12911111111111};

static  Double_t y_Spiesberger_ele_nlo_ra_fast_et[number_etbins] = {2.54852E+000,1.20789E+000,6.10730E-001,2.17412E-001};
static  Double_t y_Spiesberger_ele_nlo_ra_fast_eta[number_etabins] = {7.852, 6.2275, 5.8245, 5.5495};
static  Double_t y_Spiesberger_ele_nlo_ra_fast_Q2[number_Q2bins] = {0.23892, 0.115905, 0.0572375, 0.0272342857142857, 0.0063785};
static  Double_t y_Spiesberger_ele_nlo_ra_fast_x[number_xbins] = {3901.125, 1897.45, 406.542857142857, 37.365};
static  Double_t y_Spiesberger_ele_nlo_ra_fast_et_jet[number_et_jetbins] = {1.52786666666667, 1.23375, 0.9209, 0.6367, 0.32124, 0.029575};
static  Double_t y_Spiesberger_ele_nlo_ra_fast_eta_jet[number_eta_jetbins] = {1.66175, 3.38675, 3.925875, 3.34788888888889};

/// end of running aem

static  Double_t y_Spiesberger_pos_lo_et[number_etbins] = {2.527, 1.191, 0.593, 0.202};
static  Double_t y_Spiesberger_pos_lo_eta[number_etabins] = {7.697, 5.853, 5.382, 5.159};
static  Double_t y_Spiesberger_pos_lo_Q2[number_Q2bins] = {0.1979, 0.1002, 0.0574, 0.0286, 0.00679};
static  Double_t y_Spiesberger_pos_lo_x[number_xbins] = {3037, 1736, 462.5, 43.58};
static  Double_t y_Spiesberger_pos_lo_et_jet[number_et_jetbins] = {1.464, 1.111, 0.8830, 0.6402, 0.3149, 0.0297};
static  Double_t y_Spiesberger_pos_lo_eta_jet[number_eta_jetbins] = {2.354, 3.631, 3.348, 2.406};

static  Double_t y_Spiesberger_pos_lo_et_min[number_etbins] = {0.4450000, 0.1870000, 0.0804000, 0.0212000};
static  Double_t y_Spiesberger_pos_lo_eta_min[number_etabins] = {1.2980000, 0.9910000, 0.8500000, 0.7050000};
static  Double_t y_Spiesberger_pos_lo_Q2_min[number_Q2bins] = {0.0359000, 0.0177800, 0.0098200, 0.0040700, 0.0007080};
static  Double_t y_Spiesberger_pos_lo_x_min[number_xbins] = {525.0000000, 292.0000000, 68.4000000, 4.1100000};
static  Double_t y_Spiesberger_pos_lo_et_jet_min[number_et_jetbins] = {0.3440000, 0.2059000, 0.1338000, 0.0804000, 0.0295000, 0.0014500};
static  Double_t y_Spiesberger_pos_lo_eta_jet_min[number_eta_jetbins] = {0.4380000, 0.6210000, 0.5160000, 0.3040000};

static  Double_t y_Spiesberger_pos_lo_et_max[number_etbins] = {0.3770000, 0.1540000, 0.0650000, 0.0167000};
static  Double_t y_Spiesberger_pos_lo_eta_max[number_etabins] = {1.0930000, 0.8320000, 0.6940000, 0.5650000};
static  Double_t y_Spiesberger_pos_lo_Q2_max[number_Q2bins] = {0.0303000, 0.0148000, 0.0081600, 0.0033300, 0.0005530};
static  Double_t y_Spiesberger_pos_lo_x_max[number_xbins] = {449.0000000, 246.0000000, 55.4000000, 3.0700000};
static  Double_t y_Spiesberger_pos_lo_et_jet_max[number_et_jetbins] = {0.2890000, 0.1710000, 0.1116000, 0.0657000, 0.0238000, 0.0010600};
static  Double_t y_Spiesberger_pos_lo_eta_jet_max[number_eta_jetbins] = {0.3800000, 0.5270000, 0.4230000, 0.2310000};

static  Double_t y_Spiesberger_ele_lo_et[number_etbins] = {2.622, 1.262, 0.642, 0.227};
static  Double_t y_Spiesberger_ele_lo_eta[number_etabins] = {8.125, 6.239, 5.702, 5.416};
static  Double_t y_Spiesberger_ele_lo_Q2[number_Q2bins] = {0.2047, 0.1055, 0.0610, 0.0305, 0.00729};
static  Double_t y_Spiesberger_ele_lo_x[number_xbins] = {3115, 1839, 498.7, 46.96};
static  Double_t y_Spiesberger_ele_lo_et_jet[number_et_jetbins] = {1.568, 1.207, 0.9540, 0.6782, 0.3223, 0.0273};
static  Double_t y_Spiesberger_ele_lo_eta_jet[number_eta_jetbins] = {2.443, 3.791, 3.564, 2.610};

static  Double_t y_Spiesberger_ele_lo_et_min[number_etbins] = {0.4540000, 0.1930000, 0.0831000, 0.0228000};
static  Double_t y_Spiesberger_ele_lo_eta_min[number_etabins] = {1.3350000, 1.0180000, 0.8630000, 0.7210000};
static  Double_t y_Spiesberger_ele_lo_Q2_min[number_Q2bins] = {0.0366000, 0.0182300, 0.0100900, 0.0042000, 0.0007290};
static  Double_t y_Spiesberger_ele_lo_x_min[number_xbins] = {533.0000000, 301.0000000, 70.8000000, 4.2100000};
static  Double_t y_Spiesberger_ele_lo_et_jet_min[number_et_jetbins] = {0.3540000, 0.2127000, 0.1379000, 0.0826000, 0.0297000, 0.0013600};
static  Double_t y_Spiesberger_ele_lo_eta_jet_min[number_eta_jetbins] = {0.4430000, 0.6310000, 0.5290000, 0.3170000};

static  Double_t y_Spiesberger_ele_lo_et_max[number_etbins] = {0.3830000, 0.1580000, 0.0677000, 0.0177000};
static  Double_t y_Spiesberger_ele_lo_eta_max[number_etabins] = {1.1150000, 0.8530000, 0.7100000, 0.5760000};
static  Double_t y_Spiesberger_ele_lo_Q2_max[number_Q2bins] = {0.0309000, 0.0152000, 0.0084000, 0.0034200, 0.0005780};
static  Double_t y_Spiesberger_ele_lo_x_max[number_xbins] = {456.0000000, 252.0000000, 57.2000000, 3.1400000};
static  Double_t y_Spiesberger_ele_lo_et_jet_max[number_et_jetbins] = {0.2990000, 0.1770000, 0.1139000, 0.0655000, 0.0239000, 0.0010700};
static  Double_t y_Spiesberger_ele_lo_eta_jet_max[number_eta_jetbins] = {0.3870000, 0.5380000, 0.4370000, 0.2390000};

static  Double_t y_Zotov_et[number_etbins] = {3.495248, 1.881006, 1.075966, 0.449562};
static  Double_t y_Zotov_eta[number_etabins] = {11.360027, 10.052792, 8.805129, 7.662673};
static  Double_t y_Zotov_Q2[number_Q2bins] = {0.509965, 0.214755, 0.089220, 0.030072, 0.005290};
static  Double_t y_Zotov_x[number_xbins] = {8050.307500, 3054.110500, 537.109800, 41.495265};
static  Double_t y_Zotov_et_jet[number_et_jetbins] = {1.895937, 1.783060, 1.491320, 1.090219, 0.532862, 0.045571};
static  Double_t y_Zotov_eta_jet[number_eta_jetbins] = {4.748785, 5.722657, 4.643728, 3.400125};

static  Double_t y_Zotov_et_eyl[number_etbins] = {0.341544, 0.257052, 0.174815, 0.087082};//0.077259, 0.019416, 0.00582, 0.0000084};
static  Double_t y_Zotov_eta_eyl[number_etabins] = {1.551864, 1.387609, 1.116696, 0.899415};//{0.28564, 0.13898, 0.172415, 0.076321};
static  Double_t y_Zotov_Q2_eyl[number_Q2bins] = {0.108413, 0.029199, 0.007544, 0.001766, 0.000234};//{0.0042655, 0.0028603, 0.00343944, 0.00148881, 0.000175671};
static  Double_t y_Zotov_x_eyl[number_xbins] = {1299.269500, 404.222500, 45.759700, 2.147375};
static  Double_t y_Zotov_et_jet_eyl[number_et_jetbins] = {0.291646, 0.247816, 0.193900, 0.127446, 0.061705, 0.004948};
static  Double_t y_Zotov_eta_jet_eyl[number_eta_jetbins] = {0.597137, 0.603969, 0.613762, 0.589267};

static  Double_t y_Zotov_et_eyh[number_etbins] = {0.341544, 0.257052, 0.174815, 0.087082};//0.103759, 0.031662, 0.00958, 0.0041325};
static  Double_t y_Zotov_eta_eyh[number_etabins] = {1.551864, 1.387609, 1.116696, 0.899415};//0.09112, 0.15142, 0.109484, 0.160555};
static  Double_t y_Zotov_Q2_eyh[number_Q2bins] = {0.108413, 0.029199, 0.007544, 0.001766, 0.000234};//0.0065276, 0.0016325, 0.00262898, 0.00126356, 0.000214477};
static  Double_t y_Zotov_x_eyh[number_xbins] = {1299.269500, 404.222500, 45.759700, 2.147375};
static  Double_t y_Zotov_et_jet_eyh[number_et_jetbins] = {0.291646, 0.247816, 0.193900, 0.127446, 0.061705, 0.004948};
static  Double_t y_Zotov_eta_jet_eyh[number_eta_jetbins] = {0.597137, 0.603969, 0.613762, 0.589267};

/* static Double_t eyl_syst_et[number_etbins] = {0.330169022,0.168249286,0.087743121,0.030058075}; */
/* static Double_t eyl_syst_eta[number_etabins] = {1.189436548,0.763698993,0.687645812,0.685620475}; */
/* static Double_t eyl_syst_Q2[number_Q2bins] = {0.036086494,0.016873383,0.010891640,0.003702754,0.000990727}; */
/* static Double_t eyl_syst_x[number_xbins] = {511.936872904,207.904061741,65.059802668,8.384118618}; */
/* static Double_t eyl_syst_et_jet[number_et_jetbins] = {0.263185938,0.140872205,0.119718166,0.103524212,0.032934415,0.005581137}; */
/* static Double_t eyl_syst_eta_jet[number_eta_jetbins] = {0.279390688,0.324957743,0.433653925,0.482568274}; */

/* static Double_t eyh_syst_et[number_etbins] = {0.316130344,0.124424408,0.083197857,0.030862840}; */
/* static Double_t eyh_syst_eta[number_etabins] = {1.016515079,0.715740975,0.685390308,0.588175306}; */
/* static Double_t eyh_syst_Q2[number_Q2bins] = {0.038957003,0.018060491,0.006273332,0.003419672,0.000786238}; */
/* static Double_t eyh_syst_x[number_xbins] = {534.704606970,217.183968304,45.861072485,9.031963124}; */
/* static Double_t eyh_syst_et_jet[number_et_jetbins] = {0.236725452,0.133781457,0.119310615,0.085284157,0.052627055,0.005354525}; */
/* static Double_t eyh_syst_eta_jet[number_eta_jetbins] = {0.245547853,0.306921828,0.431799594,0.397256927}; */

static Double_t eyl_syst_et[number_etbins] = {0.221894006, 0.118678762, 0.087901800, 0.034040892};
static Double_t eyl_syst_eta[number_etabins] = {0.797297321, 0.635746988, 0.626181818, 0.651630915};
static Double_t eyl_syst_Q2[number_Q2bins] = {0.030775081, 0.015025327, 0.006774138, 0.002560945, 0.000684151};
static Double_t eyl_syst_x[number_xbins] = {457.005501488, 173.886629364, 33.941110078, 7.698774575};
static Double_t eyl_syst_et_jet[number_et_jetbins] = {0.183548169, 0.145839524, 0.121718906, 0.087550653, 0.036331115, 0.006231391};
static Double_t eyl_syst_eta_jet[number_eta_jetbins] = {0.231408677, 0.313109375, 0.359690027, 0.364821034}; 

  static Double_t eyl_syst_xgamma[number_xgamma_bins] =       {0,  0., 0., 0.,  0., 0.};
  static Double_t eyl_syst_xp[number_xp_bins] =               {0,  0., 0., 0.,  0., 0.};
  static Double_t eyl_syst_dphi[number_dphi_bins] =           {0., 0., 0., 0.,  0., 0., 0.};
  static Double_t eyl_syst_deta[number_deta_bins] =           {0,  0., 0., 0.,  0., 0.};
  static Double_t eyl_syst_dphi_e_ph[number_dphi_e_ph_bins] = {0., 0., 0., 0.,  0., 0.};
  static Double_t eyl_syst_deta_e_ph[number_deta_e_ph_bins] = {0., 0., 0., 0.,  0.};// {0., 0., 0., 0.,  0., 0, 0};

static Double_t eyh_syst_et[number_etbins] = {0.221894006, 0.118678762, 0.087901800, 0.034040892};
static Double_t eyh_syst_eta[number_etabins] = {0.797297321, 0.635746988, 0.626181818, 0.651630915};
static Double_t eyh_syst_Q2[number_Q2bins] = {0.030775081, 0.015025327, 0.006774138, 0.002560945, 0.000684151};
static Double_t eyh_syst_x[number_xbins] = {457.005501488, 173.886629364, 33.941110078, 7.698774575};
static Double_t eyh_syst_et_jet[number_et_jetbins] = {0.183548169, 0.145839524, 0.121718906, 0.087550653, 0.036331115, 0.006231391};
static Double_t eyh_syst_eta_jet[number_eta_jetbins] = {0.231408677, 0.313109375, 0.359690027, 0.364821034}; 

  static Double_t eyh_syst_xgamma[number_xgamma_bins] =     {0,  0., 0., 0.,  0., 0.};
  static Double_t eyh_syst_xp[number_xp_bins] =         {0,  0., 0., 0.,  0., 0.};
  static Double_t eyh_syst_dphi[number_dphi_bins] =       {0., 0., 0., 0.,  0., 0., 0.};
  static Double_t eyh_syst_deta[number_deta_bins] =       {0,  0., 0., 0.,  0., 0.};
  static Double_t eyh_syst_dphi_e_ph[number_dphi_e_ph_bins] = {0., 0., 0., 0.,  0., 0.};
  static Double_t eyh_syst_deta_e_ph[number_deta_e_ph_bins] = {0., 0., 0., 0.,  0.};// {0., 0., 0., 0.,  0., 0, 0};

/* static Double_t y_1stcross_et[number_etbins] = {2.41818500, 1.31228600, 0.62352600, 0.24772540}; */
/* static Double_t y_1stcross_eta[number_etabins] = {7.35327000, 6.46609000, 5.55047800, 5.12890200}; */
/* static Double_t y_1stcross_Q2[number_Q2bins] = {0.29033640, 0.12958410, 0.04843650, 0.02239158, 0.00368611}; */
/* static Double_t y_1stcross_x[number_xbins] = {4695.07, 1795.362, 271.719, 26.63652}; */
/* static Double_t y_1stcross_eta_jet[number_eta_jetbins] = {1.48587900, 2.68198900, 3.84862700, 3.68495}; */
/* static Double_t y_1stcross_et_jet[number_et_jetbins] = {1.48100300, 1.15483400, 0.97312100, 0.68623900, 0.29178700, 0.02654470}; */

//systematics
  static Double_t y_1stcross_et[number_etbins] = {2.380595, 1.281126, 0.615687, 0.2556284};
  static Double_t y_1stcross_eta[number_etabins] = {7.56376, 6.7092, 5.833198, 5.155802};
  static Double_t y_1stcross_Q2[number_Q2bins] = {0.2978664, 0.1291481, 0.0486016, 0.02243448, 0.00374546};
  static Double_t y_1stcross_x[number_xbins] = {4869.363, 1811.322, 278.144, 25.08412};
  static Double_t y_1stcross_et_jet[number_et_jetbins] = {1.397993, 1.192586, 1.008039, 0.742206, 0.315349, 0.0305625};
  static Double_t y_1stcross_eta_jet[number_eta_jetbins] = {1.529209, 2.844649, 3.910317, 3.57312}; 

  static Double_t y_1stcross_xgamma[number_xgamma_bins] =     {0,  0., 0., 0.,  0., 0.};
  static Double_t y_1stcross_xp[number_xp_bins] =         {0,  0., 0., 0.,  0., 0.};
  static Double_t y_1stcross_dphi[number_dphi_bins] =       {0., 0., 0., 0.,  0., 0., 0.};
  static Double_t y_1stcross_deta[number_deta_bins] =       {0,  0., 0., 0.,  0., 0.};
  static Double_t y_1stcross_dphi_e_ph[number_dphi_e_ph_bins] = {0., 0., 0., 0.,  0., 0.};
  static Double_t y_1stcross_deta_e_ph[number_deta_e_ph_bins] = {0., 0., 0., 0.,  0.};// {0., 0., 0., 0.,  0., 0, 0};

static Double_t y_1stcross_stat_err_et[number_etbins] = {0.176738497667854, 0.101735141240576, 0.0757867780543711, 0.0250895660515701};
static Double_t y_1stcross_stat_err_eta[number_etabins] = {0.592805180136139, 0.535135083185741, 0.552219412138907, 0.528456345547378};
static Double_t y_1stcross_stat_err_Q2[number_Q2bins] = {0.0241382608376184, 0.0116559102705781, 0.00522137006664161, 0.00231890492485138, 0.000652872072841993};
static Double_t y_1stcross_stat_err_x[number_xbins] = {333.870254537837, 139.355179833726, 31.4037400209657, 6.87083375343924};
static Double_t y_1stcross_stat_err_et_jet[number_et_jetbins] = {0.163429798765586, 0.108232420004359, 0.0976833027625499, 0.0722924934799188, 0.0320890488381956, 0.00565461683652827};
static Double_t y_1stcross_stat_err_eta_jet[number_eta_jetbins] = {0.172542804225502, 0.246457814901476, 0.32935485426664, 0.290544835217837}; 

//? corrections on hadronisations
static Double_t hadronisation_et[number_etbins] = {0.838210272, 0.889243789, 0.935850904, 0.990858773};
static Double_t hadronisation_eta[number_etabins] = {0.931679121, 0.906123134, 0.856993579, 0.812842358};
static Double_t hadronisation_Q2[number_Q2bins] = {0.925662245, 0.890410989, 0.841515777, 0.855468932, 0.877568605};
static Double_t hadronisation_x[number_xbins] = {0.960521282, 0.863584534, 0.822940578, 0.851649185};
static Double_t hadronisation_et_jet[number_et_jetbins] = {0.817530672, 0.829958746, 0.865591600, 0.997190900, 1.019911904, 1.072649139};
static Double_t hadronisation_eta_jet[number_eta_jetbins] = {0.656861428, 0.806478183, 0.961334927, 1.071164515};

//Peter's suggestion : *bin2.root jet_jet > 4 GeV FINAL!!!!
//static Double_t Q2_bin[number_Q2bins+1] = {10., 20., 40., 80., 150., 350.};
//static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
//static Double_t et_jet_bin[number_et_jetbins+1] = {4., 6., 8., 10., 15., 35.};
//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.7, 0.1, 0.9, 1.8};

//Peter's suggestion : *bin2_extend.root
//static Double_t Q2_bin[number_Q2bins+1] = {10., 40., 80., 150., 350.};
//static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
//static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 4., 6., 8., 10., 15., 35.};
//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.7, 0.1, 0.9, 1.8};

//Peter's suggestion + optimal + round: *bin3.root
//static Double_t Q2_bin[number_Q2bins+1] = {10., 40., 80., 150., 350.};
//static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
//static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 3.5, 5., 7.5, 35.};
//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.7, 0.1, 0.9, 1.8};

//optimal round binning for 5 bins for et_jet
//2.5 3.2  4.2  5.9 8.6 35.
//2.5   3.   4.   5.5    8.5   35

//optimal round binning for 6 bins for et_jet *bin4.root
//2.5 3  3.8 5  6.6  9.2 35.
//2.5   3.   4.   5.    6.5   9. 35.
//static Double_t Q2_bin[number_Q2bins+1] = {10., 40., 80., 150., 350.};
//static Double_t x_bin[number_xbins+1] = {2.e-4, 1.e-3, 3.e-3, 1.e-2, 2.e-2};
//static Double_t et_jet_bin[number_et_jetbins+1] = {2.5, 3., 4., 5., 6.5, 9., 35.};
//static Double_t eta_jet_bin[number_eta_jetbins+1] = {-1.5, -0.7, 0.1, 0.9, 1.8};
  


static Double_t m_bin[number_mbins+1] = {20., 30., 45., 65., 120.};
static Double_t xi_bin[number_xibins+1] = {-2.0, -1.6, -1.45, -1.3, -1.1, -0.0};
static Double_t xbj_bin[number_xbjbins +1] = {0.0001,0.01,0.02,0.035,0.07,0.1};


static Double_t Q2_bin_trijet[number_Q2bins_trijet+1] = {125., 250., 780., 1300., 5000., 20000.};

 //double differential cross sections
  static const Int_t xidiffq2bins2_1 = 4;
  static Double_t xi_diff_q2_bins2_1[xidiffq2bins2_1 + 1] = {-2.1, -1.65, -1.5, -1.3, -0.4};
  
  static const Int_t xidiffq2bins2_2 = 4;
   static Double_t xi_diff_q2_bins2_2[xidiffq2bins2_2 + 1] = {-2.0, -1.55, -1.4, -1.25, -0.4};
  
  static const Int_t xidiffq2bins2_3 = 4;
   static Double_t xi_diff_q2_bins2_3[xidiffq2bins2_3 + 1] = {-1.9, -1.45, -1.3, -1.15, -0.4};
  
  static const Int_t xidiffq2bins2_4 = 4;
   static Double_t xi_diff_q2_bins2_4[xidiffq2bins2_4 + 1] = {-1.7, -1.25, -1.15, -1.0, -0.25};
  
  static const Int_t xidiffq2bins2_5 = 3;
   static Double_t xi_diff_q2_bins2_5[xidiffq2bins2_5 + 1] = {-1.5, -1.0, -0.85, -0.2};
  
  static const Int_t xidiffq2bins2_6 = 3;
   static Double_t xi_diff_q2_bins2_6[xidiffq2bins2_6 + 1] = {-1.1, -0.75, -0.55, -0.0};
  
  //and et in bins of q2
  
  static const Int_t etmeandiffq2bins2_1 = 4;
   static Double_t etmean_diff_q2_bins2_1[etmeandiffq2bins2_1 + 1] = {8., 15., 22., 30., 60.};
  
  static const Int_t etmeandiffq2bins2_2 = 4;
   static Double_t etmean_diff_q2_bins2_2[etmeandiffq2bins2_2 + 1] = {8., 15., 22., 30., 60.};
  
  static const Int_t etmeandiffq2bins2_3 = 4;
   static Double_t etmean_diff_q2_bins2_3[etmeandiffq2bins2_3 + 1] = {8., 15., 22., 30., 60.};
  
  static const Int_t etmeandiffq2bins2_4 = 4;
   static Double_t etmean_diff_q2_bins2_4[etmeandiffq2bins2_4 + 1] = {8., 15., 22., 30., 60.};
  
  static const Int_t etmeandiffq2bins2_5 = 3;
   static Double_t etmean_diff_q2_bins2_5[etmeandiffq2bins2_5 + 1] = {8., 16., 28., 60.};
  
  static const Int_t etmeandiffq2bins2_6 = 3;
   static Double_t etmean_diff_q2_bins2_6[etmeandiffq2bins2_6 + 1] = {8., 16., 28., 60.};


 static const Int_t number_of_eta_bins = 14;
 static Float_t eta_bin[number_of_eta_bins+1] = {-1., -0.75, -0.5, -0.25, 0.,
                                                 0.25, 0.5, 0.75, 1., 1.25,
                                                 1.5, 1.75, 2., 2.25, 2.5};

//start-lepto-989900-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:54 CET 2010
static Double_t lepto_factor_cells_989900[number_of_eta_bins] = { 
1.02260671929943236513, 1.02260671929943236513, 1.02260671929943236513, 1.02260671929943236513,
1.01258971448100476564, 1.01258971448100476564, 1.01258971448100476564, 1.01258971448100476564, 
0.98616665068256514903, 0.98616665068256514903, 
1.00857027963993894915, 1.00857027963993894915, 
1.00704832477663375201, 1.00704832477663375201 
};
//stop-lepto-989900-hadronic-scale-correction-cells//

//start-ariadne-989900-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:53 CET 2010
static Double_t ariadne_factor_cells_989900[number_of_eta_bins] = {
  1.02004378152046326633, 1.02004378152046326633, 1.02004378152046326633, 1.02004378152046326633, 
  1.00735115902074379157, 1.00735115902074379157, 1.00735115902074379157, 1.00735115902074379157, 
  0.98306107405028775315, 0.98306107405028775315, 
  1.00724762641653353690, 1.00724762641653353690, 
  1.00867751589144805102, 1.00867751589144805102
};
//stop-ariadne-989900-hadronic-scale-correction-cells//

//start-lepto-0405e-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:56 CET 2010 
static Double_t lepto_factor_cells_0405e[number_of_eta_bins] = {
1.03292940879942185539, 1.03292940879942185539, 1.03292940879942185539, 1.03292940879942185539, 
1.01401317900758347612, 1.01401317900758347612, 1.01401317900758347612, 1.01401317900758347612, 
0.98913114824820624893, 0.98913114824820624893, 
1.00189541647120616830, 1.00189541647120616830, 
1.01037178591362075686, 1.01037178591362075686
};
//stop-lepto-0405e-hadronic-scale-correction-cells//

//start-ariadne-0405e-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:55 CET 2010 
static Double_t ariadne_factor_cells_0405e[number_of_eta_bins] = {
  1.02842505932429340021, 1.02842505932429340021, 1.02842505932429340021, 1.02842505932429340021, 
  1.00943501005649993729, 1.00943501005649993729, 1.00943501005649993729, 1.00943501005649993729, 
  0.98584902144858788109, 0.98584902144858788109, 
  1.00063178022954613056, 1.00063178022954613056, 
  1.01221602560694856621, 1.01221602560694856621  
};
//stop-ariadne-0405e-hadronic-scale-correction-cells//

//start-lepto-06e-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:58 CET 2010 
static Double_t lepto_factor_cells_06e[number_of_eta_bins] = {
1.03385279669320762075, 1.03385279669320762075, 1.03385279669320762075, 1.03385279669320762075, 
1.02126960986365378226, 1.02126960986365378226, 1.02126960986365378226, 1.02126960986365378226, 
0.99579737032633519611, 0.99579737032633519611, 
1.01109878925485219270, 1.01109878925485219270, 
1.01254325501886000005, 1.01254325501886000005
};
//stop-lepto-06e-hadronic-scale-correction-cells//

//start-ariadne-06e-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:57 CET 2010 
static Double_t ariadne_factor_cells_06e[number_of_eta_bins] = {
  1.02971973304207287114, 1.02971973304207287114, 1.02971973304207287114, 1.02971973304207287114, 
  1.01471208348225094831, 1.01471208348225094831, 1.01471208348225094831, 1.01471208348225094831, 
  0.99325067892659268942, 0.99325067892659268942, 
  1.01082830314668781213, 1.01082830314668781213, 
  1.01484917707574595624, 1.01484917707574595624
};
//stop-ariadne-06e-hadronic-scale-correction-cells//

//start-lepto-0607p-hadronic-scale-correction-cells//
// Sun Feb 14 22:04:00 CET 2010 
static Double_t lepto_factor_cells_0607p[number_of_eta_bins] = {
1.03346749369724766332, 1.03346749369724766332, 1.03346749369724766332, 1.03346749369724766332, 
1.01761860683906180114, 1.01761860683906180114, 1.01761860683906180114, 1.01761860683906180114, 
0.99321522163233877389, 0.99321522163233877389, 
1.00857514933604752727, 1.00857514933604752727, 
1.01123947400645608496, 1.01123947400645608496
};
//stop-lepto-0607p-hadronic-scale-correction-cells//

//start-ariadne-0607p-hadronic-scale-correction-cells//
// Sun Feb 14 22:03:59 CET 2010 
static Double_t ariadne_factor_cells_0607p[number_of_eta_bins] = {
1.03056105869332648162, 1.03056105869332648162, 1.03056105869332648162, 1.03056105869332648162, 
1.01312221577016092766, 1.01312221577016092766, 1.01312221577016092766, 1.01312221577016092766, 
0.99049546039429625655, 0.99049546039429625655, 
1.00705953481823806683, 1.00705953481823806683, 
1.01377495819234697549, 1.01377495819234697549 
};
//stop-ariadne-0607p-hadronic-scale-correction-cells//

static Double_t lepto_test_factor_cells_0607p[number_of_eta_bins] = { 
1.03378489028789921633, 1.03378489028789921633, 1.03378489028789921633, 1.03378489028789921633,
1.01914023300842293551, 1.01914023300842293551, 1.01914023300842293551, 1.01914023300842293551, 
0.99257106435677633360, 0.99257106435677633360, 
1.00794562280882149530, 1.00794562280882149530, 
1.01132103313399390743, 1.01132103313399390743
};

static Double_t ariadne_test_factor_cells_0607p[number_of_eta_bins] = { 
1.03378489028789921633, 1.03378489028789921633, 1.03378489028789921633, 1.03378489028789921633,
1.01914023300842293551, 1.01914023300842293551, 1.01914023300842293551, 1.01914023300842293551, 
0.99257106435677633360, 0.99257106435677633360, 
1.00794562280882149530, 1.00794562280882149530, 
1.01132103313399390743, 1.01132103313399390743
};

//start-lepto-989900-jet-energy-correction-cells//
// Mon Feb 15 11:30:27 CET 2010 
 static Double_t lepto_989900_cells[number_of_eta_bins][2] = {
{0.1849227963489093184978884, 0.7453767911160920078827985},
{-0.5040924294828398322110274, 0.8424526325756606093264622},
{-0.3981631171559820714378475, 0.8704850393574046218958529},
{-0.5218927020286086015943283, 0.8966699062901966055250114},
{-0.4989186136866772813824866, 0.8897986338151886442915384},
{-0.4636338031908088797017342, 0.8817252393722254666386107},
{-0.5886855960563245337979765, 0.8756298567376629149450196},
{-0.6829619839369843736776033, 0.8564877898490527163843922},
{-0.6102544836896163227279999, 0.7978548624799891575420929},
{-0.4919653747498594076326128, 0.8138192227403022727116877},
{-0.6497798525941992364707289, 0.8796228650415203187407087},
{-0.5204317492262586775098043, 0.9086371764882766877136078},
{-0.2135783938648975199825486, 0.9213911941146836204197257},
{0.4976661719835241415310634, 0.8981605843195622496111241}
};
//stop-lepto-989900-jet-energy-correction-cells//

//start-ariadne-989900-jet-energy-correction-cells//
// Mon Feb 15 11:30:22 CET 2010 
static Double_t ariadne_989900_cells[number_of_eta_bins][2] = {
{0.2676319721655046612873718, 0.7375115997465190309512195},
{-0.0764228623638322501765430, 0.7986109758473056707828164},
{-0.2051523253503808874587122, 0.8463196579987710688186553},
{-0.4052527387800591296418418, 0.8887852263203325220430884},
{-0.2937397452422775834257607, 0.8717228990726222370710730},
{-0.4719187706069959964416682, 0.8847704594343392958677441},
{-0.4246789592523220191466748, 0.8677806732507772347062769},
{-0.5500369921008881890500675, 0.8460266908841060917723098},
{-0.4852657970803974651907708, 0.7899836121976108183062593},
{-0.5135606701874362300941357, 0.8212965858127055929216453},
{-0.4956329129299711366840597, 0.8747844053657789498146258},
{-0.4447603909602636518272334, 0.9058135667404716429018663},
{-0.2279378735484343110329064, 0.9356753186581331371840520},
{0.1009740184152077557833138, 0.9834524310731996843415459}
};
//stop-ariadne-989900-jet-energy-correction-cells//

//start-lepto-0405e-jet-energy-correction-cells//
// Mon Feb 15 11:30:36 CET 2010 
 static Double_t lepto_0405e_cells[number_of_eta_bins][2] = {
{-0.1416814587562369165407716, 0.7586663483944117158941367},
{-0.1870919086190597213636977, 0.7814608583172878431000186},
{-0.2632766781879706807778518, 0.8290868213886158910597146},
{-0.3369666779895192920690761, 0.8619045851965633220359564},
{-0.3255572602664155157370374, 0.8602166718671052647593456},
{-0.4198513764858714525374239, 0.8717385182575667634807814},
{-0.4876778069683951888357853, 0.8669751459886041455860095},
{-0.6180187765911269170970854, 0.8565011527823898074274211},
{-0.5035736498819244344460344, 0.8050768908003671509021615},
{-0.5336730318087457014541997, 0.8381717949984434090637819},
{-0.5222782687171805271475478, 0.8865718465289715854993347},
{-0.3724291397641172896015860, 0.8996898225716721819722466},
{-0.1381130247517866338480275, 0.9275473896671319051776550},
{0.3333283642471459029543723, 0.9276179781344889851979474}
};
//stop-lepto-0405e-jet-energy-correction-cells//

//start-ariadne-0405e-jet-energy-correction-cells//
// Mon Feb 15 11:30:32 CET 2010 
static Double_t ariadne_0405e_cells[number_of_eta_bins][2] = {
{0.2225373811557372905056695, 0.7176739128227159758566245},
{-0.1652858825855437197027697, 0.7803729139662227698437391},
{-0.1497961052909715928826984, 0.8219978028053079865955510},
{-0.2854088263627053145832235, 0.8600148328514193307725577},
{-0.2583301831761816158561373, 0.8528889284722284136464054},
{-0.2843580432674828739614270, 0.8582787104801181676805300},
{-0.4090196748506415569401895, 0.8643252201904678688748618},
{-0.4938678300957760614053882, 0.8487714637885368507141948},
{-0.3813281762360967164227077, 0.7993846889444983228045771},
{-0.4155872485521025505583737, 0.8382697960595822950224942},
{-0.3467708894468328750804176, 0.8815507972845165740594098},
{-0.2311269185819042337737272, 0.9043775075920014172581318},
{-0.2478506317557906402093693, 0.9541451684030349822407402},
{0.3159286453417153306411080, 0.9580548534460318066408036}
};
//stop-ariadne-0405e-jet-energy-correction-cells//
 
//start-lepto-06e-jet-energy-correction-cells//
// Mon Feb 15 11:30:45 CET 2010 
 static Double_t lepto_06e_cells[number_of_eta_bins][2] = {
{0.1520079026669963195850244, 0.7184518809381781334622019},
{-0.1612629801835912535246109, 0.7831654138601723547097322},
{-0.2880410195074482593291521, 0.8332324152093121272244503},
{-0.3592586763345858025253676, 0.8711500389664946775170051},
{-0.4489119472646265207593785, 0.8769952991757498095282131},
{-0.3540991775998418300552828, 0.8676444157606099105706221},
{-0.5309540609179738135381399, 0.8737880824948036639909787},
{-0.6170447785516101069092088, 0.8597907420455488392008192},
{-0.3212932545751582824600234, 0.7967897268131414190150963},
{-0.5833355705783670464015245, 0.8499748195721087418164075},
{-0.5342315711713360171231102, 0.8969242723487351742761575},
{-0.4737948670496093250825709, 0.9245480570466966252496377},
{-0.1855454743060036548740754, 0.9370340304649934504155340},
{0.1906237482990523757209189, 0.9400853379537014165023834}
};
//stop-lepto-06e-jet-energy-correction-cells//

//start-ariadne-06e-jet-energy-correction-cells//
// Mon Feb 15 11:30:40 CET 2010 
static Double_t ariadne_06e_cells[number_of_eta_bins][2] = {
{0.3682730309188352868510208, 0.6951779214055960887108654},
{-0.2963365468254468004793978, 0.7922009683743140850964437},
{-0.3804472870634948433732347, 0.8526724828677524348208294},
{-0.2453025337613712886231099, 0.8589116107190382098934833},
{-0.3869890250267897235936232, 0.8708463045755555631544098},
{-0.3672895195308965576863613, 0.8719829634317454170044925},
{-0.3667271994323068606824734, 0.8637687628711911980161631},
{-0.4788065720751108478303593, 0.8520635302034017533046040},
{-0.4512892374983840104718524, 0.8086323731777917966923042},
{-0.5841506827300547710635215, 0.8593914264620228227187226},
{-0.5867108134792367035714733, 0.9140493947391423734671889},
{-0.3426638104516264715648788, 0.9169729602826580672214618},
{-0.1246472427925167381079774, 0.9436043318302111959638978},
{-0.0690347446113764862118956, 1.0013264563312345156731453}
};
//stop-ariadne-06e-jet-energy-correction-cells//

//start-lepto-0607p-jet-energy-correction-cells//
// Mon Feb 15 11:30:53 CET 2010 
 static Double_t lepto_0607p_cells[number_of_eta_bins][2] = {
{0.2697198541243220337726427, 0.7057277033628084339511588},
{-0.0969184826401890547886708, 0.7776475001082610027935971},
{-0.3128709222930441957721825, 0.8384296391249159707825811},
{-0.2467051034708179491072855, 0.8546859264280391998980235},
{-0.4165929996062189277772347, 0.8727285806726619021844726},
{-0.4484019108519092244868887, 0.8759896666769411144315427},
{-0.5729007335738222073118209, 0.8795821846030360546819793},
{-0.6282623401311580035866200, 0.8596465394238299584372953},
{-0.5337504894535839428115764, 0.8103862337445485453102378},
{-0.6822658330887703126421684, 0.8575353416995695754465601},
{-0.5029139650181817522422989, 0.8909381119353076972089411},
{-0.3377897481995944528598841, 0.9059693703212764459564710},
{-0.1046060011474696654909522, 0.9311529070130365459689870},
{0.1854596691169398403609136, 0.9358370257136565673050654}
};
//stop-lepto-0607p-jet-energy-correction-cells//

      //start-ariadne-0607p-jet-energy-correction-cells//
// Mon Feb 15 11:30:49 CET 2010 
 static Double_t ariadne_0607p_cells[number_of_eta_bins][2] = {
{-0.0858226696850253600867831, 0.7548814058667212556485993},
{-0.1531167487445077857177722, 0.7839156877907866860510921},
{-0.2667503320891648765922355, 0.8352685411147485750049668},
{-0.2632198199304686681010423, 0.8614428601881457536748599},
{-0.2977937064161247371352204, 0.8633573667573156695809189},
{-0.3775927057560864441576598, 0.8739232043299010044279385},
{-0.4092680240713346129233230, 0.8681320098926734729616328},
{-0.4578000203116136535008707, 0.8496187314212856334094681},
{-0.2850872974766703871196682, 0.7956851678217694567862850},
{-0.4307182636647894358183919, 0.8443605043108689534037126},
{-0.4088968524429608675241354, 0.8924701704004568414418941},
{-0.1824491115758038739524949, 0.9072424715069741241180168},
{-0.1567117610620280543365368, 0.9558452476745694648130325},
{0.2861932233284934512518305, 0.9644215895459269383138690}
};
//stop-ariadne-0607p-jet-energy-correction-cells//


//2nd analysis correction factors

static Double_t C_scale_lepto_989900[number_of_eta_bins] = {
  1.0225028922, 1.0225028922, 1.0225028922, 1.0225028922, 
  1.0127230592, 1.0127230592, 1.0127230592, 1.0127230592, 
  0.9860663260, 0.9860663260, 
  1.0086231693, 1.0086231693, 
  1.0077253061, 1.0077253061
};

static Double_t C_scale_ariadne_989900[number_of_eta_bins] = {
  1.0198914384, 1.0198914384, 1.0198914384, 1.0198914384,
  1.0074458532, 1.0074458532, 1.0074458532, 1.0074458532,
  0.9830311677, 0.9830311677,
  1.0072909294, 1.0072909294,
  1.0093458301, 1.0093458301
};

static Double_t C_scale_lepto_040506e[number_of_eta_bins] = {
  1.0331997763, 1.0331997763, 1.0331997763, 1.0331997763,
  1.0144512448, 1.0144512448, 1.0144512448, 1.0144512448,
  0.9884928323, 0.9884928323,
  1.0016552918, 1.0016552918,
  1.0087696414, 1.0087696414
};

static Double_t C_scale_ariadne_040506e[number_of_eta_bins] = {
  1.0331997763, 1.0331997763, 1.0331997763, 1.0331997763,
  1.0144512448, 1.0144512448, 1.0144512448, 1.0144512448,
  0.9884928323, 0.9884928323,
  1.0016552918, 1.0016552918,
  1.0087696414, 1.0087696414
};

static Double_t C_scale_lepto_0405e[number_of_eta_bins] = {
  1.0330104247, 1.0330104247, 1.0330104247, 1.0330104247, 
  1.0135996993, 1.0135996993, 1.0135996993, 1.0135996993, 
  0.9892943495, 0.9892943495, 
  1.0018700972, 1.0018700972, 
  1.0103521273, 1.0103521273
};

static Double_t C_scale_ariadne_0405e[number_of_eta_bins] = {
  1.0283467706, 1.0283467706, 1.0283467706, 1.0283467706,
  1.0095403507, 1.0095403507, 1.0095403507, 1.0095403507,
  0.9860448972, 0.9860448972,
  1.0005940244, 1.0005940244,
  1.0121803795, 1.0121803795
};

static Double_t C_scale_lepto_06e[number_of_eta_bins] = {
  1.0336222652, 1.0336222652, 1.0336222652, 1.0336222652, 
  1.0210305567, 1.0210305567, 1.0210305567, 1.0210305567, 
  0.9958293088, 0.9958293088, 
  1.0114150622, 1.0114150622,
  1.0127002674, 1.0127002674
};

static Double_t C_scale_ariadne_06e[number_of_eta_bins] = {
  1.0296597873, 1.0296597873, 1.0296597873, 1.0296597873,
  1.0143370120, 1.0143370120, 1.0143370120, 1.0143370120,
  0.9933600973, 0.9933600973,
  1.0111724002, 1.0111724002,
  1.0150863752, 1.0150863752
};

static Double_t C_scale_lepto_0607p[number_of_eta_bins] = {
  1.0336436865, 1.0336436865, 1.0336436865, 1.0336436865, 
  1.0178750245, 1.0178750245, 1.0178750245, 1.0178750245, 
  0.9933788886, 0.9933788886, 
  1.0086423855, 1.0086423855, 
  1.0119956772, 1.0119956772
};

static Double_t C_scale_ariadne_0607p[number_of_eta_bins] = {
  1.0306111292, 1.0306111292, 1.0306111292, 1.0306111292,
  1.0134667638, 1.0134667638, 1.0134667638, 1.0134667638,
  0.9906452365, 0.9906452365,
  1.0070537995, 1.0070537995,
  1.0144976676, 1.0144976676
};
static Double_t C_scale_lepto_test_0607p[number_of_eta_bins] = {
  1.0348378534, 1.0348378534, 1.0348378534, 1.0348378534, 
  1.0188460231, 1.0188460231, 1.0188460231, 1.0188460231, 
  0.9930495372, 0.9930495372, 
  1.0086470519, 1.0086470519, 
  1.0116133431, 1.0116133431
};
static Double_t C_scale_ariadne_test_0607p[number_of_eta_bins] = {
  1.0307463325, 1.0307463325, 1.0307463325, 1.0307463325,
  1.0138564216, 1.0138564216, 1.0138564216, 1.0138564216, 
  0.9902727731, 0.9902727731, 
  1.0071489727, 1.0071489727, 
  1.0133274028, 1.0133274028
};
static Double_t C_scale_lepto_0607pr[number_of_eta_bins] = {
  1.0252694383, 1.0252694383, 1.0252694383, 1.0252694383,
  1.0165938565, 1.0165938565, 1.0165938565, 1.0165938565, 
  0.9862916471, 0.9862916471, 
  1.0054365385, 1.0054365385, 
  1.0085375312, 1.0085375312
};

static Double_t C_scale_lepto_0607pl[number_of_eta_bins] = {
  1.0266258580, 1.0266258580, 1.0266258580, 1.0266258580,
  1.0180222977, 1.0180222977, 1.0180222977, 1.0180222977, 
  0.9888860710, 0.9888860710, 
  1.0052515647, 1.0052515647, 
  1.0100809525, 1.0100809525
};

static Double_t C_scale_no_corr[number_of_eta_bins] = {
  1.,  1., 1., 1.,
  1., 1., 1., 1.,
  1., 1.,
  1., 1.,
  1., 1.
};

static Double_t a0_corr_lepto_989900[number_of_eta_bins] = {
  0.371478, -0.544991, -0.437083, -0.561746, -0.423271, -0.478446, -0.503835, -0.546216, -0.495498, -0.405129, -0.429863, -0.399184, -0.174602, 0.146424
};
static Double_t a1_corr_lepto_989900[number_of_eta_bins] = {
  0.717379, 0.848159, 0.873862, 0.90077, 0.882385, 0.885539, 0.868592, 0.842674, 0.787146, 0.806415, 0.859821, 0.900416, 0.929918, 0.948647
};

static Double_t a0_corr_ariadne_989900[number_of_eta_bins] = {
0.274404, -0.0131425, -0.238328, -0.372403, -0.214161, -0.445541, -0.31771, -0.485488, -0.377258, -0.389489, -0.369316, -0.301869, -0.165053, 0.5061
};
static Double_t a1_corr_ariadne_989900[number_of_eta_bins] = {
0.730019, 0.787965, 0.851243, 0.887464, 0.864729, 0.885984, 0.858542, 0.841696, 0.780604, 0.813894, 0.86767, 0.896887, 0.935331, 0.942634
};

static Double_t a0_corr_lepto_040506e[number_of_eta_bins] = {-0.900762, -0.027889, -0.119307, -0.264639, -0.338724, -0.327753, -0.36992, -0.452522, -0.339662, -0.432949, -0.404978, -0.448368, -0.235449, 0.0701483};
static Double_t a1_corr_lepto_040506e[number_of_eta_bins] = {0.860666, 0.760438, 0.813597, 0.853678, 0.866559, 0.863476, 0.857733, 0.841143, 0.788444, 0.834551, 0.882007, 0.921459, 0.944337, 0.95925};

static Double_t a0_corr_ariadne_040506e[number_of_eta_bins] = {-0.900762, -0.027889, -0.119307, -0.264639, -0.338724, -0.327753, -0.36992, -0.452522, -0.339662, -0.432949, -0.404978, -0.448368, -0.235449, 0.0701483};
static Double_t a1_corr_ariadne_040506e[number_of_eta_bins] = {0.860666, 0.760438, 0.813597, 0.853678, 0.866559, 0.863476, 0.857733, 0.841143, 0.788444, 0.834551, 0.882007, 0.921459, 0.944337, 0.95925};

static Double_t a0_corr_lepto_0405e[number_of_eta_bins] = {
  -0.0115378, 0.0217799, -0.342834, -0.342917, -0.344993, -0.341084, -0.430182, -0.411263, -0.299054, -0.413295, -0.303871, -0.263272, -0.130439, 0.316174
};
static Double_t a1_corr_lepto_0405e[number_of_eta_bins] = {
  0.736107, 0.752704, 0.838865, 0.863251, 0.863812, 0.862795, 0.862373, 0.836555, 0.786412, 0.828142, 0.864368, 0.895586, 0.932258, 0.938019
};

static Double_t a0_corr_ariadne_0405e[number_of_eta_bins] = {
0.309831, -0.812089, -0.1367, -0.297971, -0.255224, -0.314823, -0.36343, -0.389898, -0.209286, -0.294967, -0.205577, -0.144427, -0.185809, 0.256132
};
static Double_t a1_corr_ariadne_0405e[number_of_eta_bins] = {
0.702584, 0.860734, 0.825322, 0.866615, 0.857402, 0.867267, 0.866837, 0.844046, 0.785107, 0.832084, 0.870849, 0.900414, 0.948102, 0.960532
};

static Double_t a0_corr_lepto_06e[number_of_eta_bins] = {
  -0.147735, -0.0439031, -0.298287, -0.361772, -0.448146, -0.283049, -0.574606, -0.505626, -0.221189, -0.558173, -0.19498, -0.422563, -0.0247049, 0.0838887,
};
static Double_t a1_corr_lepto_06e[number_of_eta_bins] = {
  0.753809, 0.763386, 0.834941, 0.871255, 0.875909, 0.859106, 0.88134, 0.847222, 0.786637, 0.844492, 0.865203, 0.927898, 0.92715, 0.961016
};

static Double_t a0_corr_ariadne_06e[number_of_eta_bins] = {
0.208482, -0.334897, -0.201592, -0.206557, -0.388075, -0.420187, -0.425316, -0.394716, -0.231174, -0.212338, -0.224531, -0.278957, 0.0729052, 0.258603,
};
static Double_t a1_corr_ariadne_06e[number_of_eta_bins] = {
0.71033, 0.800384, 0.832528, 0.85867, 0.880793, 0.885185, 0.877943, 0.851621, 0.788882, 0.825817, 0.879192, 0.916129, 0.92475, 0.965931
};

static Double_t a0_corr_lepto_0607p[number_of_eta_bins] = {
  0.157304, 0.0560584, -0.21093, -0.274096, -0.338353, -0.404981, -0.454319, -0.563538, -0.330781, -0.45229, -0.390409, -0.214802, -0.0354642, 0.197917,
};
static Double_t a1_corr_lepto_0607p[number_of_eta_bins] = {
  0.715862, 0.754733, 0.823202, 0.858535, 0.864347, 0.872465, 0.868016, 0.853751, 0.787941, 0.834452, 0.884165, 0.897311, 0.928416, 0.941718
};

static Double_t a0_corr_ariadne_0607p[number_of_eta_bins] = {
0.0342872, 0.0170076, -0.20871, -0.233043, -0.313787, -0.330107, -0.347063, -0.350122, -0.117197, -0.25752, -0.318401, 0.0173798, -0.0518676, 0.298698,
};
static Double_t a1_corr_ariadne_0607p[number_of_eta_bins] = {
0.736535, 0.761035, 0.829381, 0.862394, 0.871721, 0.874463, 0.868625, 0.842114, 0.779895, 0.8315, 0.888192, 0.890801, 0.943441, 0.958188
};

static Double_t a0_corr_lepto_test_0607p[number_of_eta_bins] = {
  0.13086, -0.0875975, -0.0704085, -0.21607, -0.255768, -0.29617, -0.373893, -0.550067, -0.281231, -0.466699, -0.371816, -0.237988, -0.185237, 0.1334
};
static Double_t a1_corr_lepto_test_0607p[number_of_eta_bins] = {
  0.718119, 0.77003, 0.80946, 0.853998, 0.857412, 0.862881, 0.861148, 0.852516, 0.783407, 0.834636, 0.880465, 0.897784, 0.939777, 0.94244
};

static Double_t a0_corr_ariadne_test_0607p[number_of_eta_bins] = {
0.121433, -0.0598809, -0.128602, -0.171875, -0.225122, -0.218616, -0.249222, -0.327882, -0.193571, -0.377589, -0.35448, -0.175048, -0.221599, 0.17493
};

static Double_t a1_corr_ariadne_test_0607p[number_of_eta_bins] = {
  0.722044, 0.766442, 0.819193, 0.852875, 0.857674, 0.860496, 0.854396, 0.83673, 0.782634, 0.838028, 0.891386, 0.90799, 0.959113, 0.973407
};


static Double_t a0_no_corr[number_of_eta_bins] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
static Double_t a1_no_corr[number_of_eta_bins] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

/*Double_t        C_scale[number_of_eta_bins];
  Double_t        a0_corr[number_of_eta_bins];
  Double_t        a1_corr[number_of_eta_bins];*/


#endif

