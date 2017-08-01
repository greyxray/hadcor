#ifndef  AFG_h
#define  AFG_h
#include "../inc/constants.h"

Double_t * all_bins_arrays[12] = {et_bin, eta_bin_crosssec, Q2_bin, x_bin, et_jet_bin, eta_jet_bin, xgamma_bin, xp_bin, dphi_bin, deta_bin, dphi_e_ph_bin, deta_e_ph_bin, };

//Fontannaz
	// 0 < Q2 < 350 range
		//NEW
			//pT_cut in center-of-mass frame = 0.5 GeV/c
				Double_t font_xgamma_pt05_Q2full[number_xgamma_bins] = {1.27, 3.04, 5.45, 8.86, 18.54, 49.39};
				Double_t font_xp_pt05_Q2full[number_xp_bins] = 	{380., 761., 446., 152., 29.13, 2.90};
				//Double_t font_dphi_pt05_Q2full[number_dphi_bins] = 
				Double_t font_deta_pt05_Q2full[number_deta_bins] = {0.884, 3.457, 3.743, 3.318, 2.180, 0.435};
				Double_t font_dphi_e_gamma_pt05_Q2full[number_dphi_e_ph_bins] = {0.0185, 0.024, 0.037, 0.064, 0.106, 0.133};
				//Double_t font_deta_e_gamma_pt05_Q2full[number_deta_e_ph_bins] =

				//e = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; print [ format(l[i] * (0.08), '.5f') for i in range(len(l)) ]
					Double_t font_xgamma_pos_pt05_Q2full[number_xgamma_bins] = {0.10160, 0.24320, 0.43600, 0.70880, 1.48320, 3.95120};
					Double_t font_xp_pos_pt05_Q2full[number_xp_bins] = {30.40000, 60.88000, 35.68000, 12.16000, 2.33040, 0.23200};
					Double_t font_deta_pos_pt05_Q2full[number_deta_bins] = {0.07072, 0.27656, 0.29944, 0.26544, 0.17440, 0.03480};
					Double_t font_dphi_e_gamma_pos_pt05_Q2full[number_dphi_e_ph_bins] = {0.00148, 0.00192, 0.00296, 0.00512, 0.00848, 0.01064};
						Double_t font_xgamma_neg_pt05_Q2full[number_xgamma_bins];
						Double_t font_xp_neg_pt05_Q2full[number_xp_bins];
						Double_t font_deta_neg_pt05_Q2full[number_dphi_bins];
						Double_t font_dphi_e_gamma_neg_pt05_Q2full[number_deta_bins];
			//pT_cut in center-of-mass frame = 2.5 GeV/c - ADDITIONAL 8% unc added everywhere - for now: only this 10 % included!!!!
				Double_t font_xgamma_pt25_Q2full[number_xgamma_bins] = {1.2178, 2.8124, 4.7374, 7.3692, 14.5175, 45.1021};//l=[1.2178, 2.8124, 4.7374, 7.3692, 14.5175, 45.1021]
				Double_t font_xp_pt25_Q2full[number_xp_bins] = {326.7, 605.4, 344.1, 130.3, 26.76, 2.865};//l=[326.7, 605.4, 344.1, 130.3, 26.76, 2.865]
				Double_t font_dphi_pt25_Q2full[number_dphi_bins] = {0.0181, 0.0546, 0.0720, 0.0829, 0.0951, 0.1006, 0.0778};// Before the July2017 email {0.0118, 0.0455, 0.0639, 0.0711, 0.0849, 0.113, 0.188}; // l = [0.0118, 0.0455, 0.0639, 0.0711, 0.0849, 0.113, 0.188]
				Double_t font_deta_pt25_Q2full[number_deta_bins] = {0.593, 2.121, 3.167, 3.116, 2.095, 0.461};//l=[0.593, 2.121, 3.167, 3.116, 2.095, 0.461]
				Double_t font_dphi_e_gamma_pt25_Q2full[number_dphi_e_ph_bins] = {0.0180, 0.0235, 0.0363, 0.0628, 0.0936, 0.0798};//l=[0.0180, 0.0235, 0.0363, 0.0628, 0.0936, 0.0798]
				Double_t font_deta_e_gamma_pt25_Q2full[number_deta_e_ph_bins] = {1.451, 3.545, 4.476, 3.019, 1.070}; // l = [1.451, 3.545, 4.476, 3.019, 1.070]

				//e = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; print [ format(l[i] * (0.08), '.5f') for i in range(len(l)) ]
					Double_t font_xgamma_pos_pt25_Q2full[number_xgamma_bins] = {0.09742, 0.22499, 0.37899, 0.58954, 1.16140, 3.60817};
					Double_t font_xp_pos_pt25_Q2full[number_xp_bins] = {26.13600, 48.43200, 27.52800, 10.42400, 2.14080, 0.22920};
					Double_t font_dphi_pos_pt25_Q2full[number_dphi_bins] ={0.00094, 0.00364, 0.00511, 0.00569, 0.00679, 0.00904, 0.01504};
					Double_t font_deta_pos_pt25_Q2full[number_deta_bins] = {0.04744, 0.16968, 0.25336, 0.24928, 0.16760, 0.03688};
					Double_t font_dphi_e_gamma_pos_pt25_Q2full[number_dphi_e_ph_bins] = {0.00144, 0.00188, 0.00290, 0.00502, 0.00749, 0.00638};
					Double_t font_deta_e_gamma_pos_pt25_Q2full[number_deta_e_ph_bins] = {0.11608, 0.28360, 0.35808, 0.24152, 0.08560};
						Double_t font_xgamma_neg_pt25_Q2full[number_xgamma_bins];
						Double_t font_xp_neg_pt25_Q2full[number_xp_bins];
						Double_t font_deta_neg_pt25_Q2full[number_dphi_bins];
						Double_t font_dphi_e_gamma_neg_pt25_Q2full[number_deta_bins];
						Double_t font_dphi_neg_pt25_Q2full[number_dphi_e_ph_bins];
						Double_t font_deta_e_gamma_neg_pt25_Q2full[number_deta_e_ph_bins];		

	// 0 < Q2 < 30 range
		
		//NEW
			//pT_cut in  in center-of-mass frame = 0.5 GeV/c
				// l1  = [[0.64, 1.85, 2.51, 3.47, 6.33, 20.60], [171.62, 306.68, 157.94, 52.04, 10.13, 0.861], [0.342, 1.070, 1.654, 1.464, 0.897, 0.252], [0.0106, 0.0134, 0.0195, 0.0291, 0.0356, 0.0367]]
				Double_t font_xgamma_pt05_Q2lt30[number_xgamma_bins] = {0.64, 1.85, 2.51, 3.47, 6.33, 20.60}; // l = [0.64, 1.85, 2.51, 3.47, 6.33, 20.60]
				Double_t font_xp_pt05_Q2lt30[number_xp_bins] = {171.62, 306.68, 157.94, 52.04, 10.13, 0.963}; // l = [171.62, 306.68, 157.94, 52.04, 10.13, 0.861]
				Double_t font_deta_pt05_Q2lt30[number_deta_bins] = {0.342, 1.070, 1.654, 1.464, 0.897, 0.252}; // l = [0.342, 1.070, 1.654, 1.464, 0.897, 0.252]
				Double_t font_dphi_e_gamma_pt05_Q2lt30[number_dphi_e_ph_bins] = {0.0106, 0.0134, 0.0195, 0.0291, 0.0356, 0.0367}; // l = [0.0106, 0.0134, 0.0195, 0.0291, 0.0356, 0.0367]

				//e = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; print [ format(l[i] * (0.08), '.5f') for i in range(len(l)) ]
				//dl1 = [[0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]]
					Double_t font_xgamma_pos_pt05_Q2lt30[number_xgamma_bins] = {0.05120, 0.14800, 0.20080, 0.27760, 0.50640, 1.64800};
					Double_t font_xp_pos_pt05_Q2lt30[number_xp_bins] = {13.72960, 24.53440, 12.63520, 4.16320, 0.81040, 0.06888};
					Double_t font_deta_pos_pt05_Q2lt30[number_deta_bins] = {0.02736, 0.08560, 0.13232, 0.11712, 0.07176, 0.02016};
					Double_t font_dphi_e_gamma_pos_pt05_Q2lt30[number_dphi_e_ph_bins] = {0.00085, 0.00107, 0.00156, 0.00233, 0.00285, 0.00294};
						Double_t font_xgamma_neg_pt05_Q2lt30[number_xgamma_bins];
						Double_t font_xp_neg_pt05_Q2lt30[number_xp_bins];
						Double_t font_deta_neg_pt05_Q2lt30[number_dphi_bins];
						Double_t font_dphi_e_gamma_neg_pt05_Q2lt30[number_deta_bins];
			//pT_cut in center-of-mass frame = 2.5 GeV/c
				// l2 =[[0.613, 1.80, 2.43, 3.26, 5.68, 18.78], [158.89, 281.51, 146.03, 49.92, 9.6153, 0.899], [0.296, 0.921, 1.419, 1.396, 0.871, 0.259], [0.0105, 0.0132, 0.0190, 0.0286, 0.0327, 0.0302]]
				Double_t font_xgamma_pt25_Q2lt30[number_xgamma_bins] = {0.613, 1.80, 2.43, 3.26, 5.68, 18.78}; // l = [0.613, 1.80, 2.43, 3.26, 5.68, 18.78]
				Double_t font_xp_pt25_Q2lt30[number_xp_bins] = {158.89, 281.51, 146.03, 49.92, 9.6153, 0.899}; // l = [158.89, 281.51, 146.03, 49.92, 9.6153, 0.899]
				Double_t font_dphi_pt25_Q2lt30[number_dphi_bins] = {0.00265, 0.0208, 0.0389, 0.0476, 0.0572, 0.0659, 0.0515}; // l = [0.00265, 0.0208, 0.0389, 0.0476, 0.0572, 0.0659, 0.0515]
				Double_t font_deta_pt25_Q2lt30[number_deta_bins] = {0.296, 0.921, 1.419, 1.396, 0.871, 0.259}; // l = [0.296, 0.921, 1.419, 1.396, 0.871, 0.259]
				Double_t font_dphi_e_gamma_pt25_Q2lt30[number_dphi_e_ph_bins] = {0.0105, 0.0132, 0.0190, 0.0286, 0.0327, 0.0302}; // l = [0.0105, 0.0132, 0.0190, 0.0286, 0.0327, 0.0302]
				Double_t font_deta_e_gamma_pt25_Q2lt30[number_deta_e_ph_bins] = {1.408, 2.236, 1.875, 0.512, 0.017}; // l = [1.408, 2.236, 1.875, 0.512, 0.017]

				//e = [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]; print [ format(l[i] * (0.08), '.5f') for i in range(len(l)) ]
				//dl2 = [[0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08], [0.08, 0.08, 0.08, 0.08, 0.08, 0.08]]
					Double_t font_xgamma_pos_pt25_Q2lt30[number_xgamma_bins] = {0.04904, 0.14400, 0.19440, 0.26080, 0.45440, 1.50240};
					Double_t font_xp_pos_pt25_Q2lt30[number_xp_bins] = {12.71120, 22.52080, 11.68240, 3.99360, 0.76922, 0.07192};
					Double_t font_dphi_pos_pt25_Q2lt30[number_dphi_bins] = {0.00021, 0.00166, 0.00311, 0.00381, 0.00458, 0.00527, 0.00412};
					Double_t font_deta_pos_pt25_Q2lt30[number_deta_bins] = {0.02368, 0.07368, 0.11352, 0.11168, 0.06968, 0.02072};
					Double_t font_dphi_e_gamma_pos_pt25_Q2lt30[number_dphi_e_ph_bins] = {0.00084, 0.00106, 0.00152, 0.00229, 0.00262, 0.00242};
					Double_t font_deta_pos_e_gamma_pt25_Q2lt30[number_deta_e_ph_bins] = {0.11264, 0.17888, 0.15000, 0.04096, 0.00136};
						Double_t font_xgamma_neg_pt25_Q2lt30[number_xgamma_bins];
						Double_t font_xp_neg_pt25_Q2lt30[number_xp_bins];
						Double_t font_dphi_neg_pt25_Q2lt30[number_dphi_bins];
						Double_t font_deta_neg_pt25_Q2lt30[number_dphi_bins];
						Double_t font_dphi_e_gamma_neg_pt25_Q2lt30[number_deta_bins];	
						Double_t font_deta_neg_e_gamma_pt25_Q2lt30[number_deta_e_ph_bins];
			
	// 30 < Q2 < 350 range - EMPTY and is filled later
			//pT_cut in  in center-of-mass frame = 0.5 GeV/c
				Double_t font_xgamma_pt05_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_pt05_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_pt05_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_pt05_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};

				Double_t font_xgamma_pos_pt05_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_pos_pt05_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_pos_pt05_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_pos_pt05_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};

				Double_t font_xgamma_neg_pt05_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_neg_pt05_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_neg_pt05_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_neg_pt05_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};
			//pT_cut in center-of-mass frame = 2.5 GeV/c
				Double_t font_xgamma_pt25_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_pt25_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_pt25_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_pt25_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};

				Double_t font_xgamma_pos_pt25_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_pos_pt25_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_pos_pt25_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_pos_pt25_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};

				Double_t font_xgamma_neg_pt25_Q2gt30[number_xgamma_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_xp_neg_pt25_Q2gt30[number_xp_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_deta_neg_pt25_Q2gt30[number_deta_bins] = {0, 0, 0, 0, 0, 0};
				Double_t font_dphi_e_gamma_neg_pt25_Q2gt30[number_dphi_e_ph_bins] = {0, 0, 0, 0, 0, 0};

	// Q2gt30
		//pT_cut in center-of-mass frame = 2.5 GeV/c
		Double_t * all_theory_cs_font_pt25_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pt25_Q2gt30, font_xp_pt25_Q2gt30, 0, font_deta_pt25_Q2gt30, font_dphi_e_gamma_pt25_Q2gt30, 0};
		Double_t * all_theory_cs_font_pt25_pos_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt25_Q2gt30, font_xp_pos_pt25_Q2gt30, 0, font_deta_pos_pt25_Q2gt30, font_dphi_e_gamma_pos_pt25_Q2gt30, 0};
		Double_t * all_theory_cs_font_pt25_neg_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt25_Q2gt30, font_xp_neg_pt25_Q2gt30, 0, font_deta_neg_pt25_Q2gt30, font_dphi_e_gamma_neg_pt25_Q2gt30, 0};
		//pT_cut in center-of-mass frame = 0.5 GeV/c
		Double_t * all_theory_cs_font_pt05_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pt05_Q2gt30, font_xp_pt05_Q2gt30, 0, font_deta_pt05_Q2gt30, font_dphi_e_gamma_pt05_Q2gt30, 0};
		Double_t * all_theory_cs_font_pt05_pos_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt05_Q2gt30, font_xp_pos_pt05_Q2gt30, 0, font_deta_pos_pt05_Q2gt30, font_dphi_e_gamma_pos_pt05_Q2gt30, 0};
		Double_t * all_theory_cs_font_pt05_neg_Q2gt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt05_Q2gt30, font_xp_neg_pt05_Q2gt30, 0, font_deta_neg_pt05_Q2gt30, font_dphi_e_gamma_neg_pt05_Q2gt30, 0};
		
	// Q2lt30
		//pT_cut in center-of-mass frame = 2.5 GeV/c
		Double_t * all_theory_cs_font_pt25_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pt25_Q2lt30, font_xp_pt25_Q2lt30, font_dphi_pt25_Q2lt30, font_deta_pt25_Q2lt30, font_dphi_e_gamma_pt25_Q2lt30, font_deta_e_gamma_pt25_Q2lt30};
		Double_t * all_theory_cs_font_pt25_pos_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt25_Q2lt30, font_xp_pos_pt25_Q2lt30, font_dphi_pos_pt25_Q2lt30, font_deta_pos_pt25_Q2lt30, font_dphi_e_gamma_pos_pt25_Q2lt30, font_deta_pos_e_gamma_pt25_Q2lt30};
		Double_t * all_theory_cs_font_pt25_neg_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt25_Q2lt30, font_xp_neg_pt25_Q2lt30, font_dphi_neg_pt25_Q2lt30, font_deta_neg_pt25_Q2lt30, font_dphi_e_gamma_neg_pt25_Q2lt30, font_deta_neg_e_gamma_pt25_Q2lt30};
		//pT_cut in center-of-mass frame = 0.5 GeV/c
		Double_t * all_theory_cs_font_pt05_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pt05_Q2lt30, font_xp_pt05_Q2lt30, 0, font_deta_pt05_Q2lt30, font_dphi_e_gamma_pt05_Q2lt30, 0};
		Double_t * all_theory_cs_font_pt05_pos_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt05_Q2lt30, font_xp_pos_pt05_Q2lt30, 0, font_deta_pos_pt05_Q2lt30, font_dphi_e_gamma_pos_pt05_Q2lt30, 0};
		Double_t * all_theory_cs_font_pt05_neg_Q2lt30[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt05_Q2lt30, font_xp_neg_pt05_Q2lt30, 0, font_deta_neg_pt05_Q2lt30, font_dphi_e_gamma_neg_pt05_Q2lt30, 0};
		
	// Full region
		//pT_cut in center-of-mass frame = 2.5 GeV/c
		Double_t * all_theory_cs_font_pt25[12] 	   = {0, 0, 0, 0, 0, 0, font_xgamma_pt25_Q2full, font_xp_pt25_Q2full, font_dphi_pt25_Q2full, font_deta_pt25_Q2full, font_dphi_e_gamma_pt25_Q2full, font_deta_e_gamma_pt25_Q2full};
		Double_t * all_theory_cs_font_pt25_pos[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt25_Q2full, font_xp_pos_pt25_Q2full, font_dphi_pos_pt25_Q2full, font_deta_pos_pt25_Q2full, font_dphi_e_gamma_pos_pt25_Q2full, font_deta_e_gamma_pos_pt25_Q2full};
		Double_t * all_theory_cs_font_pt25_neg[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt25_Q2full, font_xp_neg_pt25_Q2full, font_dphi_neg_pt25_Q2full, font_deta_neg_pt25_Q2full, font_dphi_e_gamma_neg_pt25_Q2full, font_deta_e_gamma_neg_pt25_Q2full};
		//pT_cut in center-of-mass frame = 0.5 GeV/c
		Double_t * all_theory_cs_font_pt05[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pt05_Q2full, font_xp_pt05_Q2full, 0, font_deta_pt05_Q2full, font_dphi_e_gamma_pt05_Q2full, 0};
		Double_t * all_theory_cs_font_pt05_pos[12] = {0, 0, 0, 0, 0, 0, font_xgamma_pos_pt05_Q2full, font_xp_pos_pt05_Q2full, 0, font_deta_pos_pt05_Q2full, font_dphi_e_gamma_pos_pt05_Q2full, 0};
		Double_t * all_theory_cs_font_pt05_neg[12] = {0, 0, 0, 0, 0, 0, font_xgamma_neg_pt05_Q2full, font_xp_neg_pt05_Q2full, 0, font_deta_neg_pt05_Q2full, font_dphi_e_gamma_neg_pt05_Q2full, 0};
#endif