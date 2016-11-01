#define selector_c
// The class definition in selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts ,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("selector.C")
// Root > T->Process("selector.C","some options")
// Root > T->Process("selector.C+")
//
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
//#include "selector.h"
//#ifndef constants_h
//#error include constants.h in selector.c
#include "constants.h"
//#endif
//#ifdef selector_cxx
#include "selector.h"
#include "runinfo.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TSelector.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include "Dijets.h"
#include "BCR_CalibTools.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"


using KtJet::KtLorentzVector;
using KtJet::KtEvent;

extern Bool_t excludedrunlist(Int_t);
extern Bool_t take_run(Int_t);

extern Bool_t poltake(Int_t, Int_t, Int_t, Int_t, RunInfo runinfo[], const int);
extern Double_t polcorr0405e(Double_t );
extern Double_t polcorr0607p(Double_t );
extern Double_t polcorr06e(Double_t );
//extern Double_t jet_en_corr(Double_t eta, Double_t et);

extern "C"
{
	float pcele_(int &Imod,float &Emes,float &Emip,float &The,float &Phi);
}
/*turned off as usual
	//extern Double_t reweighting(Double_t , TString period, TString mc_type);
	//extern void fill_trigger_bits(TH1D* h_dis, Int_t r_Fltw[], Int_t n);
	//extern Bool_t flt_bit_is_1(UInt_t bit, Int_t r_Fltw[], Int_t n);
	//extern Double_t zvtx_reweighting_0607p(const Double_t x);

	//Int_t findMaxEt(Int_t sorted[], Double_t et_jets[], Int_t N);
	//void sortEt(TObjArray *arrayNotSorted, TObjArray *arraySorted);
*/
Bool_t selector::Process()
{	
	/*
		// The Process() function is called for each entry in the tree (or possibly
		// keyed object in the case of PROOF) to be processed. The entry argument
		// specifies which entry in the currently loaded tree is to be processed.
		// It can be passed to either selector::GetEntry() or TBranch::GetEntry()
		// to read either all or the required parts of the data. When processing
		// keyed objects with PROOF, the object is already loaded and is available
		// via the fObject pointer.
		//
		// This function should contain the "body" of the analysis. It can contain
		// simple or elaborate selection criteria, run algorithms on the data
		// of the event and typically fill histograms.
		//
		// The processing can be stopped by calling Abort().
		//
		// Use fStatus to set the return value of TTree::Process().
		//
		// The return value is currently not used.
		//  cout <<"geting tree... " ;//<< endl;
	*/
	
	Long64_t nentries = fChain->GetEntries();
	if (!nodebugmode) cout <<"Number of events to process: " << nentries << endl;

	//parameter initialising
		Bool_t take_had_event = kFALSE;
		Bool_t take_part_event = kFALSE;

		Bool_t breakindex = kFALSE;
		Bool_t take_event = kTRUE;
		Bool_t take_event_trig = kTRUE;
		here_is_prph = kFALSE;
		Bool_t here_is_jet = kFALSE;
		Bool_t take_jet = kTRUE;
		Bool_t print_interesting_events = kFALSE;
		Bool_t study_jets = kFALSE;
		phot_count = 0;

	//debugging tools
		ofstream list_of_runs("runs" + period);
		Int_t Runnr_prev = 0;
		Bool_t debugcontinue = kTRUE;
		Int_t missed = 0;
		int num_had(0), num_part(0);
	check_cuts = kTRUE;
	for(Long64_t entry = 0; entry < nentries -1  && debugcontinue; entry++)
	{
		en_mom_conservation = true;
		// To test specific entry
			// if (entry < 4 ) continue;
			// else 
			// if (entry > 4) 
			// {
			// 	cout << "doh 0 " << endl;
			// 	exit(0);
			// 	cout << "doh" << endl;
			// }
		//
		
			// It goes by the end by:
			// 			entry:2  Eventnr:13
			// In the entry:0 :
			// 			take_event = 1, here_is_jet = 0, here_is_prph = 0
			// In the entry:1 :
			// 			eta = -0.605987, wtx_eta = 1.33793, et = 3.42134, wtx_et = 1
			// 			q2 = 197.658, wtx_q2 = 0.431328
			// 			electron cand 0: E = 6.92512 2.13603 177.263
			// 			rejected by Siecorr[0][2] = 6.92512
			// 			rejected by Sith = 122.385
			// 			here_is_prph = 0, take_event = 0
			// 			send to prpktjetb -1
						  
			// 			 CALL TO CELE SUBROUTINE with option 20
			// 			 ======================================
			// 			 BPRE corrects MC  without non-uniform
			// 			 ======================================
			// 			here_is_prph = 0, take_event = 0

		
		fChain->GetEntry(entry);
		//if (Eventnr != 26) continue;
		//if (Eventnr > 26 || entry > 26) exit(1);
		cout << "entry: " << entry << "; Eventnr: " << Eventnr << "; Runnr_prev: " << Runnr_prev << endl;

		//save list of runnumbers
			if (Runnr != Runnr_prev) 
			{
				list_of_runs << Runnr << endl;
				Runnr_prev = Runnr;
			}
			//hist.hist_runnr->Fill(Runnr);
		
		wtx = 1.;

		//Missed compare to full code - but I added
			// if (!Data) 
			// {
			// 	Double_t et_true_photon;
			// 	Int_t index_true_photon = -1;
			// 	TLorentzVector true_photon;

			// 	if (mc_type == "mc_prph")    
			// 	{
			// 		for (Int_t i = 0; i < Npart; i++)       
			// 		{
			// 			if (Part_id[i]==12)     
			// 			{
			// 				index_true_photon = i;
			// 				break;
			// 			}
			// 		}
			// 		if (index_true_photon >= 0) true_photon.SetPxPyPzE(Part_p[index_true_photon][0], Part_p[index_true_photon][1], Part_p[index_true_photon][2], Part_p[index_true_photon][3]);
			// 		else cerr << "photon not found" << endl; exit(-1);

			// 		wtx *= etaph_reweighting(true_photon.Eta(), mc_type) ;
			// 		//  wtx *= etph_reweighting(true_photon.Et(), mc_type);
			// 		wtx *= eph_reweighting(true_photon.E(), mc_type) ;
			// 	}
			// }
		
			// commented for test
			if (!Data && mc_type != "mc_bg_rad")
			{
				//wtx *= q2_reweighting(Mc_q2, mc_type); //warning
			}
		if (wtx != 1.0) cout << "wtx: " << wtx << endl;
		{     
			//    if(Eventnr == 21334) check_cuts = kTRUE;//52258 22172
			
			//parameter initialising
				take_had_event = kFALSE;
				take_part_event = kFALSE;

				take_event = kTRUE;
				take_event_trig = kTRUE;
				here_is_prph = kFALSE;
				here_is_jet = kFALSE;
				hist.zvtx->Fill(Zvtx, wtx);
			//detector lev selection - was commented as not needed
				// 	det_cross_sec_bin_et = -1;
				// 	det_cross_sec_bin_eta = -1;
				// 	det_cross_sec_bin_q2 = -1;
				// 	det_cross_sec_bin_x = -1;
				// 	det_cross_sec_bin_et_jet = -1;
				// 	det_cross_sec_bin_eta_jet = -1;

				// 	det_cross_sec_bin_xgamma = -1;
				// 	det_cross_sec_bin_xp = -1;
				// 	det_cross_sec_bin_dphi = -1;
				// 	det_cross_sec_bin_deta = -1;
				// 	det_cross_sec_bin_dphi_e_ph = -1;
				// 	det_cross_sec_bin_deta_e_ph = -1;
				// 	deltaRtrack = 999.;   

				// if (Runnr >= 47010 && Runnr <= 51245) period = "04p";
				// if (Runnr >= 52258 && Runnr <= 57123) period = "0405e";
				// if (Runnr >= 58207 && Runnr <= 59947) period = "06e";
				// if (Runnr >= 60005 && Runnr <= 62639) period = "0607p";

				// //electron selection
				// 	Int_t sinistra_electron_number = -1;//
				// 	Bool_t take_electron_event = kTRUE;

				// 	if (Sincand < 1)
				// 	{
				// 		take_event = kFALSE;
				// 		if (check_cuts) cout << "rejected by Sincand = " << Sincand << endl;
				// 	}
				// 	else if (check_cuts)
				// 		for(Int_t i = 0; i < Sincand; i++)
				// 			cout << "electron cand " << i << ": E = " << Siecorr[i][2] << " " << Sith[i] << " " << Siq2el[i] << endl;
					
				// 	//find electron candidate passing cuts - sinistra_electron_number
				// 		for(Int_t i = 0; i < 1; i++) 
				// 		{
				// 			take_electron_event = kTRUE;
				// 			if (Siq2el[i] < 10. || Siq2el[i] > 350.) take_electron_event = kFALSE;
				// 			if (Siecorr[i][2] < 10.)
				// 			{
				// 				take_electron_event = kFALSE;
				// 				if (check_cuts) cout << "rejected by Siecorr[0][2] = " << Siecorr[i][2] << endl;
				// 			}  
				// 			if (Sith[i] * 180.0 / TMath::Pi() < 139.8 || Sith[i] * 180.0 / TMath::Pi() > 180.0) 
				// 			{
				// 				take_electron_event = kFALSE;
				// 				if (check_cuts) cout << "rejected by Sith = " << Sith[0] * 180./TMath::Pi() << endl;
				// 			}
				// 			if ((TMath::Abs(Sipos[i][0]) < 14.8) &&
				// 					(Sipos[i][1]  >-14.6)  &&
				// 					(Sipos[i][1]  < 12.5))  
				// 			{
				// 				take_electron_event = kFALSE;
				// 				if (check_cuts) cout << "rejected by Sipos[0][0] = " << Sipos[0][0] << ", Sipos[0][1] = " << Sipos[0][1] << endl;
				// 			}
							
				// 			if (take_electron_event)
				// 			{
				// 				sinistra_electron_number = i;
				// 				break;
				// 			}
				// 		}

				// 		if (!take_electron_event) 
				// 		{
				// 			cout << "rejected by e investigation, take_electron_event: " << take_electron_event << endl;
				// 			take_event = kFALSE;
				// 		}
				// 		else 
				// 		{
				// 			if  (check_cuts)  cout << " =============>passed e selection" << endl;
				// 			//cout << entry << endl;
				// 			//exit(-1);
				// 		}

						

				// //kinematical cuts
				// 	if (Zvtx < -40. || Zvtx > 40.) 
				// 	{
				// 		take_event = kFALSE;
				// 		if (check_cuts) cout << "rejected by Zvtx = " << Zvtx << endl;
				// 	}

				// 	if (Cal_empz < 35. || Cal_empz > 65.) 
				// 	{
				// 		take_event = kFALSE;
				// 		if (check_cuts) cout << "rejected by Cal_empz = " << Cal_empz << endl;
				// 	}

				// 	Int_t nVertTracks = 0;
				// 	for(Int_t ii = 0;ii < Trk_ntracks;ii++)
				// 	{
				// 		if (Trk_prim_vtx[ii] != 1) continue;
				// 		TVector3 tr(Trk_px[ii],Trk_py[ii],Trk_pz[ii]);
				// 		if (tr.Mag() <= 0.25) continue;
				// 		if (tr.Theta() >= 2.44) continue;
				// 		nVertTracks++;
				// 	}

				// 	if (nVertTracks < 1) 
				// 	{
				// 		take_event = kFALSE;
				// 		if (check_cuts) cout << "rejected by nVertTracks = " << nVertTracks << endl;
				// 	}

				// 	Bool_t flt_trigger = kTRUE;
				// 	if ((Fltw[0] & (1 << 27)) == 0 &&
				// 			(Fltw[0] & (1 << 29)) == 0 &&
				// 			(Fltw[1] & (1 << 3)) == 0 &&
				// 			(Fltw[1] & (1 << 6)) == 0 &&
				// 			(Fltw[1] & (1 << 7)) == 0 &&
				// 			(Fltw[1] & (1 << 8)) == 0 &&
				// 			(Fltw[1] & (1 << 10)) == 0 &&
				// 			(Fltw[1] & (1 << 11)) == 0 &&
				// 			(Fltw[1] & (1 << 13)) == 0 &&
				// 			(Fltw[1] & (1 << 14)) == 0 )
				// 		flt_trigger = kFALSE;

				// 	Bool_t slt_trigger = kTRUE;
				// 	if ((Sltw[5] & (1 << 12)) == 0 &&
				// 			(Sltw[5] & (1 << 13)) == 0 &&
				// 			(Sltw[5] & (1 << 14)) == 0 &&
				// 			(Sltw[5] & (1 << 0)) == 0 &&
				// 			(Sltw[5] & (1 << 6)) == 0 )
				// 		slt_trigger = kFALSE;

				// 	Bool_t tlt_trigger = kTRUE;
				// 	if (period == "0405e" || period == "04p") 
				// 		if ((Tltw[3] & (1 << 2)) == 0 &&
				// 				(Tltw[2] & (1 << 1)) == 0 &&
				// 				(Tltw[2] & (1 << 8)) == 0 &&
				// 				(Tltw[13] & (1 << 0)) == 0)
				// 			tlt_trigger = kFALSE;
				// 	else if (period == "0607p" || period == "06e")
				// 		if ((Tltw[3] & (1 << 2)) == 0 &&
				// 				(Tltw[2] & (1 << 1)) == 0 &&
				// 				(Tltw[2] & (1 << 8)) == 0 &&
				// 				(Tltw[13] & (1 << 0)) == 0)
				// 			tlt_trigger = kFALSE;

				// 	Bool_t spp01trigger = kTRUE;
				// 	if ((Tltw[2] & (1 << 0)) == 0) spp01trigger = kFALSE;

				// 	Bool_t spp02trigger = kTRUE;
				// 	if ((Tltw[2] & (1 << 1)) == 0) spp02trigger = kFALSE;

				// 	Bool_t spp09trigger = kTRUE;
				// 	if ((Tltw[2] & (1 << 8)) == 0) spp09trigger = kFALSE;

				// 	Bool_t spp_trigger;
				// 	if (period == "0405e" || period == "04p") spp_trigger = spp02trigger;  ///WARNING !!! MUST be 02!
				// 	else if (period == "06e" || period == "0607p") spp_trigger = spp09trigger; ///WARNING !!! MUST be 09!

				// 	//      if (mc_type == "mc_prph")
				// 	//    spp_trigger = spp02trigger;  ///WARNING !!! MUST be 02!

				// 	//
				// 	//trigger efficiency study
				// 	//

				// 	//      if (!(flt_trigger && slt_trigger && tlt_trigger)) { //warning
				// 	if (!spp_trigger) 
				// 	{
				// 		take_event_trig = kFALSE;
				// 		if (check_cuts) cout << "rejected by trigger condition " << endl;
				// 	}

				// 	//      cout << "dis sel: " << Runnr << " " << Eventnr << " " << Siq2el[0] << " " << Zvtx << " " << Cal_empz << " " << Sith[0] << " " 
				// 	//       << Sipos[0][0] << " " << Sipos[0][1] << " " << spp02trigger << endl;
				
				// //parameter initialising
				// 	TLorentzVector v_electron; 
				// 	Int_t electron_number = -1;
				// 	max_et_candidate_number = -1;
				// 	Double_t dr_elec_min = 999.;
				// 	Double_t de_elec_min = 999.;
				// 	Double_t n_close_tracks = 0.;
				// 	max_et_candidate_number = -1;
				// 	candidate_jet_number = -1;
				// 	radiated_candidate_number = -1;
				// 	max_et_candidate = -999.;
				// 	max_et_candidate_uncorr = -999.;
				// 	jet_energy_frac = 0.;
				// 	cell_energy_frac = 0.;
				// 	Double_t pt_of_event = 0.;
				// 	Double_t px_of_event = 0.;
				// 	Double_t py_of_event = 0.;
				// 	Double_t pz_of_event = 0.;
				// 	Double_t e_of_event = 0.;
				// 	Double_t et_of_event = 0.;
				// 	Double_t e_cal_event = 0.;
				// 	Double_t e_had_event = 0.;
				// 	Double_t e_emc_event = 0.;
				// 	Double_t pt_sqrt_et;

				// 	if (!take_event_trig) take_event = kFALSE;

				// 	if (take_event) 
				// 	{
				// 		cout << "========>passed kinematics:" << endl;
				// 		//cout << entry << endl;
				// 		//exit(-1);
				// 	}
				// 	else 
				// 	{
				// 		cout << "		didn't passed kinematics:" << endl;
				// 		continue;
				// 		//exit(-1);
				// 	}

				// //check whether electron excluded from zufos , yjb
				// 	if (take_event)
				// 	{
				// 		//finding electron among ZUFO
				// 			v_electron.SetPtEtaPhiE( Sipt[sinistra_electron_number], 
				// 				-TMath::Log( TMath::Tan( Sith[sinistra_electron_number] / 2.)), 
				// 				Siph[sinistra_electron_number], 
				// 				Siecorr[sinistra_electron_number][2]);

				// 			for(Int_t zloop = 0; zloop < Nzufos; zloop++)
				// 			{
				// 				hist.dis_zufo_type->Fill(Tufo[zloop][0], wtx);
				// 				e_cal_event += Zufoecal[zloop];
				// 				e_had_event += Zufoecal[zloop] - Zufoeemc[zloop];
				// 				e_emc_event += Zufoeemc[zloop];
				// 				pt_of_event += TMath::Sqrt(Zufo[zloop][0]*Zufo[zloop][0] + Zufo[zloop][1]*Zufo[zloop][1]);
				// 				px_of_event += Zufo[zloop][0];
				// 				py_of_event += Zufo[zloop][1];
				// 				pz_of_event += Zufo[zloop][2];
				// 				e_of_event += Zufo[zloop][3];
				// 				myKzufodeltaz[zloop] = Kzufodeltaz[zloop] / 5.45;

				// 				TLorentzVector v(Zufo[zloop][0], Zufo[zloop][1], Zufo[zloop][2], Zufo[zloop][3]);
				// 				Double_t et = v.Et();
				// 				et_of_event += et;
				// 				Double_t dr = v.DeltaR(v_electron);

				// 				if (v.E() > 5. && dr < dr_elec_min)
				// 				{
				// 					dr_elec_min = dr;
				// 					de_elec_min = v.E() - v_electron.E();
				// 					if (dr < 0.5) electron_number = zloop;
				// 				}
				// 				if (dr < 0.2) n_close_tracks++;
				// 			}
				// 			if (check_cuts) cout << "electron found: zufo # " << electron_number << endl;
				// 			if (electron_number < 0) 
				// 			{
				// 				take_event = kFALSE;
				// 				if (check_cuts) cout << "event rejected because was not found zufo that is close to sinistra electron" << endl;
				// 			}

				// 		//yjb for x_gamma variable
				// 			Empz = 0.;
				// 			for(Int_t zloop=0; zloop<Nzufos; zloop++)
				// 			{
				// 				if (zloop != electron_number)
				// 					Empz += Zufo[zloop][3] - Zufo[zloop][2];
				// 			}     
				// 			yjb = Empz / (2. * E_e);

				// 		// filling inclusive DIS histograms
							
				// 			hist.dis_siyjb->Fill(Siyjb[0], wtx);
				// 			hist.dis_siyjb_cell->Fill(Siyjb_cell[0], wtx);
				// 			hist.dis_handmade_yjb->Fill(yjb, wtx);

				// 			hist.dis_Q2_x->Fill(Siq2el[0], Sixel[0], wtx);
				// 			hist.dis_Q2->Fill(Siq2el[0], wtx);
				// 			hist.dis_Q2_el->Fill(Siq2el[0], wtx);
				// 			hist.dis_Q2_da->Fill(Siq2da[0], wtx);
				// 			hist.dis_Q2_da_cell->Fill(Siq2da_cell[0], wtx);
				// 			hist.dis_Q2_jb->Fill(Siq2jb[0], wtx);
				// 			hist.dis_Q2_jb_cell->Fill(Siq2jb_cell[0], wtx);
				// 			hist.dis_px_event->Fill(px_of_event, wtx);
				// 			hist.dis_py_event->Fill(py_of_event, wtx);
				// 			hist.dis_pz_event->Fill(pz_of_event, wtx);
				// 			hist.dis_pt_event->Fill(pt_of_event, wtx);
				// 			hist.dis_e_event->Fill(e_of_event, wtx);
				// 			hist.dis_et_event->Fill(et_of_event, wtx);
				// 			hist.dis_n_zufos->Fill(Nzufos, wtx);
				// 			hist.dis_emc->Fill(e_emc_event, wtx);
				// 			hist.dis_had->Fill(e_had_event, wtx);
				// 			hist.dis_cal->Fill(e_cal_event, wtx);
				// 			hist.dis_had_tot_ratio->Fill(e_had_event/e_cal_event, wtx);
							
				// 			hist.dis_elec_zufo_dr_min->Fill(dr_elec_min, wtx);
				// 			hist.dis_elec_zufo_de_min->Fill(de_elec_min, wtx);
				// 			hist.dis_zufos_near_electron->Fill(n_close_tracks, wtx);
				// 			hist.dis_zvtx->Fill(Zvtx, wtx);
				// 			hist.dis_n_vtx_tracks->Fill(nVertTracks, wtx);
				// 			hist.dis_cal_empz->Fill(Cal_empz, wtx);
				// 			hist.dis_zufo_empz->Fill(Empz, wtx);
				// 			pt_sqrt_et = pt_of_event/TMath::Sqrt(et_of_event);
				// 			hist.dis_pt_sqrtet->Fill(pt_sqrt_et, wtx);
				// 			hist.dis_x_el->Fill(Sixel[0], wtx);
				// 			hist.dis_x_el_log->Fill(Sixel[0], wtx);
				// 			hist.dis_y_el->Fill(Siyel[0], wtx);
				// 			hist.dis_electron_e->Fill(Siecorr[0][2], wtx);
				// 			hist.dis_electron_theta->Fill(Sith[0]*180.0/TMath::Pi(), wtx);
				// 			hist.dis_electron_probability->Fill(Siprob[0], wtx);

				// 		//printout interesting events
				// 			if(print_interesting_events) 
				// 			{
				// 				if(pt_sqrt_et > 10.) {
				// 					if (nodebugmode) cout << endl;
				// 					if (nodebugmode) cout << "============ INTERESTING EVENT ================" << endl;
				// 					if (nodebugmode) cout << "========= " << Runnr << " " << Eventnr << " ================" << endl;
				// 					if (nodebugmode) cout << " seems to cosmic: pt_sqrt_et = " << pt_sqrt_et << endl;
				// 					if (nodebugmode) cout << "===============================================" << endl << endl;      
				// 				}
				// 				if(pz_of_event > 500.) {
				// 					if (nodebugmode) cout << endl;
				// 					if (nodebugmode) cout << "============ INTERESTING EVENT ================" << endl;
				// 					if (nodebugmode) cout << "========= " << Runnr << " " << Eventnr << " ================" << endl;
				// 					if (nodebugmode) cout << " proton fucked up: pz_of_event = " << pz_of_event << endl;
				// 					if (nodebugmode) cout << "===============================================" << endl << endl;      
				// 				}
				// 			}
				// 	}
				
				// 	if (check_cuts) cout << "here_is_prph = " << here_is_prph << ", take_event = " << take_event << endl;
				
				// //finding photon
				// 	//finding some photon index by ID in MC --> index_true_rad_photon
				// 		Int_t index_true_rad_photon = -1;
				// 		if (!Data) 
				// 		{	
				// 			cout << "Npart = " << Npart << endl;
				// 			for(Int_t i = 0; i < Npart; i++)
				// 			{
				// 				cout << "Fmckin instance: " << i << " :: Part_id =" << Part_id[i]  << ", Part_prt = "<< Part_prt[i] << endl;
				// 				if (Part_prt[i] == 29)
				// 				{
				// 					index_true_rad_photon = i;
				// 					cout << "true rad photon found: " << index_true_rad_photon << endl;
				// 					break;
				// 				}
				// 			}
				// 			// OLd not working and commented
				// 				// cout << "Photn = " << Photn << endl;
				// 				// for(Int_t ploop = 0; ploop < Photn; ploop++)
				// 				// {
				// 				// 	cout << "photid[" << ploop << "] = " << Photid[ploop] << endl;
				// 				// 	if (Photid[ploop] == 5)
				// 				// 	{
				// 				// 		index_true_rad_photon = ploop;
				// 				// 		cout << "true rad photon found: " << index_true_rad_photon << endl;
				// 				// 		break;
				// 				// 	}
				// 				// }
				// 		}

				// 	//finding photon in Data and MC -->  max_et_candidate_number, max_et_candidate, v_uncorr_prompt_photon, v_corr_prompt_photon, glob_deltaz
				// 	if (use_ktjetb_prph)//kTrue
				// 	{
				// 		if (check_cuts && !nodebugmode) cout << "send to prpktjetb " << index_true_rad_photon << endl;
				// 		if (SelectPrPhKtJetB(index_true_rad_photon, electron_number)) here_is_prph = kTRUE;
				// 		//Missed compared to prev code - not presented in prev code
				// 		if (!here_is_prph) take_event = kFALSE;
				// 		if (!here_is_prph) continue;
				// 	}
				// if(check_cuts && !nodebugmode) cout << "here_is_prph = " << here_is_prph << ", take_event = " << take_event << endl;
				
				
				// //
				// if (take_event)
				// {
				// 	cout << "max_et_candidate_number: "<< max_et_candidate_number << endl;
				// 	cout << "constructing vector on/of zufos without e and photon" <<endl;
				// 	//constructing vector on/of zufos without e and photon - input_vec, input_vec_to_zufos
				// 		vector<KtLorentzVector> input_vec;
				// 		vector<Int_t> input_vec_to_zufos;
				// 		for(Int_t zloop = 0; zloop < Nzufos; zloop++)
				// 		{
				// 			if (electron_number > -1 && zloop == electron_number)
				// 			{
				// 				if (false && check_cuts) cout << "rejecting zufo number " << electron_number << " with e = " << Zufo[zloop][3] << " because it is most probably electron " << endl;
				// 				continue;
				// 			}
				// 			if (use_ktjetb_prph && here_is_prph && zloop == max_et_candidate_number) 
				// 			{
				// 				// Missed compared to full code
				// 				//but I added
				// 					index_phot_vector = input_vec_to_zufos.size();
				// 					//    continue;
								
				// 				if (check_cuts) 
				// 					cout << "rejecting zufo number " << max_et_candidate_number 
				// 						<< " with e = " << Zufo[zloop][3] 
				// 						<< " because it is prph candidate " << endl;
				// 				continue;
				// 			}
				// 			KtLorentzVector p(Zufo[zloop][0], Zufo[zloop][1], Zufo[zloop][2], Zufo[zloop][3]);
				// 			input_vec.push_back(p);
				// 			input_vec_to_zufos.push_back(zloop);// check if this is 1 2 3 4 5 ...

				// 		}
				// 		if (check_cuts) cout << Nzufos << " " << input_vec.size() << " hadrons added" << endl;
					
				// 	//constructing jets
				// 	int type  = 3; // pe
				// 	int angle = 2; // deltaR
				// 	int recom = 1; // !pt
				// 	double rparameter = 1.0;

				// 	// Construct the KtEvent object 
				// 	KtEvent ev(input_vec, type, angle, recom, rparameter);//Missed compare to full code

				// 	// Retrieve the final state jets from KtEvent sorted by Et
				// 	// before it was in form of:  vector<KtLorentzVector> jets;
				// 	vector<KtLorentzVector> jets = ev.getJetsEt();

				// 	//Since it is not used - Commented out
						
				// 			// // find jet containing prompt photon from HepForge 
				// 			// if (!use_ktjetb_prph) //()=kFalse -> candidate_jet_number, glob_deltaz, v_uncorr_prompt_photon, v_corr_prompt_photon, v_prompt_photon_jet, max_et_candidate, max_et_candidate_number,
				// 			// {
								
				// 			// 	// find closest by distance to radiated photon zufo - radiated_candidate_number
				// 			// 		if (!Data && index_true_rad_photon > -1) 
				// 			// 		{
				// 			// 			TVector3 v_true_rad_photon(Photp[index_true_rad_photon][0], Photp[index_true_rad_photon][1], Photp[index_true_rad_photon][2]);
				// 			// 			Double_t dr, dr_min = 999.;
				// 			// 			for(Int_t zloop = 0; zloop < input_vec.size(); zloop++)
				// 			// 			{
				// 			// 				Double_t zufo_px = input_vec[zloop].px(), zufo_py = input_vec[zloop].py(), zufo_pz = input_vec[zloop].pz();
				// 			// 				if (TMath::Sqrt((zufo_px*zufo_px + zufo_py*zufo_py) / (zufo_pz*zufo_pz)) < 1.e-10 
				// 			// 					|| zufo_px!=zufo_px || zufo_py!=zufo_py || zufo_pz!=zufo_pz 
				// 			// 					|| zufo_px*zufo_px+zufo_py*zufo_py+zufo_pz*zufo_pz<1.e-200)
				// 			// 					continue;
				// 			// 				TVector3 v_zufo(input_vec[zloop].px(), input_vec[zloop].py(), input_vec[zloop].pz());
				// 			// 				dr = v_zufo.DeltaR(v_true_rad_photon);
				// 			// 				if (dr < dr_min && dr < 0.2) 
				// 			// 				{
				// 			// 					dr_min = dr;
				// 			// 					radiated_candidate_number = zloop;
				// 			// 				}
				// 			// 			}
				// 			// 		}

				// 			// 	//find prph candidate - its energy 
				// 			// 		for(Int_t jloop = 0; jloop < jets.size(); jloop++)
				// 			// 		{
				// 			// 			for(Int_t zloop = 0; zloop < input_vec.size(); zloop++)
				// 			// 			{
				// 			// 				// if this is LL MC then only consider radiated photon
				// 			// 					if (!Data && mc_type=="mc_bg_rad" 
				// 			// 						&& radiated_candidate_number != zloop)
				// 			// 						continue;

				// 			// 				// if this is bg MC then do NOT consider radiated photon
				// 			// 					if (!Data && mc_type=="mc_bg_norad" 
				// 			// 						&& radiated_candidate_number == zloop)
				// 			// 						continue;

				// 			// 				// if data or QQ MC
				// 			// 				if ( jets[jloop].contains(input_vec[zloop]) ) 
				// 			// 				{           
				// 			// 					Int_t ntrack = 0;
				// 			// 					Bool_t is_prph_candidate = kTRUE;

				// 			// 					//calculating number of close tracks for zufo
				// 			// 						for (Int_t i = 0; i < Trk_ntracks; i++)
				// 			// 						{
				// 			// 							TVector3 tr;          
				// 			// 							tr.SetXYZ(Trk_px[i],Trk_py[i],Trk_pz[i]);
				// 			// 							if (tr.Mag() < 0.25) continue;             
				// 			// 							// Calculating phi, eta of track
				// 			// 							Double_t ddr = (tr.PseudoRapidity() - input_vec[zloop].eta())
				// 			// 										  *(tr.PseudoRapidity() - input_vec[zloop].eta());
				// 			// 							Double_t dphi = tr.Phi()-input_vec[zloop].phi();
				// 			// 							if (dphi > TMath::Pi())  dphi -= 2.*TMath::Pi();
				// 			// 							if (dphi < -TMath::Pi()) dphi += 2.*TMath::Pi();
				// 			// 							ddr = TMath::Sqrt(ddr + dphi*dphi);
				// 			// 							if (ddr < 0.2) ntrack++; 
				// 			// 						}
				// 			// 						//        if (ntrack > 0) take_event = kFALSE;
				// 			// 						//        cout << "jet: " << jets[jloop].px() << " " << jets[jloop].py() << " " << jets[jloop].pz() << " " << jets[jloop].et() << " " <<jets[jloop].eta() << endl;
				// 			// 						//        cout << "zufo: " << Zufo[zloop][0] << " " << Zufo[zloop][1] << " " << Zufo[zloop][2] << " " << input_vec[zloop].et() << " " <<input_vec[zloop].eta() << " " << input_vec[zloop].phi() << endl;
				// 			// 						//        cout << "matched " << jloop <<"th jet and " << zloop << "th zufo" << endl;
												
				// 			// 					//calculate distance between jet and zufo
				// 			// 						Double_t   dR = (jets[jloop].eta() - input_vec[zloop].eta()) 
				// 			// 						              * (jets[jloop].eta() - input_vec[zloop].eta());
				// 			// 						Double_t dphi = jets[jloop].phi() - input_vec[zloop].phi();
				// 			// 						if (dphi > TMath::Pi())  dphi = dphi - 2.*TMath::Pi();
				// 			// 						if (dphi < -TMath::Pi()) dphi = dphi + 2.*TMath::Pi();
				// 			// 						dR = TMath::Sqrt(dR + dphi*dphi);
				// 			// 						hist.phcand_jet_deltaR->Fill(dR, wtx);
				// 			// 						hist.phcand_jet_deltaEta->Fill(TMath::Abs(jets[jloop].eta() - input_vec[zloop].eta()), wtx);
				// 			// 						hist.phcand_jet_deltaPhi->Fill(TMath::Abs(jets[jloop].phi() - input_vec[zloop].phi()), wtx);
				// 			// 						/ *cout << Runnr << " " << Eventnr << " " << Tufo[zloop][0] << " "<< Zufoeemc[zloop]/Zufoecal[zloop] 
				// 			// 						  << " "<< input_vec[zloop].et() << " "
				// 			// 						  << jets[jloop].et() << " " << input_vec[zloop].et()/jets[jloop].et() << " " << dR << " " 
				// 			// 						  << input_vec[zloop].e() << " " << jets[jloop].e() << " " << input_vec[zloop].e()/jets[jloop].e() << endl;
				// 			// 						  * /

				// 			// 					// Photon(zufo) energy correction
				// 			// 						Float_t photon_eta   = input_vec[zloop].eta();
				// 			// 						Float_t photon_theta = input_vec[zloop].theta();
				// 			// 						Float_t photon_phi   = input_vec[zloop].phi();
				// 			// 						Float_t photon_e     = input_vec[zloop].e();
				// 			// 						TVector3 v_photon; v_photon.SetXYZ(input_vec[zloop].px(), input_vec[zloop].py(), input_vec[zloop].pz());
				// 			// 						Float_t bpres_mips = (Float_t)BPRES_mips(v_photon, bpres_conerad);
				// 			// 						//            cout << "calling correction " << pcele_mode << " " << photon_e << " " << bpres_mips << " " << photon_theta << " " << photon_phi << endl;
				// 			// 						Double_t photon_corr_e = pcele_(pcele_mode, photon_e, bpres_mips, photon_theta, photon_phi);// energy of the photon
				// 			// 						photon_corr_et = photon_corr_e * TMath::Sin(photon_theta);
				// 			// 						//            cout << "corrected et, e: " << photon_corr_et << " " << photon_corr_e << endl;
				// 			// 						cout << "prph candidate " << zloop  << " of " << input_vec.size() << ", eta = " << input_vec[zloop].eta() 
				// 			// 							<< ", et_corr = " << photon_corr_et << endl;
				// 			// 						//          if (Kzufofmax[input_vec_to_zufos[zloop]] < 0.05)
				// 			// 						//            {
				// 			// 						//              is_prph_candidate = kFALSE;
				// 			// 						//              if (check_cuts)
				// 			// 						//                cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by Kzufofmax[" << input_vec_to_zufos[zloop] << "] = " << Kzufofmax[input_vec_to_zufos[zloop]] << endl;
				// 			// 						//            }

				// 			// 					// cuts on zufo(photon)
				// 			// 						if (input_vec[zloop].eta() < -0.7 || input_vec[zloop].eta() > 0.9)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by input_vec[" << zloop << "].eta() = " << input_vec[zloop].eta() << endl;
				// 			// 						}
				// 			// 						if (photon_corr_et < 4.)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by input_vec[" << zloop << "].et() = " << input_vec[zloop].et() << endl;
				// 			// 						}
				// 			// 						if (photon_corr_et > 15.)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by input_vec[" << zloop << "].et() = " << input_vec[zloop].et() << endl;
				// 			// 						}
				// 			// 						if (Zufoecal[input_vec_to_zufos[zloop]] != 0 && Zufoeemc[input_vec_to_zufos[zloop]]/Zufoecal[input_vec_to_zufos[zloop]] < 0.9)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by Zufoeemc[" << input_vec_to_zufos[zloop] << "]/Zufoecal[" << input_vec_to_zufos[zloop] << "] = " << Zufoeemc[input_vec_to_zufos[zloop]]/Zufoecal[input_vec_to_zufos[zloop]] << endl;
				// 			// 						}
				// 			// 						if (Tufo[input_vec_to_zufos[zloop]][0] != 31)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by Tufo[" << input_vec_to_zufos[zloop] << "][0] = " << Tufo[input_vec_to_zufos[zloop]][0] << endl;
				// 			// 						}
				// 			// 						if (ntrack > 0)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by ntrack = " << ntrack << endl;
				// 			// 						}
				// 			// 						if (Zufoecal[input_vec_to_zufos[zloop]] == 0.)
				// 			// 						{
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by Zufoecal[" << input_vec_to_zufos[zloop] << "] = " << Zufoecal[input_vec_to_zufos[zloop]] << endl;
				// 			// 						}
				// 			// 						/ *           
				// 			// 							if (input_vec[zloop].et()/jets[jloop].et() < 0.9)
				// 			// 							{
				// 			// 								is_prph_candidate = kFALSE;
				// 			// 								if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by input_vec[" << zloop << "].et()/jets[" << jloop << "].et() = "
				// 			// 								<< input_vec[zloop].et()/jets[jloop].et() << endl;
				// 			// 							}
				// 			// 						* /
				// 			// 						if (photon_corr_e / jets[jloop].e() < 0.9)
				// 			// 						{
				// 			// 							//            is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by input_vec[" << zloop << "].e()/jets[" << jloop << "].e() = " 
				// 			// 									<< input_vec[zloop].e()/jets[jloop].e() << endl;
				// 			// 						}
				// 			// 						//if (input_vec[zloop].e()/jets[jloop].e() < 0.9)
				// 			// 						if (input_vec[zloop].e()/Kpjets[Kzufoidjet[input_vec_to_zufos[zloop]]-1][3] < 0.9) {//tak u Natashi 
				// 			// 							is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() 
				// 			// 									<< " rejected by input_vec[" << zloop << "].e()/Kpjets[" << Kzufoidjet[input_vec_to_zufos[zloop]]-1 << "][3] = " 
				// 			// 									<< input_vec[zloop].e()/Kpjets[Kzufoidjet[input_vec_to_zufos[zloop]]-1][3] << endl;
				// 			// 						}
				// 			// 						if (dR > 0.2)
				// 			// 						{
				// 			// 							//is_prph_candidate = kFALSE;
				// 			// 							if (check_cuts) cout << "candidate with eta " << input_vec[zloop].eta() << " rejected by dR = " << dR << endl;
				// 			// 						}

				// 			// 						/ *if (ntrack==0 && (Tufo[input_vec_to_zufos[zloop]][0] == 31) && (Zufoeemc[input_vec_to_zufos[zloop]]/Zufoecal[input_vec_to_zufos[zloop]] > 0.9) && 
				// 			// 						  input_vec[zloop].et()>4. && input_vec[zloop].et()<15. && input_vec[zloop].et()/jets[jloop].et()>0.9 && 
				// 			// 						  input_vec[zloop].e()/jets[jloop].e()>0.9 && dR<0.2 && input_vec[zloop].eta()>-0.7 &&
				// 			// 						  input_vec[zloop].eta()<0.9 && Kzufofmax[input_vec_to_zufos[zloop]] > 0.05)* /
												
				// 			// 					// if passed the photons cuts
				// 			// 						if (is_prph_candidate)
				// 			// 						{
				// 			// 							here_is_prph = kTRUE;
				// 			// 							jet_energy_frac = input_vec[zloop].e() / jets[jloop].e();
				// 			// 							cell_energy_frac = Zufoeemc[input_vec_to_zufos[zloop]] / Zufoecal[input_vec_to_zufos[zloop]];
				// 			// 							//if this was the biggest prph candidate ever
				// 			// 								if (photon_corr_et > max_et_candidate)
				// 			// 								{
				// 			// 									max_et_candidate = photon_corr_et;
				// 			// 									max_et_candidate_uncorr = input_vec[zloop].et();
				// 			// 									max_et_candidate_number = zloop;
				// 			// 									candidate_jet_number = jloop;
				// 			// 									v_uncorr_prompt_photon->SetPxPyPzE(input_vec[zloop].px(), 
				// 			// 											input_vec[zloop].py(), 
				// 			// 											input_vec[zloop].pz(), 
				// 			// 											input_vec[zloop].e());
				// 			// 									Double_t corfac = photon_corr_et / input_vec[zloop].et();
				// 			// 									v_corr_prompt_photon->SetPxPyPzE(corfac * input_vec[zloop].px(),
				// 			// 											corfac * input_vec[zloop].py(),
				// 			// 											corfac * input_vec[zloop].pz(),
				// 			// 											corfac * input_vec[zloop].e());
				// 			// 									v_prompt_photon_jet->SetPxPyPzE(jets[jloop].px(),
				// 			// 											jets[jloop].py(),
				// 			// 											jets[jloop].pz(),
				// 			// 											jets[jloop].e());
				// 			// 									glob_fmax = Kzufofmax[input_vec_to_zufos[max_et_candidate_number]];
				// 			// 									//              glob_fmax = stretch_calib(Kzufofmax[input_vec_to_zufos[max_et_candidate_number]],
				// 			// 									//                            "zu", "fmax", l);
				// 			// 									glob_deltaz = myKzufodeltaz[input_vec_to_zufos[max_et_candidate_number]];
				// 			// 								}
				// 			// 							if (check_cuts) cout << "====\n HERE IS PRPH CANDIDATE WITH ETA = " << input_vec[zloop].eta() << "\n====" << endl;
				// 			// 							phot_count++;
				// 			// 							/ *          cout << "---->>>>" << phot_count << " " << Runnr << " " << Eventnr << endl;
				// 			// 								    cout << Siq2el[0] << endl;
				// 			// 								    cout << "zufo et, e: " << input_vec[zloop].et() << ", " << input_vec[zloop].e() << "\n"
				// 			// 								    << "jet et, e: " << jets[jloop].et() << ", "<< jets[jloop].e() << "\n" 
				// 			// 								    << input_vec[zloop].e()/jets[jloop].e()  << " " << input_vec[zloop].et()/jets[jloop].et() << endl;
				// 			// 								    cout << Zufoeemc[input_vec_to_zufos[zloop]]/Zufoecal[input_vec_to_zufos[zloop]] << endl;
				// 			// 								    cout << "its type: " << Tufo[input_vec_to_zufos[zloop]][0] << endl;
				// 			// 								    cout << "number of close tracks: " << ntrack << endl;
				// 			// 								    cout << endl;* /
				// 			// 						}

				// 			// 					cout << "prph candidate " << zloop  << " of " << input_vec.size() << ", eta = " << input_vec[zloop].eta() 
				// 			// 						<< ", et_corr = " << photon_corr_et << ": analysed successfuly" << endl;
				// 			// 				}// if jet contains the zufo
				// 			// 			}//loop over zufos
				// 			// 		}//loop over jets
				// 			// }// if (!use_ktjetb_prph)
						

				// 	//Missing compared to prev code - because only for data prob
				// 	//but I added 
				// 	// find jet containing prompt photon
				// 			if (use_ktjetb_prph)
				// 			{
				// 				Int_t candidate_jet_number = -1;
				// 				for(Int_t jloop = 0; jloop < jets.size(); jloop++)
				// 				{
				// 					if(check_cuts) 
				// 					{
				// 						if (nodebugmode) cout << "jet " << jloop << ": et = " << jets[jloop].et() << ", eta = " << jets[jloop].eta() << endl;
				// 						if (nodebugmode) cout << "max_et_cand (zufos) = "<< max_et_candidate_number << ", max_et_cand (input_vec) = " << input_vec_to_zufos[max_et_candidate_number] 
				// 							<< ", nzufos = " << Nzufos << ", knzufos = " << Knzufos << ", input_vec.size() = " << input_vec_to_zufos.size() << " " << input_vec.size() 
				// 							<< ", zufo_e = " << Zufo[max_et_candidate_number][3] << ", input vec e = " << input_vec[input_vec_to_zufos[max_et_candidate_number]].e() << " " << input_vec[input_vec_to_zufos[max_et_candidate_number]].et() << endl;
				// 						if (nodebugmode) cout << "[][][].et = " << (input_vec[input_vec_to_zufos[max_et_candidate_number]]).et() << endl;
				// 						if (nodebugmode) cout << "index_phot_vector = " << index_phot_vector << endl;
				// 					}

				// 					if( jets[jloop].contains(input_vec[index_phot_vector]) )
				// 					{
				// 						candidate_jet_number = jloop;	
				// 						//break;						
				// 					}
				// 				}
				// 				if (nodebugmode) cout << "candidate_jet_number: " << candidate_jet_number  << endl;
				// 			}

				// 	cout << "After search for PrPh candidate candidate_jet_number: " << candidate_jet_number<< endl;

				// 	// searching for second jet
				// 	Int_t max_et_accomp_jet_candidate_number = -1;
				// 	Double_t max_et_accomp_jet = -999.;

				// 	// Missed compared to prev code - but we don't need this actually here
				// 	// 	//    wtx *= etaph_reweighting(v_corr_prompt_photon->Eta()); //it is as function of background detector level (we do not need cross section for the BG) warning
				// 	// 	if(!Data && mc_type == "mc_bg_norad") 
				// 	// 	{
				// 	// 		wtx *= etaph_reweighting(v_corr_prompt_photon->Eta(), mc_type) * systAcc - (systAcc - 1); //it is as function of background detector level (we do not need cross section for the BG)      
				// 	// 		wtx *= eph_reweighting(v_corr_prompt_photon->E(), mc_type) * systAcc - (systAcc - 1); //it is as function of background detector level (we do not need cross section for the BG)
				// 	// 		if (!nodebugmode) cout << "eta = " << v_corr_prompt_photon->Eta() << ", wtx_eta = " << etaph_reweighting(v_corr_prompt_photon->Eta(), mc_type) << ", e = " << v_corr_prompt_photon->E() << "wtx_e = " << eph_reweighting(v_corr_prompt_photon->E(), mc_type) << endl;
				// 	// 	}
				// 	// 	if(wtx < 0) 
				// 	// 	{
				// 	// 		cerr << "wtx = " << wtx << endl; 
				// 	// 		exit(-1);
				// 	// 	}
					

				// 	if (take_event && here_is_prph)
				// 	{
				// 		for(Int_t jloop = 0; jloop < jets.size(); jloop++)
				// 		{
				// 			if (check_cuts) cout << "jet " << jloop << ": calling correction with parameters: " << jets[jloop].eta() 
				// 					<< ", " << jets[jloop].et() << ", " << period << ", " << mc_type << endl;
							
				// 			//Missing compared to prev code
				// 			//Double_t et_jet_corr = jet_en_corr_1stana(jets[jloop].eta(), jets[jloop].et(), period, mc_type);
				// 			Double_t et_jet_corr = jet_en_corr(jets[jloop].eta(), jets[jloop].et(), period, mc_type);
							
				// 			if (check_cuts) cout << "jet " << jloop << " of " << jets.size() 
				// 					<< ", eta = " << jets[jloop].eta() 
				// 					<< ", et_uncorr = " << jets[jloop].et() << ", et_corr = " << et_jet_corr << endl;
							
				// 			take_jet = kTRUE;
				// 			//Missing compared to prev code: if((jloop == candidate_jet_number)) 
				// 			if ((!use_ktjetb_prph) && (jloop == candidate_jet_number)) 
				// 			{
				// 				take_jet = kFALSE;
				// 				if (check_cuts) cout << "jet with eta " << jets[jloop].eta() << " rejected because it is prompt photon jet " << endl;
				// 			}
				// 			if (et_jet_corr < ET_JET_CUT) 
				// 			{
				// 				take_jet = kFALSE;
				// 				if (check_cuts) cout << "jet with eta " << jets[jloop].eta() << " rejected by et_jet_corr = "<< et_jet_corr << endl;
				// 			}
				// 			if (jets[jloop].eta() < -1.5 || jets[jloop].eta() > 1.8) 
				// 			{
				// 				take_jet = kFALSE;
				// 				if (check_cuts) cout << "jet with eta " << jets[jloop].eta() << " rejected by jets[" << jloop << "].eta() = "<< jets[jloop].eta() << endl;
				// 			}
							
				// 			if (take_jet)
				// 			{
				// 				// cout << "jet " << jloop << ": " << jets[jloop].eta() << " " << jets[jloop].et() << endl;
				// 				here_is_jet = kTRUE;
				// 				// max_et_accomp_jet_candidate_number = jloop;
				// 				if (et_jet_corr > max_et_accomp_jet)
				// 				{
				// 					max_et_accomp_jet = et_jet_corr;
				// 					max_et_accomp_jet_candidate_number = jloop;
				// 					Double_t Px = jets[jloop].px();
				// 					Double_t Py = jets[jloop].py();
				// 					Double_t Pz = jets[jloop].pz();
				// 					Double_t E = jets[jloop].e();
				// 					v_accomp_uncorr_jet->SetPxPyPzE(Px, Py, Pz, E);
				// 					/*Missing compared to prev code:
				// 					Px = jet_en_corr_1stana(jets[jloop].eta(), jets[jloop].px(), period, mc_type);
				// 					Py = jet_en_corr_1stana(jets[jloop].eta(), jets[jloop].py(), period, mc_type);
				// 					Pz = jet_en_corr_1stana(jets[jloop].eta(), jets[jloop].pz(), period, mc_type);
				// 					E = jet_en_corr_1stana(jets[jloop].eta(), jets[jloop].e(), period, mc_type);
				// 					*/
				// 					Px = jet_en_corr(jets[jloop].eta(), jets[jloop].px(), period, mc_type);
				// 					Py = jet_en_corr(jets[jloop].eta(), jets[jloop].py(), period, mc_type);
				// 					Pz = jet_en_corr(jets[jloop].eta(), jets[jloop].pz(), period, mc_type);
				// 					E = jet_en_corr(jets[jloop].eta(), jets[jloop].e(), period, mc_type);
				// 					v_accomp_corr_jet->SetPxPyPzE(Px, Py, Pz, E);
				// 				}

				// 				if (study_jets) 
				// 				{
				// 					Double_t n_subzufos = 0;
				// 					for(Int_t zloop = 0; zloop < input_vec.size(); zloop++)
				// 						if (jets[max_et_accomp_jet_candidate_number].contains(input_vec[zloop]) )
				// 							n_subzufos++;
				// 					hist.h2d_jet_npart_et->Fill(n_subzufos, et_jet_corr, wtx);
				// 				}//if study_jets

				// 			}//if take_jet
				// 		}// for jloop over jets
				// 	}
					

				// 	cout << "max_et_accomp_jet_candidate_number: " << max_et_accomp_jet_candidate_number << endl;
	 
				// 	// fill inclusive prompt photon histograms
				// 		if (take_event && here_is_prph) 
				// 		{
				// 			hist.incl_phot_fmax->Fill(glob_fmax, wtx);
				// 			hist.incl_phot_deltaz->Fill(glob_deltaz, wtx);
				// 			hist.incl_phot_elphot_dphi->Fill(delta_phi(v_corr_prompt_photon->Phi(), v_electron.Phi()), wtx);
				// 			hist.incl_phot_elphot_deta->Fill(v_corr_prompt_photon->Eta() - v_electron.Eta(), wtx);
				// 			hist.incl_phot_phi->Fill(v_corr_prompt_photon->Phi(), wtx);
				// 			hist.incl_phot_eta->Fill(v_corr_prompt_photon->Eta(), wtx);
				// 		}
				
				// 	// fill prompt photon + jet
				// 	Double_t hardest_jet_et_corr = 0;
				// 	Double_t hardest_jet_et = 0;
				// 	Double_t hardest_jet_eta = -999.;
				// 	Double_t hardest_jet_phi = -999.;
				// 	if (here_is_jet && !inclusive_prph_analysis) {
				// 		hardest_jet_et = jets[max_et_accomp_jet_candidate_number].et();
				// 		hardest_jet_eta = jets[max_et_accomp_jet_candidate_number].eta();
				// 		hardest_jet_phi = jets[max_et_accomp_jet_candidate_number].phi();
				// 		//Missing compared to prev code: hardest_jet_et_corr = jet_en_corr_1stana(hardest_jet_eta, hardest_jet_et, period, mc_type);
				// 		hardest_jet_et_corr = jet_en_corr(hardest_jet_eta,
				// 				hardest_jet_et, period, mc_type);
				// 	}
				// 	if (inclusive_prph_analysis)
				// 		here_is_jet = kTRUE;
					
				// 	if(check_cuts && !nodebugmode) cout << "take_event = " << take_event << ", here_is_jet = " << here_is_jet << ", here_is_prph = " << here_is_prph << endl;
		
				// 	// Prev
				// 	// 	x_gamma = v_uncorr_prompt_photon->E() - v_uncorr_prompt_photon->Pz();
				// 	// 	for(Int_t zloop=0; zloop<input_vec.size(); zloop++)
				// 	// 		if ( jets[max_et_accomp_jet_candidate_number].contains(input_vec[zloop]) )
				// 	// 			x_gamma += input_vec[zloop].e() - input_vec[zloop].pz();
						
				// 	// 	x_gamma /= (2. * E_e * yjb);
				// 	// 	if (x_gamma > 0.75) take_event = kFALSE;
					

				// 	if (take_event && here_is_jet && here_is_prph) 
				// 	{

				// 		cout << "max_et_accomp_jet_candidate_number: " << max_et_accomp_jet_candidate_number << endl;
	 
				// 		//
				// 		// x_gamma
				// 		//
				// 		Double_t x_pomeron = (v_uncorr_prompt_photon->E() + v_uncorr_prompt_photon->Pz()
				// 							+ jets[max_et_accomp_jet_candidate_number].e() + jets[max_et_accomp_jet_candidate_number].pz())
				// 							/ (2. * E_p);
				// 		x_gamma = (v_uncorr_prompt_photon->E() - v_uncorr_prompt_photon->Pz() 
				// 				+ jets[max_et_accomp_jet_candidate_number].e() - jets[max_et_accomp_jet_candidate_number].pz())
				// 		 		/ (2. * E_e * yjb);

				// 		if (!nodebugmode) cout << "v_uncorr_prompt_photon->E() : " << v_uncorr_prompt_photon->E() << endl;
				// 		if (!nodebugmode) cout << "v_uncorr_prompt_photon->Pz() : " << v_uncorr_prompt_photon->Pz() << endl;
				// 		if (!nodebugmode) cout << "E_e : " << E_e << endl;
				// 		if (!nodebugmode) cout << "yjb : " << yjb << endl;
				// 		if (!nodebugmode) cout << "v_uncorr_prompt_photon->E() - v_uncorr_prompt_photon->Pz(): " << Runnr << " " << Eventnr << " " << x_gamma<<endl;
				// 		Double_t empz_jet = jets[max_et_accomp_jet_candidate_number].e() - jets[max_et_accomp_jet_candidate_number].pz();
						
				// 		Double_t empz_jet_plus = jets[max_et_accomp_jet_candidate_number].e() + jets[max_et_accomp_jet_candidate_number].pz();
						
				// 		if (!nodebugmode) cout<<"JETS[]"<<endl;
				// 		if (!nodebugmode) cout<<"!    jets[max_et_accomp_jet_candidate_number].getNConstituents() = "<< jets[max_et_accomp_jet_candidate_number].getNConstituents()<<endl;
				// 		if (!nodebugmode) cout<<"!    max_et_accomp_jet_candidate_number = "<< max_et_accomp_jet_candidate_number <<endl;
				// 		if (!nodebugmode && (x_pomeron>=1 )) cout << "x_pomeron: " << x_pomeron <<endl;
						
				// 		if (check_cuts && !nodebugmode) cout << "event accepted" << endl;
				// 		hist.eta_gamma_tef[1]->Fill(v_corr_prompt_photon->Eta(), wtx);
				// 		hist.e_jet_tef[1]->Fill(hardest_jet_et_corr, wtx);
				// 		hist.q2_tef[1]->Fill(Siq2el[0], wtx);
				// 		hist.e_electron_tef[1]->Fill(Siecorr[0][2], wtx);
				// 		hist.e_gamma_tef[1]->Fill(v_corr_prompt_photon->Et(), wtx);
				// 		hist.theta_electron_tef[1]->Fill(Sith[0], wtx);

				// 		if (take_event_trig)
				// 		{
				// 			//cout<< "THIS EVENT PASSED THROUGH THE END"<< endl;
				// 			event_list->Enter(entry);
				// 			hist.eta_gamma_tef[0]->Fill(v_corr_prompt_photon->Eta(), wtx);
				// 			hist.e_jet_tef[0]->Fill(hardest_jet_et_corr, wtx);
				// 			hist.q2_tef[0]->Fill(Siq2el[0], wtx);
				// 			hist.e_electron_tef[0]->Fill(Siecorr[0][2], wtx);
				// 			hist.e_gamma_tef[0]->Fill(v_corr_prompt_photon->Et(), wtx);
				// 			hist.theta_electron_tef[0]->Fill(Sith[0], wtx);

				// 			fill_trigger_bits(hist.dstw, Dstw);
				// 			fill_trigger_bits_general(hist.fltw, Fltw, 2);
				// 			fill_trigger_bits_general(hist.sltw, Sltw, 6);
				// 			fill_trigger_bits_general(hist.tltw, Tltw, 15);

				// 			hist.fmax->Fill(glob_fmax, wtx);
				// 			hist.deltaz->Fill(glob_deltaz, wtx);
				// 			hist.fmax_deltaz->Fill(glob_fmax, glob_deltaz, wtx);
				// 			hist.prph_energy->Fill(v_corr_prompt_photon->Et(), wtx);
				// 			hist.prph_eta->Fill(v_corr_prompt_photon->Eta(), wtx);
				// 			hist.prph_phi->Fill(v_corr_prompt_photon->Phi(), wtx);
				// 			hist.prph_jet_energy_frac->Fill(jet_energy_frac, wtx);
				// 			hist.prph_cell_energy_frac->Fill(cell_energy_frac, wtx);

				// 			//for cross sections:
				// 				hist.det_cross_et->Fill(v_corr_prompt_photon->Et(), wtx);
				// 				hist.det_cross_eta->Fill(v_corr_prompt_photon->Eta(), wtx);
				// 				hist.det_cross_Q2->Fill(Siq2el[0], wtx);
				// 				hist.det_cross_x->Fill(Sixel[0], wtx);
				// 				hist.det_cross_et_jet->Fill(hardest_jet_et_corr, wtx);
				// 				hist.det_cross_eta_jet->Fill(hardest_jet_eta, wtx);
				// 					Double_t temp_dphi = delta_phi(hardest_jet_phi, v_corr_prompt_photon->Phi()) * 180.0/TMath::Pi();
				// 					Double_t temp_deta = hardest_jet_eta - v_corr_prompt_photon->Eta();
				// 					Double_t temp_dphi_e_ph = delta_phi(Siph[sinistra_electron_number], v_corr_prompt_photon->Phi()) * 180.0/TMath::Pi();
				// 					Double_t temp_deta_e_ph = -TMath::Log(TMath::Tan(Sith[sinistra_electron_number]/2.)) - v_corr_prompt_photon->Eta();
				// 					hist.det_cross_xgamma->Fill(x_gamma, wtx); 
				// 					hist.det_cross_xp->Fill(x_pomeron, wtx); 
				// 					hist.det_cross_dphi->Fill( temp_dphi, wtx);
				// 					hist.det_cross_deta->Fill( temp_deta, wtx);
				// 					hist.det_cross_dphi_e_ph->Fill(temp_dphi_e_ph, wtx); 
				// 					hist.det_cross_deta_e_ph->Fill(temp_deta_e_ph, wtx); 

				// 				hist.prof_det_cross_et->Fill(v_corr_prompt_photon->Et(), v_corr_prompt_photon->Et(), wtx);
				// 				hist.prof_det_cross_eta->Fill(v_corr_prompt_photon->Eta(), v_corr_prompt_photon->Eta(), wtx);
				// 				hist.prof_det_cross_Q2->Fill(Siq2el[0], Siq2el[0], wtx);
				// 				hist.prof_det_cross_x->Fill(Sixel[0], Sixel[0], wtx);
				// 				hist.prof_det_cross_et_jet->Fill(hardest_jet_et_corr, hardest_jet_et_corr, wtx);
				// 				hist.prof_det_cross_eta_jet->Fill(hardest_jet_eta, hardest_jet_eta, wtx);
				// 					hist.prof_det_cross_xgamma->Fill(x_gamma, x_gamma, wtx); 
				// 					hist.prof_det_cross_xp->Fill(x_pomeron, x_pomeron, wtx); 
				// 					hist.prof_det_cross_dphi->Fill(temp_dphi, temp_dphi, wtx);
				// 					hist.prof_det_cross_deta->Fill(temp_deta, temp_deta, wtx);
				// 					hist.prof_det_cross_dphi_e_ph->Fill(temp_dphi_e_ph, temp_dphi_e_ph, wtx); 
				// 					hist.prof_det_cross_deta_e_ph->Fill(temp_deta_e_ph, temp_deta_e_ph, wtx); 

				// 				det_cross_sec_bin_et = hist.det_cross_et->FindBin(v_corr_prompt_photon->Et());
				// 				det_cross_sec_bin_eta = hist.det_cross_eta->FindBin(v_corr_prompt_photon->Eta());
				// 				det_cross_sec_bin_q2 = hist.det_cross_Q2->FindBin(Siq2el[0]);
				// 				det_cross_sec_bin_x = hist.det_cross_x->FindBin(Sixel[0]);
				// 				det_cross_sec_bin_et_jet = hist.det_cross_et_jet->FindBin(hardest_jet_et_corr);
				// 				det_cross_sec_bin_eta_jet = hist.det_cross_eta_jet->FindBin(hardest_jet_eta);
				// 					det_cross_sec_bin_xgamma = hist.det_cross_xgamma->FindBin(x_gamma); 
				// 					det_cross_sec_bin_xp = hist.det_cross_xp->FindBin(x_pomeron); 
				// 					det_cross_sec_bin_dphi = hist.det_cross_dphi->FindBin(temp_dphi);
				// 					det_cross_sec_bin_deta = hist.det_cross_deta->FindBin(temp_deta);
				// 					det_cross_sec_bin_dphi_e_ph = hist.det_cross_dphi_e_ph->FindBin(temp_dphi_e_ph); 
				// 					det_cross_sec_bin_deta_e_ph = hist.det_cross_deta_e_ph->FindBin(temp_deta_e_ph); 
								
				// 			//check for bin migrations
				// 				for(Int_t i = 0; i < number_etbins; i++)
				// 				{
				// 					if (v_corr_prompt_photon->Et() > et_bin[i] && v_corr_prompt_photon->Et() < et_bin[i+1]) {
				// 						hist.fmax_et[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_et[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i < number_Q2bins; i++)
				// 				{
				// 					if (Siq2el[0] > Q2_bin[i] && Siq2el[0] < Q2_bin[i+1]) {
				// 						hist.fmax_q2[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_q2[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i < number_xbins; i++)
				// 				{
				// 					if (Sixel[0] > x_bin[i] && Sixel[0] < x_bin[i+1]) {
				// 						hist.fmax_x[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_x[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i < number_etabins; i++)
				// 				{
				// 					if (v_corr_prompt_photon->Eta() > eta_bin_crosssec[i] && v_corr_prompt_photon->Eta() < eta_bin_crosssec[i+1]) {
				// 						hist.fmax_eta[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_eta[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i < number_et_jetbins; i++)
				// 				{
				// 					if (hardest_jet_et_corr > et_jet_bin[i] && 
				// 							hardest_jet_et_corr < et_jet_bin[i+1]) {
				// 						hist.fmax_et_jet[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_et_jet[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i < number_eta_jetbins; i++)
				// 				{
				// 					if (hardest_jet_eta > eta_jet_bin[i] && 
				// 							hardest_jet_eta < eta_jet_bin[i+1]) {
				// 						hist.fmax_eta_jet[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_eta_jet[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				Int_t flag = 0;
				// 				for(Int_t i = 0; i<number_xgamma_bins; i++)
				// 				{
				// 					if(x_gamma > xgamma_bin[i] && x_gamma < xgamma_bin[i+1]) 
				// 					{
				// 						hist.fmax_xgamma[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_xgamma[i]->Fill(glob_deltaz, wtx);
				// 						flag++;
				// 					}
				// 				}
				// 				if (flag != 0) missed++;
				// 				for(Int_t i = 0; i<number_xp_bins; i++)
				// 				{
				// 					if(x_pomeron > xp_bin[i] && x_pomeron < xp_bin[i+1]) 
				// 					{
				// 						hist.fmax_xp[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_xp[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i<number_dphi_bins; i++)
				// 				{
				// 					if(temp_dphi > dphi_bin[i] && temp_dphi < dphi_bin[i+1]) 
				// 					{
				// 						hist.fmax_dphi[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_dphi[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i<number_deta_bins; i++)
				// 				{
				// 					if(temp_deta > deta_bin[i] && temp_deta < deta_bin[i+1]) 
				// 					{
				// 						hist.fmax_deta[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_deta[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i<number_dphi_e_ph_bins; i++)
				// 				{
				// 					if(temp_dphi_e_ph > dphi_e_ph_bin[i] && temp_dphi_e_ph < dphi_e_ph_bin[i+1]) 
				// 					{
				// 						hist.fmax_dphi_e_ph[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_dphi_e_ph[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}
				// 				for(Int_t i = 0; i<number_deta_e_ph_bins; i++)
				// 				{
				// 					if(temp_deta_e_ph > deta_e_ph_bin[i] && temp_deta_e_ph < deta_e_ph_bin[i+1]) 
				// 					{
				// 						hist.fmax_deta_e_ph[i]->Fill(glob_fmax, wtx);
				// 						hist.deltaz_deta_e_ph[i]->Fill(glob_deltaz, wtx);
				// 					}
				// 				}

				// 			hist.x_gamma->Fill(x_gamma, wtx);
				// 			hist.x_pomeron->Fill(x_pomeron, wtx);
				// 			//
				// 			// Sinistra candidates and photon candidate
				// 			//
				// 			Double_t photon_eta = v_corr_prompt_photon->Eta();
				// 			Double_t photon_phi = v_corr_prompt_photon->Phi();
				// 			Double_t min_sin_phot_dr = 999.;
				// 			Int_t closest_sincand = -1;
				// 			for(Int_t sloop=0; sloop<Sincand; sloop++)
				// 			{
				// 				Double_t dphi = delta_phi(photon_phi, Siph[sloop]);
				// 				Double_t deta = photon_eta+TMath::Log(TMath::Tan(Sith[sloop]/2.));
				// 				Double_t dr = TMath::Sqrt(dphi*dphi + deta*deta);
				// 				hist.phjet_photon_sincand_dr[sloop]->Fill(dr, wtx);
				// 				//          cout << "sloop = " << sloop << ": " << dphi << " " << deta << " "  << dr << endl;
				// 				if (dr < min_sin_phot_dr) {

				// 					min_sin_phot_dr = dr;
				// 					closest_sincand = sloop;
				// 				}
				// 			}
				// 			hist.phjet2d_photon_sincand_dr_nr->Fill((Double_t)closest_sincand, min_sin_phot_dr, wtx);
				// 			//          cout << Runnr << " " << Eventnr << " " << closest_sincand << " " << min_sin_phot_dr << endl;
				// 			//
				// 			//  Fill histograms
				// 			//
				// 			Double_t photon_jet_deta = hardest_jet_eta - v_corr_prompt_photon->Eta();
				// 			Double_t photon_jet_dphi = delta_phi(hardest_jet_phi, photon_phi);

				// 			Double_t photon_electron_dphi_true = delta_phi(Siph[sinistra_electron_number], photon_phi);
				// 			Double_t photon_electron_deta_true = -TMath::Log(TMath::Tan(Sith[sinistra_electron_number]/2.)) - v_corr_prompt_photon->Eta();
				// 			Double_t photon_electron_deta = -TMath::Log(TMath::Tan(Sith[0]/2.)) - v_corr_prompt_photon->Eta();
				// 			Double_t photon_electron_dphi = delta_phi(Siph[0], photon_phi);

				// 			Double_t jet_electron_deta = -TMath::Log(TMath::Tan(Sith[0]/2.)) - hardest_jet_eta;
				// 			Double_t jet_electron_dphi = delta_phi(hardest_jet_phi, Siph[0]);

				// 			hist.comp_jet_e->Fill(max_et_accomp_jet, wtx);
				// 			hist.comp_prph_e->Fill(max_et_candidate, wtx);
				// 			hist.comp_Q2->Fill(Siq2el[0], wtx);
				// 			hist.comp_x->Fill(Sixel[0], wtx);
				// 			hist.comp_y->Fill(Siyel[0], wtx);
				// 			if (Sixel[0]*Siyel[0] != 0)
				// 				hist.comp_s->Fill(Siq2el[0]/(Sixel[0]*Siyel[0]));
				// 			hist.phjet_electron_e->Fill(Siecorr[0][2], wtx);
				// 			/*          Int_t MCflag = 0;
				// 				    Float_t m_Siecorr[4][3];
				// 				    if (Runnr == 1) MCflag = 1;
				// 				    m_Siecorr[0][2] = BCR_CalibTools::ecorr_newscheme(Sicalene[0], Six0[0],
				// 				    Sicalpos[0][0],Sicalpos[0][1],
				// 				    MCflag);
				// 				    hist.phjet_electron_e_corr_burkard->Fill(m_Siecorr[0][2], wtx);*/
				// 			hist.phjet_electron_theta->Fill(Sith[0]*180.0/TMath::Pi(), wtx);
				// 			hist.phjet_electron_phi->Fill(Siph[0], wtx);
				// 			hist.phjet_electron_probability->Fill(Siprob[0], wtx);
				// 			hist.phjet_dphi_deta->Fill(photon_jet_dphi, photon_jet_deta, wtx);
				// 			hist.phjet_jet_eta->Fill(hardest_jet_eta, wtx);
				// 			hist.phjet_jet_phi->Fill(hardest_jet_phi, wtx);

				// 			hist.phjet_deta->Fill(photon_jet_deta, wtx);
				// 			hist.phjet_dphi->Fill(photon_jet_dphi * 180.0/TMath::Pi(), wtx);
				// 			hist.phjet_deta_el_ph_true->Fill(photon_electron_deta_true, wtx);
				// 			hist.phjet_dphi_el_ph_true->Fill(photon_electron_dphi_true * 180.0/TMath::Pi(), wtx);

				// 			hist.phjet_deta_el_ph->Fill(photon_electron_deta, wtx);
				// 			hist.phjet_dphi_el_ph->Fill(photon_electron_dphi * 180.0/TMath::Pi(), wtx);
				// 			hist.phjet_deta_el_jet->Fill(jet_electron_deta, wtx);
				// 			hist.phjet_dphi_el_jet->Fill(jet_electron_dphi, wtx);
				// 			hist.et_jet_photon_ratio->Fill(max_et_accomp_jet/max_et_candidate, wtx);
				// 			hist.h2d_et_jet_et_photon->Fill(max_et_accomp_jet, max_et_candidate, wtx);
				// 			Double_t et_3system_ev_ratio = (max_et_accomp_jet + max_et_candidate + Siecorr[0][2]*TMath::Sin(Sith[0]))/et_of_event;
				// 			hist.et_3system_event_ratio->Fill(et_3system_ev_ratio, wtx);
				// 			hist.h2d_phijet_phigamma->Fill(hardest_jet_phi, v_corr_prompt_photon->Phi(), wtx);
				// 			hist.h2d_dphijet_dphigamma->Fill(hardest_jet_phi - 
				// 					Siph[0], delta_phi(v_corr_prompt_photon->Phi(), Siph[0]), wtx);
				// 			Double_t delta_phi_jet_el = delta_phi(hardest_jet_phi, Siph[0]);
				// 			Double_t delta_phi_photon_el = delta_phi(v_corr_prompt_photon->Phi(), Siph[0]);
				// 			hist.h2d_dphijet_dphigamma_sharp->Fill(delta_phi_jet_el, delta_phi_photon_el, wtx);
				// 			hist.phjet_cal_empz->Fill(Cal_empz, wtx);
				// 			hist.phjet_zufo_empz->Fill(Empz, wtx);
				// 			hist.phjet_zvtx->Fill(Zvtx, wtx);

				// 			//
				// 			// histograms for Peter+Misha reweighting procedure
				// 			//
				// 			if (glob_deltaz > 0.35)
				// 				hist.phjet_q2_dz_gt_035->Fill(Siq2el[0], wtx);
				// 			if (glob_deltaz < 0.35)
				// 				hist.phjet_q2_dz_lt_035->Fill(Siq2el[0], wtx);

				// 			hist.phjet_q2_dz_full->Fill(Siq2el[0], wtx);

				// 			if (!Data) 
				// 			{
				// 				hist.mccorel_q2[0]->Fill(Siq2el[0] - Mc_q2, Mc_q2, wtx);
				// 				hist.mccorel_q2[1]->Fill(Siq2da[0] - Mc_q2, Mc_q2, wtx);
				// 				hist.mccorel_q2[2]->Fill(Siq2jb[0] - Mc_q2, Mc_q2, wtx);

				// 				hist.mccorel_x[0]->Fill(Sixel[0] - Mc_x, Mc_x, wtx);
				// 				hist.mccorel_x[1]->Fill(Sixda[0] - Mc_x, Mc_x, wtx);
				// 				hist.mccorel_x[2]->Fill(Sixjb[0] - Mc_x, Mc_x, wtx);

				// 				hist.mccorel_y[0]->Fill(Siyel[0] - Mc_y, Mc_y, wtx);
				// 				hist.mccorel_y[1]->Fill(Siyda[0] - Mc_y, Mc_y, wtx);
				// 				hist.mccorel_y[2]->Fill(Siyjb[0] - Mc_y, Mc_y, wtx);
				// 			}

				// 			selectedEventsCount++;

				// 			//printout interesting events

				// 			if (print_interesting_events && et_3system_ev_ratio > 0.95) 
				// 			{
				// 				cout << endl;
				// 				cout << "============ INTERESTING EVENT ================" << endl;
				// 				cout << "========= " << Runnr << " " << Eventnr << " ================" << endl;
				// 				cout << " seems to be clean: et_3system_ev_ratio = " << et_3system_ev_ratio << endl;
				// 				cout << "===============================================" << endl << endl;        
				// 			}
				// 		}//if take_event_trig


				// 	}
				// }
				// //if take_event for jet searching
			
			//Some output and some hists fill for Had and Part selections
				
				if (!Data && SelectHadronLevel(take_event && here_is_jet && here_is_prph && take_event_trig)) 
				{
					cout << "SelectHadronLevel selected" << endl;
					take_had_event = kTRUE;
					num_had += 1;
					if (take_event && here_is_jet && here_is_prph && take_event_trig) 
					{
						hist.h2d_phot_et_true_det->Fill(v_true_prompt_photon->Et(), v_corr_prompt_photon->Et());
						hist.h2d_phot_eta_true_det->Fill(v_true_prompt_photon->Eta(), v_corr_prompt_photon->Eta());
						hist.h2d_phot_phi_true_det->Fill(v_true_prompt_photon->Phi(), v_corr_prompt_photon->Phi());
						hist.h2d_corr_accjet_et_true_det->Fill(v_true_acc_jet->Et(), v_accomp_corr_jet->Et());
						hist.h2d_corr_accjet_eta_true_det->Fill(v_true_acc_jet->Eta(), v_accomp_corr_jet->Eta());
						hist.h2d_corr_accjet_phi_true_det->Fill(v_true_acc_jet->Phi(), v_accomp_corr_jet->Phi());
						hist.h2d_uncorr_accjet_et_true_det->Fill(v_true_acc_jet->Et(), v_accomp_uncorr_jet->Et());
						hist.h2d_uncorr_accjet_eta_true_det->Fill(v_true_acc_jet->Eta(), v_accomp_uncorr_jet->Eta());
						hist.h2d_uncorr_accjet_phi_true_det->Fill(v_true_acc_jet->Phi(), v_accomp_uncorr_jet->Phi());
						//      hist.h2d_phot_en_true_det->Fill();
					}
					//  event_list->Enter(entry);
				}
				if (!Data && SelectPartonLevel(take_event && here_is_jet && here_is_prph && take_event_trig, take_had_event)) 
				{
					num_part += 1;
					cout << "SelectPartonLevel" << endl;
					take_part_event = kTRUE;
					cout << "Parton level event " << Eventnr << " is selected and N of hadron level jets is " << hadron_level_jets.size() << endl;
					// for(Int_t i = 0; i < hadron_level_jets.size(); i++) 
					// 	cout << "hadr. level jet #" << i << ": et = " << hadron_level_jets[i].et() << ", eta = " << hadron_level_jets[i].eta() << endl;
					// cout << "prompt photon jet is jet #" << hadron_level_jet_cont_photon_index << " with e ratio = " << hadron_level_ephoton_over_ejet << endl;
					cout << " === THE END === " << endl;
				}
			if (!en_mom_conservation) 
			{
				cout <<"EM NOT PRESERVED"<< endl;
				continue;
			}
			//part nohad || part_bin != had_bin
				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_et->FindBin(part_et) != hist.had_cross_et->FindBin(had_et)) ) 
					hist.part_nohad_cross_et->Fill(part_et);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_eta->FindBin(part_eta) != hist.had_cross_eta->FindBin(had_eta)) ) 
					hist.part_nohad_cross_eta->Fill(part_eta);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_Q2->FindBin(part_Q2) != hist.had_cross_Q2->FindBin(had_Q2)) ) 
					hist.part_nohad_cross_Q2->Fill(part_Q2);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_x->FindBin(part_x) != hist.had_cross_x->FindBin(had_x)) ) 
					hist.part_nohad_cross_x->Fill(part_x);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_et_jet->FindBin(part_et_jet) != hist.had_cross_et_jet->FindBin(had_et_jet)) ) 
				{
					hist.part_nohad_cross_et_jet->Fill(part_et_jet);
					hist.part_nohad_cross_et_jet2->Fill(part_et_jet);
				}

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_eta_jet->FindBin(part_eta_jet) != hist.had_cross_eta_jet->FindBin(had_eta_jet)) ) 
					hist.part_nohad_cross_eta_jet->Fill(part_eta_jet);
				//new vars
				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_xgamma->FindBin(part_xgamma) != hist.had_cross_xgamma->FindBin(had_xgamma)) ) 
					hist.part_nohad_cross_xgamma->Fill(part_xgamma);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_xp->FindBin(part_xp) != hist.had_cross_xp->FindBin(had_xp)) ) 
					hist.part_nohad_cross_xp->Fill(part_xp);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_dphi->FindBin(part_dphi) != hist.had_cross_dphi->FindBin(had_dphi)) ) 
					hist.part_nohad_cross_dphi->Fill(part_dphi);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_deta->FindBin(part_deta) != hist.had_cross_deta->FindBin(had_deta)) ) 
					hist.part_nohad_cross_deta->Fill(part_deta);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_dphi_e_ph->FindBin(part_dphi_e_ph) != hist.had_cross_dphi_e_ph->FindBin(had_dphi_e_ph)) ) 
				{
					hist.part_nohad_cross_dphi_e_ph->Fill(part_dphi_e_ph);
				}

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_deta_e_ph->FindBin(part_deta_e_ph) != hist.had_cross_deta_e_ph->FindBin(had_deta_e_ph)) ) 
					hist.part_nohad_cross_deta_e_ph->Fill(part_deta_e_ph);

			//had nopart || part_bin != had_bin
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_et->FindBin(part_et) != hist.had_cross_et->FindBin(had_et)) ) 
					hist.had_nopart_cross_et->Fill(had_et);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_eta->FindBin(part_eta) != hist.had_cross_eta->FindBin(had_eta)) ) 
					hist.had_nopart_cross_eta->Fill(had_eta);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_Q2->FindBin(part_Q2) != hist.had_cross_Q2->FindBin(had_Q2)) ) 
					hist.had_nopart_cross_Q2->Fill(had_Q2);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_x->FindBin(part_x) != hist.had_cross_x->FindBin(had_x)) ) 
					hist.had_nopart_cross_x->Fill(had_x);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_et_jet->FindBin(part_et_jet) != hist.had_cross_et_jet->FindBin(had_et_jet)) ) 
				{
					hist.had_nopart_cross_et_jet->Fill(had_et_jet);
					hist.had_nopart_cross_et_jet2->Fill(had_et_jet);
				}
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_eta_jet->FindBin(part_eta_jet) != hist.had_cross_eta_jet->FindBin(had_eta_jet)) ) 
					hist.had_nopart_cross_eta_jet->Fill(had_eta_jet);
				//new vars

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_xgamma->FindBin(part_xgamma) != hist.had_cross_xgamma->FindBin(had_xgamma)) ) 
					hist.had_nopart_cross_xgamma->Fill(had_xgamma);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_xp->FindBin(part_xp) != hist.had_cross_xp->FindBin(had_xp)) ) 
					hist.had_nopart_cross_xp->Fill(had_xp);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_dphi->FindBin(part_dphi) != hist.had_cross_dphi->FindBin(had_dphi)) ) 
					hist.had_nopart_cross_dphi->Fill(had_dphi);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_deta->FindBin(part_deta) != hist.had_cross_deta->FindBin(had_deta)) ) 
					hist.had_nopart_cross_deta->Fill(had_deta);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_dphi_e_ph->FindBin(part_dphi_e_ph) != hist.had_cross_dphi_e_ph->FindBin(had_dphi_e_ph)) ) 
				{
					hist.had_nopart_cross_dphi_e_ph->Fill(had_dphi_e_ph);
				}
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && hist.part_cross_deta_e_ph->FindBin(part_deta_e_ph) != hist.had_cross_deta_e_ph->FindBin(had_deta_e_ph)) ) 
					hist.had_nopart_cross_deta_e_ph->Fill(had_deta_e_ph);

		}// if exact event
		hadron_level_jets.clear();
	}// for entry over entries

	Terminate();
	cout << "num_had = " << num_had << " num_part = " << num_part << endl; 
	return kTRUE;
}

void selector::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.

}


	/*
		void sortEt(TObjArray *arrayNotSorted, TObjArray *arraySorted)
		  {
		//  Double_t maxEt = 0.;// (static_cast<TLorentzVector*> (arrayNotSorted->At(0)))->Et();
		Double_t et[100];

		Int_t sortedIndex[100] = {0};
		Int_t counter = 0;
		for(Int_t i = 0; i < arrayNotSorted->GetEntries(); i++)
		{
		et[i] = (static_cast<TLorentzVector*> (arrayNotSorted->At(i)))->Et();
		}
		while(counter<arrayNotSorted->GetEntries())
		{
		Int_t max;
		Int_t N = arrayNotSorted->GetEntries();
		max = findMaxEt(sortedIndex, et, N);
		Double_t px_jet = (static_cast<TLorentzVector*>(arrayNotSorted->At(max)))->Px();
		Double_t py_jet = (static_cast<TLorentzVector*>(arrayNotSorted->At(max)))->Py();
		Double_t pz_jet = (static_cast<TLorentzVector*>(arrayNotSorted->At(max)))->Pz();
		Double_t e_jet = (static_cast<TLorentzVector*>(arrayNotSorted->At(max)))->E();
		arraySorted->Add(new TLorentzVector(px_jet, py_jet, pz_jet, e_jet) );
		sortedIndex[max] = 1;
		counter++;
		}
		}

		Int_t findMaxEt(Int_t sorted[], Double_t et_jets[], Int_t N)
		{
		Double_t et;
		Double_t maxEt = 0;
		Int_t sort = -1;
		for(Int_t j=0; j<N; j++)
		{
		et = et_jets[j];
		if (et>maxEt && sorted[j]==0)
		{
		maxEt = et;
		sort = j;
		}
		}
		return sort;
		}
	*/
	//#endif
