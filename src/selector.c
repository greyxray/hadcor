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
		bool test_on_entry = false;
		int test_entry = 4;
		check_cuts = kTRUE;
		nodebugmode = kFALSE;
		Bool_t debugcontinue = kTRUE;

  
	Int_t Runnr_prev = 0;
	Int_t missed = 0;
	int num_had(0), num_part(0);
	check_en_mom_conservation = true;
	check_en_mom_conservation_on_parton_level = true;
	part_lev_from_fmckin2 = true;
	had_lev_from_fmckin2 = true;

	for( entry = 0; entry < nentries - 1  && debugcontinue; entry++)
	{
		wtx = 1;
		//if (entry > 10) exit(1);
		en_mom_conservation = true;

		// To test specific entry
			if (test_on_entry)
			{
				if (entry < 4 ) continue;
				else if (entry > 4) 
				{
					cout << "not interested in restof the entries " << endl;
					exit(0);
				}
			}
		 
		fChain->GetEntry(entry);

		//cout << "maybe This is ev " << entry << " " <<Fmck_e[11] <<  endl;
		// if (abs(Fmck_e[11] - 5.80215) > 0.00001) continue;
		// else cout << "This is ev " << entry << " " <<Fmck_e[11] <<  endl;

		//cout << "entry: " << entry << "; Eventnr: " << Eventnr << "; Runnr_prev: " << Runnr_prev << endl;

		take_had_event = kFALSE;
		take_part_event = kFALSE;

		take_event = kTRUE;
		take_event_trig = kTRUE;
		here_is_prph = kFALSE;
		here_is_jet = kFALSE;

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
				}
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
				cout << "EM NOT PRESERVED" << endl;
				continue;
			}

			//part nohad || part_bin != had_bin
				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_et->FindBin(part_et) != hist.had_cross_et->FindBin(had_et)) ) 
					hist.part_nohad_cross_et->Fill(part_et);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_eta->FindBin(part_eta) != hist.had_cross_eta->FindBin(had_eta)) ) 
					hist.part_nohad_cross_eta->Fill(part_eta);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_Q2->FindBin(part_Q2) != hist.had_cross_Q2->FindBin(had_Q2)) ) 
					hist.part_nohad_cross_Q2->Fill(part_Q2);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_x->FindBin(part_x) != hist.had_cross_x->FindBin(had_x)) ) 
					hist.part_nohad_cross_x->Fill(part_x);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_et_jet->FindBin(part_et_jet) != hist.had_cross_et_jet->FindBin(had_et_jet)) ) 
				{
					hist.part_nohad_cross_et_jet->Fill(part_et_jet);
					hist.part_nohad_cross_et_jet2->Fill(part_et_jet);
				}

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_eta_jet->FindBin(part_eta_jet) != hist.had_cross_eta_jet->FindBin(had_eta_jet)) ) 
					hist.part_nohad_cross_eta_jet->Fill(part_eta_jet);
				//new vars
				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_xgamma->FindBin(part_xgamma) != hist.had_cross_xgamma->FindBin(had_xgamma)) ) 
					hist.part_nohad_cross_xgamma->Fill(part_xgamma);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_xp->FindBin(part_xp) != hist.had_cross_xp->FindBin(had_xp)) ) 
					hist.part_nohad_cross_xp->Fill(part_xp);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_dphi->FindBin(part_dphi) != hist.had_cross_dphi->FindBin(had_dphi)) ) 
					hist.part_nohad_cross_dphi->Fill(part_dphi);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_deta->FindBin(part_deta) != hist.had_cross_deta->FindBin(had_deta)) ) 
					hist.part_nohad_cross_deta->Fill(part_deta);

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_dphi_e_ph->FindBin(part_dphi_e_ph) != hist.had_cross_dphi_e_ph->FindBin(had_dphi_e_ph)) ) 
				{
					hist.part_nohad_cross_dphi_e_ph->Fill(part_dphi_e_ph);
				}

				if ( (take_part_event && !take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_deta_e_ph->FindBin(part_deta_e_ph) != hist.had_cross_deta_e_ph->FindBin(had_deta_e_ph)) ) 
					hist.part_nohad_cross_deta_e_ph->Fill(part_deta_e_ph);

			//had nopart || part_bin != had_bin
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_et->FindBin(part_et) != hist.had_cross_et->FindBin(had_et)) ) 
					hist.had_nopart_cross_et->Fill(had_et);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_eta->FindBin(part_eta) != hist.had_cross_eta->FindBin(had_eta)) ) 
					hist.had_nopart_cross_eta->Fill(had_eta);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_Q2->FindBin(part_Q2) != hist.had_cross_Q2->FindBin(had_Q2)) ) 
					hist.had_nopart_cross_Q2->Fill(had_Q2);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_x->FindBin(part_x) != hist.had_cross_x->FindBin(had_x)) ) 
					hist.had_nopart_cross_x->Fill(had_x);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_et_jet->FindBin(part_et_jet) != hist.had_cross_et_jet->FindBin(had_et_jet)) ) 
				{
					hist.had_nopart_cross_et_jet->Fill(had_et_jet);
					hist.had_nopart_cross_et_jet2->Fill(had_et_jet);
				}
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_eta_jet->FindBin(part_eta_jet) != hist.had_cross_eta_jet->FindBin(had_eta_jet)) ) 
					hist.had_nopart_cross_eta_jet->Fill(had_eta_jet);
				//new vars

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_xgamma->FindBin(part_xgamma) != hist.had_cross_xgamma->FindBin(had_xgamma)) ) 
					hist.had_nopart_cross_xgamma->Fill(had_xgamma);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_xp->FindBin(part_xp) != hist.had_cross_xp->FindBin(had_xp)) ) 
					hist.had_nopart_cross_xp->Fill(had_xp);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_dphi->FindBin(part_dphi) != hist.had_cross_dphi->FindBin(had_dphi)) ) 
					hist.had_nopart_cross_dphi->Fill(had_dphi);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_deta->FindBin(part_deta) != hist.had_cross_deta->FindBin(had_deta)) ) 
					hist.had_nopart_cross_deta->Fill(had_deta);

				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_dphi_e_ph->FindBin(part_dphi_e_ph) != hist.had_cross_dphi_e_ph->FindBin(had_dphi_e_ph)) ) 
				{
					hist.had_nopart_cross_dphi_e_ph->Fill(had_dphi_e_ph);
				}
				if ( (!take_part_event && take_had_event) || 
					 (take_part_event && take_had_event && 
					 	hist.part_cross_deta_e_ph->FindBin(part_deta_e_ph) != hist.had_cross_deta_e_ph->FindBin(had_deta_e_ph)) ) 
					hist.had_nopart_cross_deta_e_ph->Fill(had_deta_e_ph);
		
		hadron_level_jets.clear();
		
	}// for entry over entries
	ofs.close();
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
